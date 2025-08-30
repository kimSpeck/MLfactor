# all simulation parameters are in this list
# there are sublists for ...
#   ... controlling the simulation design in general (params$sim)
#   ... the data generation (params$dgp) 
#   ... fitting the data (params$fit)
params <- list()

################################################################################
# utils 
################################################################################
findFLs <- function(alpha, nItems, varT, varItems = 1) {
  # find factor loadings for an essential tau-equivalent model for ...
  # ... fixed Cronbach's Alpha
  # ... fixed number of items
  # ... and identical variances for all items in the model (default value = 1)
  # we also assume uncorrelated scores, no cross-loadings, and no covariances of residuals
  
  if (length(varItems) != 1) {
    stop("This function is only implemented for equal variances across items!")
  } 
  
  return(sqrt(alpha/((nItems - (nItems - 1) * alpha) * varT)))
} 

# calculate the error variance of the item with MI violation
getErrorVar <- function(loadingVec, varianceVec, targetVar) {
  ###
  # calculate the error variance for items in which loadings on the true score changes
  ## input:
  # loadingVec  - [vector] vector of factor loadings on all factors that affect the item
  # varianceVec - [vector] vector of factor variances of all factors that affect the item
  # targetVar   - [scalar] target variance of the item 
  ## output:
  #             - [scalar] amount of error variance to arrive at target variance of the item
  ###
  
  targetVar - sum(loadingVec**2 * varianceVec)
}

# calculate the coefficient implement scalar MI violations of specific variance
getB0 <- function(flT, varT, varItems, varE, varG) {
  ## input:
  # flT       - [scalar] factor loading on the common factor
  # varT      - [scalar] variance of the common factor
  # varItems  - [scalar] variance of the manifest variable
  # varE      - [scalar] variance of the residual in the item
  # varG      - [scalar] variance of the binary and balanced grouping variable
  ## output:
  #           - [scalar] coefficient that indicates mean difference in manifest variable between groups
  
  sqrt((varItems - flT**2 * varT - varE) / varG)
}

getB1 <- function(flT, varT, varItems, varE, varG, meanG) {
  ## input:
  # flT       - [scalar] factor loading on the common factor
  # varT      - [scalar] variance of the common factor
  # varItems  - [scalar] variance of the manifest variable
  # varE      - [scalar] variance of the residual in the item
  # varG      - [scalar] variance of the binary and balanced grouping variable
  # meanG     - [scalar] expected values of the binary and balanced grouping variable
  ## output:
  #           - [scalar] coefficient that indicates difference in factor loading between groups

  a = (meanG**2 * varT + varG * varT)
  b = (flT * varT)
  c = (-1 + flT**2 * varT + varE)
  
  est1 <- (-b - sqrt(b**2 - 4 * a * c)) / (2 * a)
  est2 <- (-b + sqrt(b**2 - 4 * a * c)) / (2 * a)
  
  return(max(c(est1, est2)))
}

getBy <- function(R2, sigmaEy, varT, cAlpha, nFactors = 1) {
  ## input:
  # flT       - [scalar] factor loading on the common factor
  # varT      - [scalar] variance of the common factor
  # varItems  - [scalar] variance of the manifest variable
  # varE      - [scalar] variance of the residual in the item
  # varG      - [scalar] variance of the binary and balanced grouping variable
  # meanG     - [scalar] expected values of the binary and balanced grouping variable
  ## output:
  #           - [scalar] coefficient that indicates difference in factor loading between groups
  
  return(sqrt( (R2 * sigmaEy) / (varT * (1 - R2))))
}



################################################################################
# simulation design 
################################################################################
params$sim$nTrain <- 1000 # number of training samples
params$sim$nTest <- 1 # number of holdout test samples to evaluate models on 
params$sim$nSamples <- params$sim$nTrain + params$sim$nTest

# !if 50% of training sample size is used as test sample, test sample sizes vary across N observation conditions
#   as a result, the empirical SE of Rsquared changes according to N 
#   -> fix test sample size 
params$sim$testN <- 1000            # fixed test sample size across all N conditions 

################################################################################
# data generating process
################################################################################
##### measurement model #####
#this is training sample size; test/holdout sample is set to 1000 inside function

# sample size 
#   - McClure et al. c(50, 100, 250, 500, 1000)
#   - Pratiwi et al. c(100, 300, 500)
params$dgp$trainN <- c(50, 100, 250, 500, 1000) # sample size

# number of items per component/score/factor
#   - McClure et al. c(9, 21): eher heterogene Items in clinical scales
#   - Pratiwi et al. c(4, 8): eher homogene Items in personality inventories
#   - personality inventories in the literature: 
#     + IPIP-NEO-120 (Johnson, 2014): 4 items per facet.  
#       (30 facets × 4 items each = 120 total items)  
#     + NEO-PI-R (Costa & McCrae, 1992): 8 items per facet.  
#       (5 domains × 6 facets per domain × 8 items = 240 total items)  
#     + BFI-2 (Soto & John, 2017): 4 items per facet.  
#       (5 domains × 3 facets per domain × 4 items = 60 total items)
params$dgp$P <- c(4, 8) 

# number of components/factors/scores
params$dgp$Q <- c(1) # maybe more later! but uncorrelated in any way!

# outcome parameters
params$dgp$nOutcome <- 1 # number of outcome variables
params$dgp$mOut <- 0  # mean value of the outcome (without manipulation)
params$dgp$sigmaEy <- 1 # irreducible error; measurement error in the outcome

# item variance is fixed at 1 (constant item variance; 
#   avoid alternative explanation mechanisms of MI violations e.g., increase in 
#   predictor or component variance)
params$dgp$varItems <- 1
# params$dgp$mS <- 0
params$dgp$mItems <- 0 # item means

# variance decomposition (variances that constitute the Item)
params$dgp$varT <- c(0.5, 0.25) # variance of the true factor (c(0.5, 0.33, 0.25))
params$dgp$varStable <- c(0.25, 0.5) # variance of the specific item effect  (c(0.5, 0.33, 0.25))
params$dgp$flS <- 1 # factor loadings of the specific item variance (if specific item variance exists)
params$dgp$flStable <- 1 # specific item variance

# Cronbach's alpha of the components/factors scores
# see Soto & John (2017), Maples et al. (2014), Johnson (2014)
# however, simulating reliabilities which are really high means that there is not
#   enough variance left in the item to distribute on specific item variance and 
#   residual variance! Thus, we arrive at negative residual variances which cannot
#   exist in simulated or real data!
# params$dgp$cAlpha <- c(0.6, 0.75, 0.9) 
# params$dgp$cAlpha <- c(0.55, 0.7, 0.85) 
# params$dgp$cAlpha <- c(0.6, 0.7, 0.8) # oder c(0.6, 0.8)
params$dgp$cAlpha <- c(0.6, 0.8)

# for c(0.6, 0.7, 0.8) we get two interesting boundary conditions
#   - one condition that closely matches McClure et al. (2024):
#     4 items, alpha = 0.8, varT = 0.50, varS = 0.25, flT = 1.0000000, flS = 1, varE = 0.2500000
#   - one condition in which we only have true factor and specific item variance
#     4 items, alpha = 0.8, varT = 0.25, varS = 0.5, flT = 1.4142136, flS = 1, varE = 0.0000000

# measurement invariance violations
# ... maybe add resdiual or stick to the technical argument (see Millsap, 1997)
# ... tested configural MI violations (flT = 0) in piloting and coding but not as interesting(?)
params$dgp$MIvio <- c("specific", "scalar", "metric") 

# # MI weight does not exist anymore since MI weights are controlled by beta0 and beta1
# configural MI by using 0 as a weight
# metric MI by using -0.8 as a weight
# scalar MI by using unit weight (= identical weights) but instead change mean of specific item effect
# params$dgp$miWeight = c(0, -0.8, 1) # weight of measurement variant items


# manipulation measurement invariance violations
params$dgp$nGroups = 2  # number of groups (MI due to categorical variable)
params$dgp$varG = 0.25 # variance of a binary grouping variable with p = 0.5
params$dgp$meanG = 0.5 # expected value of a binary grouping variable with p = 0.5

##### structure model #####
params$dgp$R2 <- c(0.3, 0.8) # how much variance should be explained in the outcome
# params$dgp$R2 <- c(0.2, 0.4, 0.6) 
# correlation between the grouping variable and the outcome
params$dgp$rGroupOut <- c(0, 0.2, 0.5)

# we assume ...
# ... uncorrelated factors and uncorrelated item residuals 
# ... an essential tau-equivalent model
#   this means that factor loadings of every item are identical
#   and that the height of the factor loadings control how much variance is explained
#     by the factor vs. unique item and error variance 
params$dgp$grid <- expand.grid(P = params$dgp$P, # nItems per factor
                               Q = params$dgp$Q, # nFactors
                               cAlpha = params$dgp$cAlpha, # Cronbach's Alpha
                               varT = params$dgp$varT, # true factor variance
                               # miItems = seq_len(3), # number of MI violating items
                               miItems = 0:2, # number of MI violating items
                               corGO = params$dgp$rGroupOut, # correlation group -- outcome
                               R2 = params$dgp$R2) # explained variance outcome

# number of items with deviation in factor loadings
# here: include miItems in grid; but for now only simulate one MI violating item
# params$dgp$miItems = c(1, 2, 3)  # number of items with deviation in factor loadings
# params$dgp$miItems = c(1, P/2, (3/4)*P)
# params$dgp$grid$miItems = ifelse(params$dgp$grid$miItems == 1, 1, 
#                                  ifelse(params$dgp$grid$miItems == 2, 
#                                         params$dgp$grid$P/2, 3*(params$dgp$grid$P/4)))  # number of items with deviation in factor loadings
# 0 items to model no invariance violation or no "effect" in the respective condition
params$dgp$grid$miItems = ifelse(params$dgp$grid$miItems == 2, params$dgp$grid$P/2, params$dgp$grid$miItems) 

# here: differentiate between var{specific, metric, scalar} and actually varS 
params$dgp$grid$varStable <- ifelse(params$dgp$grid$varT == 0.5, 0.25, 0.5)
params$dgp$tauFL <- mapply(FUN = findFLs, alpha = params$dgp$grid$cAlpha, nItems = params$dgp$grid$P, varT = params$dgp$grid$varT)

params$dgp$grid$flT <- params$dgp$tauFL
params$dgp$grid$flS <- 1

# assuming every item is composed of true factor variance, specific item variance and error variance 
params$dgp$varE <- sapply(seq_len(dim(params$dgp$grid)[1]), function(iRow) {
  getErrorVar(loadingVec = c(params$dgp$grid$flT[iRow], params$dgp$grid$flS[iRow]),
              varianceVec = c(params$dgp$grid$varT[iRow], params$dgp$grid$varStable[iRow]),
              targetVar = params$dgp$varItems)                        
}) 

params$dgp$grid$varE <- params$dgp$varE

# complete grid with coefficient values
# add all MI violations to the grid including the respective parameters for the matrices
params$dgp$grid <- params$dgp$grid[rep(seq_len(nrow(params$dgp$grid)), 
                                       times = length(params$dgp$MIvio)), ]
params$dgp$grid$MIvio <- rep(params$dgp$MIvio, 
                                  each = dim(params$dgp$grid)[1]/length(params$dgp$MIvio)) 

# get coefficients for grouping variable on loadings or on item mean
#   ... keep the variance of the specific item effect, the mean difference and the 
#       change in factor loadings on the factor constant! 
#   ... do I need to standardize something?! test it! 
#   ... we assume two groups, which means that the grouping variable has var(G) = 0.25
params$dgp$grid$b0 <- sapply(seq_len(dim(params$dgp$grid)[1]), function(iRow) {
  if (params$dgp$grid$MIvio[iRow] == "scalar") {
    return(getB0(flT = params$dgp$grid$flT[iRow], varT = params$dgp$grid$varT[iRow], 
                 varItems = params$dgp$varItems, varE = params$dgp$grid$varE[iRow], 
                 varG = params$dgp$varG))
  } else {
    return(0)
  }
})

params$dgp$grid$b1 <- sapply(seq_len(dim(params$dgp$grid)[1]), function(iRow) {
  if (params$dgp$grid$MIvio[iRow] == "metric") {
    return(getB1(flT = params$dgp$grid$flT[iRow], varT = params$dgp$grid$varT[iRow], 
                 varItems = params$dgp$varItems, varE = params$dgp$grid$varE[iRow], 
                 varG = params$dgp$varG, meanG = params$dgp$meanG))
  } else {
    return(0)
  }
})

params$dgp$grid$betaT <- sapply(seq_len(dim(params$dgp$grid)[1]), function(iRow) {
  return(getBy(R2 = params$dgp$grid$R2[iRow], sigmaEy = params$dgp$sigmaEy,
               varT = params$dgp$grid$varT[iRow], cAlpha = params$dgp$grid$cAlpha[iRow]))
})
  
# adjust parameters to MIvio to correctly fill in matrices
# specific item variance only gets factor loadings if specific item variance comes from specific
params$dgp$grid$flS <- ifelse(params$dgp$grid$MIvio != "specific", 0, params$dgp$flS)

# varS is only varStable for MIvio == "specific" otherwise it is zero! 
params$dgp$grid$varS <- ifelse(params$dgp$grid$MIvio != "specific", 0, params$dgp$grid$varStable)

# if miItems == 0, remove all weights that could possibly apply any manipulation
# strictly speaking, this creates duplicate conditions. 
# However, it is a sanity check for confounding conditions in different MIvio processes
# ... for MIvio == "specific": betaS = 0
# ... for MIvio == "scalar": b0 = 0
# ... for MIvio == "metric": b1 = 0
params$dgp$grid[which(params$dgp$grid$miItems == 0 & 
                        params$dgp$grid$MIvio == "specific"), "betaS"] <- 0

params$dgp$grid[which(params$dgp$grid$miItems == 0 & 
                        params$dgp$grid$MIvio == "scalar"), "b0"] <- params$dgp$mItems 

params$dgp$grid[which(params$dgp$grid$miItems == 0 & 
                        params$dgp$grid$MIvio == "metric"), "b1"] <- 0

# # remove conditions with multiple MI violating items if MIvio == "none"
# idxRmConds <- which(params$dgp$grid$MIvio == "none" & params$dgp$grid$miItems > 1)
# params$dgp$grid <- params$dgp$grid[-idxRmConds, ]

# calculate regression coefficients of the true factor(s) accordingly
#   the specific item variance adds to the fixed R² (based on the grid)
# presumably, betaT changes depending on the variance in the true factor and Cronbach's Alpha
#   and hence, the reliability of the predictor
# here: Herleitung der Regressionsgewichte für den true factor! 
params$dgp$betaT <- 0.5 

# specicic item effects do not have an effect in every MIvio manipulation except specific!
params$dgp$betaS <- 1
params$dgp$grid$betaS <- ifelse(params$dgp$grid$MIvio != "specific", 0, params$dgp$betaS)

# # I don't need the nGroups parameter since I always simulate two groups but with identical parameter values 
#   if the true DGP has only one group
# params$dgp$varFLgrid$nGroups <- ifelse(params$dgp$varFLgrid$MIvio == "specific", 1, params$dgp$nGroups)

################################################################################
# model fitting
################################################################################
##### Hyperparameter Tuning #####
params$fit$nfolds <- 10

##### ENET #####
# fixed set of lambda parameters (only used if warmStart = FALSE)
params$fit$lambda <- 10^seq(-1, 1, length = 100) 
# use warmStart?
#  pro: this is how caret operates (i.e., get lambda values for alpha = 0.5)
#         -> thus, lots of users use warm start for lambda tuning
#       with warmStart exact model is chosen more frequently
#  con: less control over lambda grid to be searched      
params$fit$warmStart <- TRUE  # thus, lambda parameter are not used
params$fit$alpha <- seq(0, 1, length.out = 20) # alpha with at least 20 steps

# the code can handle setParam$fit$lambdaCrit = {"1se", "min"} or every subset
# "lambda.min" = the lambda at which the smallest MSE is achieved.
# "lambda.1se" = the largest lambda at which the MSE is within one SE of the smallest MSE (default).
# here: "one-standard-error" rule for choosing lambda (Hastie et al. 2009)
#   Friedman et al. 2010. Regularization Paths for Generalized Linear Models via Coordinate Descent.
params$fit$lambdaCrit <- "1se" # 1se
