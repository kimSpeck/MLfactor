# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R") # general functions for simulation management
source("utils/simTools.R") # general functions for data simulation

# functions to generate data
source("utils/sampleData.R")

require(mvtnorm)
library(magic)
library(lavaan)

# test it
condVec <- which(params$dgp$grid$P == 4 &
                params$dgp$grid$cAlpha == 0.8 &
                params$dgp$grid$miItems == 2 &
                  params$dgp$grid$corGO == 0 &
                  params$dgp$grid$varE == 0 &
                params$dgp$grid$MIvio == "metric")
params$dgp$grid[condVec,]
iCond <- condVec[length(condVec)]
iCond <- condVec[1]
params$dgp$grid[iCond,]

##### measurement model #####
# ... always simulate nGroups = 2 in every condition but keep parameters identical! 

# generate lambda matrix
lambdaSEM <- defineLambda(flT = params$dgp$grid$flT[iCond], flS = params$dgp$grid$flS[iCond], 
                          nGroups = params$dgp$nGroups, nItems = params$dgp$grid$P[iCond], 
                          miItems = params$dgp$grid$miItems[iCond], 
                          # for MIvio %in% c("specific", "scalar") miWeight is automatically flT
                          miWeight = params$dgp$grid$flT[iCond] + params$dgp$grid$b1[iCond], 
                          nOutcome = params$dgp$nOutcome, nFactors = params$dgp$grid$Q[iCond])

lambdaSEM

nLatentSEM <- ncol(lambdaSEM)

# generate covariance matrix for the factors
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
# here! adapt function to model multiple group variables to generate correlation
#       between the grouping variable and the outcome
MIflag <- if (params$dgp$grid$MIvio[iCond] != "specific") T else F
covSEM <- definePhi(nLatentSEM, 
                    varT = params$dgp$grid$varT[iCond], 
                    varStable = params$dgp$grid$varStable[iCond],
                    varS = params$dgp$grid$varS[iCond],
                    nItems = params$dgp$grid$P[iCond], MIvio = MIflag,
                    miItems = params$dgp$grid$miItems[iCond], #here: this will become part of the grid!
                    # miItems = 0, # check warning message (;
                    varNames = colnames(lambdaSEM))
covSEM

# generate residual variance matrix
# I remove the specific item variance and instead model the MI violation with exactly
#   the same variance as for the specific item variance
#   thus, I do not need to adjust the residual variance of the MI violating item, right?
thetaSEM <- defineTheta(varE = params$dgp$grid$varE[iCond], 
                        nItems = params$dgp$grid$P[iCond], 
                        nGroups = params$dgp$nGroups, nOutcome = params$dgp$nOutcome)
thetaSEM
##### structure model #####
# generate beta matrix
betaSEM <- defineBeta(nLatentSEM, 
                      betaT = params$dgp$grid$betaT[iCond], 
                      betaS = params$dgp$grid$betaS[iCond], 
                      miItems = params$dgp$grid$miItems[iCond], 
                      nItems = params$dgp$grid$P[iCond], 
                      varNames = colnames(lambdaSEM), nFactors = 1)
betaSEM
##### complete SEM #####
# calculate model implied covariance matrix for ...
#     simulating from SEM by sampling from multivariate normal distribution 
sigmaSEM <- getSigma(lambdaSEM, betaSEM, covSEM, thetaSEM)

# generate mean vector for data sampling
# simulate group-outcome correlations by mean differences in the outcome
#   -> calculate difference in outcome for target correlation assuming balanced design with 2 groups! 
if (params$dgp$grid$corGO[iCond] != 0) {
  Y1g1Idx <- which(colnames(sigmaSEM) == "Y1_g1")
  Y1g2Idx <- which(colnames(sigmaSEM) == "Y1_g2")
  varWithinY <- 0.5 * sigmaSEM[Y1g1Idx, Y1g1Idx] + 0.5 * sigmaSEM[Y1g2Idx, Y1g2Idx]
  meanYg2 <- getMeanYg2(targetRho = params$dgp$grid$corGO[iCond], params$dgp$varG, 
                        varWithinY)
} else {
  # depending on correlation between group and outcome; if cor = 0: mYg2 = mOut
  meanYg2 <- params$dgp$mOut
}

M <- defineMeanVec(mItems = params$dgp$mItems, mOut = params$dgp$mOut, 
                   nItems = params$dgp$grid$P[iCond], 
                   nGroups = params$dgp$nGroups,
                   nOutcome = params$dgp$nOutcome, 
                   mXg2 = params$dgp$grid$b0[iCond], 
                   mYg2 = meanYg2,
                   miItems = params$dgp$grid$miItems[iCond])


set.seed(2892) # here! manage seeds later on!
# dataWide <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, 
#                                        mean = M, sigma = sigmaSEM)
tryN <- 20000
dataWide <- mvtnorm::rmvnorm((tryN)/params$dgp$nGroups, 
                             mean = M, sigma = sigmaSEM)
# name columns of the data
colnames(dataWide) <- c(paste0(rep(c(sapply(seq(params$dgp$grid$P[iCond]), 
                                            function(iP) {paste0("X", iP)})), times = params$dgp$nGroups), 
                               "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond])), 
                        "Y_g1", "Y_g2")

# here: calculate sum score to put it in the model as a predictor?!
# wide to long
data <- data.frame(dataWide)
data <- tidyr::pivot_longer(data, cols = !c(Y_g1, Y_g2), names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
data$Y <- ifelse(data$group == "g1", data$Y_g1, data$Y_g2)

# check variances of the manifest variables
data$group <- ifelse(data$group == "g1", 0, 1)
(varVector <- apply(data, 2, var))

# check correlations between group and outcome
cor(data["Y"], data["group"])

# check SEM data
m0 <- '# measurement model
       # T =~ NA * X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8
       T =~ NA * X1 + X2 + X3 + X4

       # all latent representation of manifest outcome 
       #    loading = 1, residual variance of manifest variable = 0
       fY =~ 1 * Y
       Y ~~ 0*Y
       
       T ~~ 0.5*T
       #T ~~ 0.25*T
       
       # structure model
       fY ~ T'

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
#fit_singleGroup@Model@GLIST
summary(fit_singleGroup, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)


# # for the multiple group analysis, the data looks fine in both groups (see first test above)
# fit_multipleGroups <- lavaan::sem(m0, data = data[,-c(1, 2)], group = "group")
# # fit_multipleGroups@Model@GLIST
# summary(fit_multipleGroups, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

# check Cronbach's Alpha
params$dgp$grid[iCond, ]
vars <- sapply(seq_len(params$dgp$grid$P[iCond]), function(x) paste0("X", x))
(CronbachsAlpha <- getCronbachsAlpha(data[,vars]))

semTools::reliability(fit_singleGroup)["alpha",]

