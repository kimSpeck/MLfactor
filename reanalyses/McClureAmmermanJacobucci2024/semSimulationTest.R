################################################################################
# simulate data from SEM (and get identical results?!)
################################################################################
# simulate from SEM:
#   - specific effects as latent variables with only 1 item ("all latent")
#   - simulation does not assume standardized data, thus, the variance of y/outcome != 1
#     variances of every other variable is known in the simulation
#     if we simulate only factor structure without the outcome however, we cannot 
#       simulate path coefficients from the specific item effect on the outcome! 
#     -> therefore we need to simulate the outcome within the SEM model to be able 
#         to simulate effects of the specific item effect on the outcome!
#     !!! in the cov Matrix of the SEM there are the residual variance of y as endognuous 
#         variable; thus variance in the matrix is residual variance (params$sigmaErrY)

# variance in Y is way smaller in SEM simulation compared to Jacobucci simulation
#   betaS.cond == 1: var(Y) = 12.31 vs. var(data$Y) = 2.16154
#   betaS.cond == 4: var(Y) = 11.80889 vs. var(data$Y) = 1.178903
#   ->  this is due to the rather "formative" conception of the true score! 
#       the true score enters y multiple times in the original Jacobucci simulation
#   -> now that we simulate a reflective model in the original simulation the variance
#      matches with the sem formulation of the data 

# step by step approach
# 1. simulate specific item effect without measurement invariance violations
# 2. simulate from multiple group SEM without any MI violation
# 3. simulate from multiple group SEM with configural MI violation
# 4. simulate from multiple group SEM with metric measurement invariance violation
# 5. simulate from multiple group SEM with scalar measurement invariance violation
# 6. add correlation between MI violating group variable and the outcome
#     -> start simulating Y for each group
# 7. 

# to do: 
# test multiple group SEM with grouping variable for all MI violation conditions
# test MI violations for multiple Items
# write functions for data generating model 

################################################################################
# load packages
################################################################################
setwd("reanalyses/McClureAmmermanJacobucci2024")
source('simTools.R')
source('setParameters.R')
source('original_Linear-Regression-Bias-Variance_functions.R')

source("../../simTools.R")
# source("../../setParams.R")

require(lavaan)
require(tidyr)
require(mvtnorm)
library(magic)

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# test it
N.train <- 1000
J <- 9
betaT <- 0.5
seed <- params$subGrid$seed[1]
measVar.cond <- 1
betaS.cond <- 1
beta_S <- setBetaS(betaS.cond, params$sigma[measVar.cond,"sigmaS"], J, params)
betaS <- params$betaS[betaS.cond]/params$sigma[measVar.cond,"sigmaS"]
params$nOutcome <- 1

N = params$holdout + N.train
W = popWeight(params$M, J) #creates a weight matrix that partitions J in M~equal unit-weighted composites

# set variances 
if (params$M == 1) { # one component case = one value
  Sigma_T <- params$sigma[measVar.cond,"sigmaT"]
} else { # multiple component case = matrix with values
  Sigma_T = matrix(params$sigma[measVar.cond,"sigmaT"], 
                   nrow=params$M , ncol=params$M)
  diag(Sigma_T) <- params$sigma[measVar.cond,"sigmaT"]
}

Sigma_S = diag(params$sigma[measVar.cond,"sigmaS"],nrow = J,ncol = J)
Sigma_e = diag(params$sigma[measVar.cond,"sigmaE"],nrow = J,ncol = J)

##### simulate only measurement model #####
# factor loading matrix without outcome; factors in columns, items in rows
# -> to simulate specific item effects, the outcome has to be simulated within the
#     SEM context
# lambda <- cbind(trueScore = rep(1, J),
#                 specificEffects = t(diag(J)))
# nFactors <- ncol(lambda)
# covMat <- diag(nFactors) * c(Sigma_T, diag(Sigma_S))
# theta <- Sigma_e
# # for simulating only from measurement model 
# sigmaMat <- lambda %*% covMat %*% t(lambda) + theta

################################################################################
# simulate from SEM 
# (with lots of explanations)
################################################################################
# factor loading matrix with outcome; but variance of the outcome is unknwon
#   therefore: factor loading matrix without outcome below

##### simulate specific item effect with only 1 item #####
# every item loads on the true score and on the specific effect; but specific effect
#   is not identified (only 1 item); therefore we cannot check if the variance is correct

lambdaSEM <- defineLambda(1, J, miItems = 0, miWeight = NA, params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)

# lambdaSEM <- cbind(trueScore = c(rep(1, J), rep(0, params$dgp$nOutcome)),
#                    specificEffects = t(cbind(diag(J), 0)),
#                    Y = c(rep(0, J), rep(1, params$dgp$nOutcome)))
# colnames(lambdaSEM) <- c("T", 
#                          sapply(seq(J), function(iJ) {paste0("S", iJ)}),
#                          "Y")


# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
# factors are not correlated 
#   T, S and e are sampled from different multivariate normal distributions, thus,
#   factors need to be independent/uncorrelated
#   relating the model to substantive research: nuance information is independent 
#   from facet information (this is the dgm, i.e., extra information in nuances 
#   that has nothing to do with what is already explained by facet)
# variables from structure model need to be included in factor correlation matrix 
# for the SEM simulation 
covSEM <- diag(nFactors) * c(Sigma_T, diag(Sigma_S), params$sigmaErrY)
covSEM

# path coefficient matrix (effect of factors on the target/criterion)
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- matrix(rep(0, 
                      (J + params$M + params$dgp$nOutcome) + # diagonal
                        ((J + params$M + params$dgp$nOutcome)*(J + params$M + params$dgp$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(betaT, beta_S, 0)

# error - Kovarianzmatrix der Residualvariablen mit p x p (= Anzahl der manifesten Variablen)
# outcome as manifest variable but all latent! (therefore the error variance is zero?)
thetaSEM <- diag(c(diag(Sigma_e), rep(0, params$nOutcome))) # without helper items

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
M <- c(rep(params$mS, J), rep(0, params$nOutcome)) # without helper items

# calculate model implied covariance matrix
# for simulating from SEM
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
set.seed(seed)
data <- mvtnorm::rmvnorm(N.train + params$holdout, mean = M, sigma = sigmaSEM)
colnames(data) <- c(paste0("X", seq_len(J)), "Y") # without helper items

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) 
# 0.228212 0.2486797 0.2208203 0.2284228 0.2128141 0.2395708 0.2245688 0.2379814 0.2142187

cor(data[,"Y"], rowMeans(data[,-outcomeIdx]))
# 0.3

######
# check SEM data
######
m0 <- '# measurement model
       T =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9

       # all latent representation of manifest outcome 
       #    loading = 1, residual variance of manifest variable = 0
       fY =~ 1 * Y
       Y ~~ 0*Y
       
       # structure model
       fY ~ T'

# without specific item effect on outcome
apply(data, 2, var)
fit <- lavaan::sem(m0, data = data)
fit@Model@GLIST

# lambda: all factor loadings in lambda are close to 1 as simulated;
# theta: error variances are 0.5 for every item with 0.25 due to specific effect and error respectively
# psi: the variance of the true score is .5; the (residual) variance of the outcome is 1
# beta: the estimated coefficient for true score on the outcome is 0.5

# # if I would simulate data with effect of specific item and named it dataSpecific
# apply(dataSpecific, 2, var)
# fit <- lavaan::sem(m0, data = dataSpecific)
# fit@Model@GLIST

##### simulate from multiple group SEM without any MI violation ##### 
# block-diagonal binding of matrices 
nGroups <- 2
lambda_g1 <- cbind(TS = rep(1, J),
                   specificEffects = diag(J))
lambda_g2 <- cbind(TS = rep(1, J),
                   specificEffects = diag(J))
# without group
lambdaSEM <- rbind(cbind(magic::adiag(lambda_g1, lambda_g2), 
                         rep(0, params$nOutcome)),
                   c(rep(0, (J + 1)*nGroups), rep(1, params$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(J), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J+1)), "Y")

# ! defineLambda function cannot model only one outcome despite multiple groups

# # with group 
# lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
#                          # group = rep(c(0, 1), each = J), 
#                          group = rep(c(0, 0), each = J), 
#                          rep(0, params$nOutcome)),
#                    c(rep(c(0, 1), each = (J + 1)), 1, rep(0, params$nOutcome)),
#                    c(rep(0, (J + 1)*nGroups + 1), rep(1, params$nOutcome)))
# colnames(lambdaSEM) <- c(paste0(rep(c("T", 
#                                       sapply(seq(J), 
#                                              function(iJ) {paste0("S", iJ)})), times = 2), 
#                                 "_g", rep(seq_len(nGroups), each = J+1)),
#                          "group", "Y")

nFactors <- ncol(lambdaSEM)

# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
# without group
covSEM <- diag(nFactors) * c(
  rep(c(Sigma_T, diag(Sigma_S)), times = 2),  
  params$sigmaErrY) # residual variance!

# # with group
# covSEM <- diag(nFactors) * c(
#   rep(c(Sigma_T, diag(Sigma_S)), times = 2),  
#   0, # residual variance?
#   params$sigmaErrY)
covSEM

# path coefficient matrix (effect of factors on the target/criterion)
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
# without group
betaMat <- matrix(rep(0, 
                      ((J + params$M)*nGroups + params$nOutcome) + # diagonal
                        (((J + params$M)*nGroups + params$nOutcome)*((J + params$M)*nGroups + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0)

# # with group
# betaMat <- matrix(rep(0, 
#                       ((J + params$M)*nGroups + 1 + params$nOutcome) + # diagonal
#                         (((J + params$M)*nGroups + 1 + params$nOutcome)*((J + params$M)*nGroups + 1 + params$nOutcome - 1))), # off diagonal
#                   ncol = nFactors, byrow = T)
# betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0, 0)
betaMat

# error - Kovarianzmatrix der Residualvariablen mit p x p (= Anzahl der manifesten Variablen)
# outcome as manifest variable but all latent! (therefore the error variance is zero?)
# without group
thetaSEM <- diag(c(rep(diag(Sigma_e), nGroups), rep(0, params$nOutcome))) # without helper items

# # with group
# thetaSEM <- diag(c(rep(diag(Sigma_e), nGroups), 0, rep(0, params$nOutcome))) # without helper items
thetaSEM

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
M <- c(rep(params$mS, J*nGroups), rep(0, params$nOutcome)) # without helper items

# # with group
# M <- c(rep(params$mS, J*nGroups), 0, rep(0, params$nOutcome)) # without helper items

# calculate model implied covariance matrix
# for simulating from SEM
# simulating one outcome in multiple group scenario (block matrices) increases the
#   variance of the outcome in the model implied variance-covarince matrix from 1.125 to 1.25
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(J), 
                                        function(iJ) {paste0("X", iJ)})), times = 2), 
                           "_g", rep(seq_len(nGroups), each = J)), "Y")

# checḱ variances of all items
#   -> variances in all items (independent of group) are 1
apply(data_multipleGroup, 2, var)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data_multipleGroup) == "Y")
cor(data_multipleGroup[,"Y"], data_multipleGroup[,-outcomeIdx]) # ~ around 0.2 per item

cor(data_multipleGroup[,"Y"], rowMeans(data_multipleGroup[,-outcomeIdx]))
# 0.39


# wide to long
data <- data.frame(data_multipleGroup)
data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)

##### check data #####
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
groupIdx <- which(colnames(data) == "group")
cor(data$group, data[,-groupIdx]) # correlations with group  are all zero
#  Y = 0, X1 = 0.02289168, -0.0003050068, 0.0154861, -0.03704896, 0.01114973, 0.01600971, 0.002800397, 0.005275015, 0.003546012

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) # group is uncorrelated with outcome; other variables around ~0.2

groupIdx <- which(colnames(data) == "group")
cor(data[,"Y"], rowMeans(data[,-c(outcomeIdx, groupIdx)])) # 0.29

# for the multiple group analysis, the data looks fine in both groups (see first test above)
fit_multipleGroups <- lavaan::sem(m0, data = data, group = "group")
fit_multipleGroups@Model@GLIST

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
fit_singleGroup@Model@GLIST

##### simulate from multiple group SEM with configural measurement invariance violation #####
# when removing the loading of 1 item on the true score in one of the groups, this
#   item decreases in variance
# introduce MI violation in last item to avoid using it as marker item in SEM
nGroups <- 2
miItems <- 1
miWeight <- 0
lambda_g1 <- cbind(TS = rep(1, J),
                   specificEffects = diag(J))
lambda_g2 <- cbind(TS = rep(c(1, miWeight), c(J-miItems, miItems)),
                   specificEffects = diag(J))

lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                         rep(0, params$nOutcome)),
                   c(rep(0, (J + 1)*nGroups), rep(1, params$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(J), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J+1)), "Y")

nFactors <- ncol(lambdaSEM)

# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
covSEM <- diag(nFactors) * c(
  rep(c(Sigma_T, diag(Sigma_S)), times = 2),  
  params$sigmaErrY) # residual variance!

# path coefficients
betaMat <- matrix(rep(0, 
                      ((J + params$M)*nGroups + params$nOutcome) + # diagonal
                        (((J + params$M)*nGroups + params$nOutcome)*((J + params$M)*nGroups + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0)

thetaSEM <- diag(c(diag(Sigma_e),  # group 1
                   rep(c(Sigma_e[1, 1], Sigma_e[1, 1]+Sigma_T), c(J-miItems, miItems)), # group 2
                   rep(0, params$nOutcome))) # outcome

sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(J), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(nGroups), each = J)), "Y")

# checḱ variances of all items
#   -> variances in all items (independent of group) are 1
apply(data_multipleGroup, 2, var)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data_multipleGroup) == "Y")
cor(data_multipleGroup[,"Y"], data_multipleGroup[,-outcomeIdx]) 
# ~ around 0.2 per item; except X9_g2 which is obviously uncorrelated with the outcome

cor(data_multipleGroup[,"Y"], rowMeans(data_multipleGroup[,-outcomeIdx])) # 0.41

# wide to long
data <- data.frame(data_multipleGroup)
data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)

##### check data #####
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
groupIdx <- which(colnames(data) == "group")
cor(data$group, data[,-groupIdx]) # correlations with group  are all zero

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) # group is uncorrelated with outcome; 
                                    # other variables around ~0.2
                                    # X9 only around .1

groupIdx <- which(colnames(data) == "group")
cor(data[,"Y"], rowMeans(data[,-c(groupIdx, outcomeIdx)]))
# 0.298

# for the multiple group analysis, the data looks fine in both groups (see first test above)
fit_multipleGroups <- lavaan::sem(m0, data = data, group = "group")
fit_multipleGroups@Model@GLIST

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
fit_singleGroup@Model@GLIST

##### simulate from multiple group SEM with metric measurement invariance violation #####
nGroups <- 2
miItems <- 1
miWeight <- -0.8
lambda_g1 <- cbind(TS = rep(1, J),
                   specificEffects = diag(J))
lambda_g2 <- cbind(TS = rep(c(1, miWeight), c(J-miItems, miItems)),
                   specificEffects = diag(J))

lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                         rep(0, params$nOutcome)),
                   c(rep(0, (J + 1)*nGroups), rep(1, params$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(J), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J+1)), "Y")
rownames(lambdaSEM) <- c(paste0(rep(c(sapply(seq(J), 
                                             function(iJ) {paste0("X", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J)), "Y")

lambdaSEM

nFactors <- ncol(lambdaSEM)

# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
covSEM <- diag(nFactors) * c(
  rep(c(Sigma_T, diag(Sigma_S)), times = 2),  
  params$sigmaErrY) # residual variance!

# path coefficients
betaMat <- matrix(rep(0, 
                      ((J + params$M)*nGroups + params$nOutcome) + # diagonal
                        (((J + params$M)*nGroups + params$nOutcome)*((J + params$M)*nGroups + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0)


# error variance
itemVar <- 1
adjustedError <- getErrorVar(loadingVec = c(miWeight, W[1]),
                             varianceVec = c(Sigma_T, Sigma_S[1, 1]),
                             targetVar = itemVar)

thetaSEM <- diag(c(diag(Sigma_e),  # group 1
                   rep(c(Sigma_e[1, 1], adjustedError), c(J-miItems, miItems)), # group 2
                   rep(0, params$nOutcome))) # outcome

# introducing negative factor loadings changes correlation with all items of the factor
# variance is reduced, too
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
sigmaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(J), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(nGroups), each = J)), "Y")

# checḱ variances of all items
#   -> variances in all items (independent of group) are 1
apply(data_multipleGroup, 2, var)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data_multipleGroup) == "Y")
cor(data_multipleGroup[,"Y"], data_multipleGroup[,-outcomeIdx]) 
# ~ around 0.2 per item; except X9_g2 which correlates around -0.15

cor(data_multipleGroup[,"Y"], rowMeans(data_multipleGroup[,-outcomeIdx]))
# 0.45

# wide to long
data <- data.frame(data_multipleGroup)
data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)

##### check data #####
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
groupIdx <- which(colnames(data) == "group")
cor(data$group, data[,-groupIdx]) # correlations with group  are all zero

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) 
# group is uncorrelated with outcome; other variables around ~0.2; X9 only around 0

groupIdx <- which(colnames(data) == "group")
cor(data[,"Y"], rowMeans(data[,-c(groupIdx, outcomeIdx)]))
# 0.33

# for the multiple group analysis, the data looks fine in both groups (see first test above)
fit_multipleGroups <- lavaan::sem(m0, data = data, group = "group")
fit_multipleGroups@Model@GLIST

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
fit_singleGroup@Model@GLIST

##### simulate from multiple group SEM with scalar measurement invariance violation #####
nGroups <- 2
miItems <- 1
miWeight <- 1
lambda_g1 <- cbind(TS = rep(1, J),
                   specificEffects = diag(J))
lambda_g2 <- cbind(TS = rep(c(1, miWeight), c(J-miItems, miItems)),
                   specificEffects = diag(J))

lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                         rep(0, params$nOutcome)),
                   c(rep(0, (J + 1)*nGroups), rep(1, params$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(J), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J+1)), "Y")
rownames(lambdaSEM) <- c(paste0(rep(c(sapply(seq(J), 
                                             function(iJ) {paste0("X", iJ)})), times = 2), 
                                "_g", rep(seq_len(nGroups), each = J)), "Y")

lambdaSEM

nFactors <- ncol(lambdaSEM)

# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
covSEM <- diag(nFactors) * c(
  rep(c(Sigma_T, Sigma_S[1, 1], 0), c(1, J-miItems, miItems)), # group 1
  rep(c(Sigma_T, Sigma_S[1, 1], 0), c(1, J-miItems, miItems)), # group 2
  params$sigmaErrY) # residual variance!

# path coefficients
betaMat <- matrix(rep(0, 
                      ((J + params$M)*nGroups + params$nOutcome) + # diagonal
                        (((J + params$M)*nGroups + params$nOutcome)*((J + params$M)*nGroups + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0)

thetaSEM <- diag(c(diag(Sigma_e),  # group 1
                   # diag(Sigma_e), # group 2
                   rep(c(Sigma_e[1, 1], Sigma_e[1, 1]), c(J-miItems, miItems)), # group 2
                   rep(0, params$nOutcome))) # outcome

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
mSpecific <- 1
M <- c(rep(params$mS, J), # group 1 
       rep(c(params$mS, mSpecific), c(J-miItems, miItems)), # group 2 
       rep(0, params$nOutcome)) # without helper items

# introducing negative factor loadings changes correlation with all items of the factor
# variance is reduced, too
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
sigmaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(J), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(nGroups), each = J)), "Y")

# checḱ variances of all items
#   -> variance in item 9 (with scalar measurement invariance violation is only 0.75) in each group
#   -> true effect with variance of 0.5 and error with variance of 0.25; no variance of specific effect
#       specific effect is either 0 in g1 or 1 in g2 -> therefore only across groups var(x9) == 1
apply(data_multipleGroup, 2, var)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data_multipleGroup) == "Y")
cor(data_multipleGroup[,"Y"], data_multipleGroup[,-outcomeIdx]) # ~ around 0.25 per item

cor(data_multipleGroup[,"Y"], rowMeans(data_multipleGroup[,-outcomeIdx]))
# 0.46

# wide to long
data <- data.frame(data_multipleGroup)
data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)

##### check data #####
# without specific item effect on outcome
# check variances across groups now
# -> across group variances of all items are indeed 1
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
groupIdx <- which(colnames(data) == "group")
cor(data$group, data[,-groupIdx]) # due to scalar measurement invariance violation, 
                                  # correlation between group and item X9!

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) 
# group is uncorrelated with outcome; other variables around ~0.25

groupIdx <- which(colnames(data) == "group")
cor(data[,"Y"], rowMeans(data[,-c(groupIdx, outcomeIdx)]))
# 0.33

# for the multiple group analysis, the data looks fine in both groups (see first test above)
fit_multipleGroups <- lavaan::sem(m0, data = data, group = "group")
fit_multipleGroups@Model@GLIST

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
fit_singleGroup@Model@GLIST

##### correlate confounding group with outcome #####
# in every MI violating condition, group is not correlated with the outcome
# correlation between group and items:
#   - configural and metric MI violations do not create correlations between items and group
#   - scalar MI violation creates correlation between MI violating item and the group

# simulate from SEM without any MI violation
nGroups <- 2

# factor loadings
lambdaSEM <- defineLambda(nGroups, nItems = J, miItems = 0, miWeight = list(c(NA)), params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)

# factor variances
covSEM <- definePhi(nFactors, Sigma_T, Sigma_S[1, 1], nItems = J, scalarMI = F, 
                    varNames = colnames(lambdaSEM))
covSEM

# path coefficient matrix (effect of factors on the target/criterion)
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- defineBeta(nFactors, betaT, rep(0, J), J, colnames(lambdaSEM))

# residual variance matrix of the items 
# nVariables <- dim(lambdaSEM)[1]
rownames(lambdaSEM)[grep("^X", rownames(lambdaSEM))]
itemNameIdxVec <- grep("^X", rownames(lambdaSEM))

sapply(itemNameIdxVec, function(iVar) {
  iVarName <- rownames(lambdaSEM)[iVar]  
  facVec <- names(which(lambdaSEM[iVarName,] != 0))
  getErrorVar(loadingVec = lambdaSEM[iVarName,facVec], 
              varianceVec = diag(covSEM[facVec, facVec]),
              targetVar = 1)
})

# check for one variable if error matches to meet the targetVar
facVec <- names(which(lambdaSEM["X9_g2",] != 0))
getErrorVar(loadingVec = lambdaSEM["X9_g2",facVec], 
            varianceVec = diag(covSEM[facVec, facVec]),
            targetVar = 1)

thetaSEM <- diag(c(rep(diag(Sigma_e), nGroups), rep(0, params$nOutcome), rep(0, params$nOutcome))) 
thetaSEM


# calculate model implied covariance matrix
# for simulating from SEM
# when simulating two outcomes instead of one from the block matrix, the outcome
#   for each block has the same variance as for one group (var(Y1) = var(Y2) = 1.125)
sigmaSEM <- getSigma(lambdaSEM, betaMat, covSEM, thetaSEM)

# sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
#   t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM


# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
Y1g1Idx <- which(colnames(sigmaSEM) == "Y1_g1")
Y1g2Idx <- which(colnames(sigmaSEM) == "Y1_g2")

# we assume balance design with 2 groups! 
varWithinY <- 0.5 * sigmaSEM[Y1g1Idx, Y1g1Idx] + 0.5 * sigmaSEM[Y1g2Idx, Y1g2Idx]
targetRho <- 0.5
meanYg2 <- getMeanYg2(targetRho, params$dgp$varGroup, varWithinY)

M <- c(rep(params$mS, J*nGroups), rep(0, params$nOutcome), rep(meanYg2, params$nOutcome)) 
M

# this function needs to enable ...
#     ... mean differences for scalar MI violations and 
#     ... mean differences between groups in the outcome to simulate their correlation
defineMeanVec <- function(scalarMI = F) {
  
}

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(J), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(nGroups), each = J)), "Y_g1", "Y_g2")

# checḱ variances of all items
#   -> variances in all items (independent of group) are 1
apply(data_multipleGroup, 2, var)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data_multipleGroup) == "Y_g1" | colnames(data_multipleGroup) == "Y_g2")
cor(data_multipleGroup[,"Y_g1"], data_multipleGroup[,-outcomeIdx]) # ~ around 0.2 per item
cor(data_multipleGroup[,"Y_g2"], data_multipleGroup[,-outcomeIdx]) # ~ around 0.2 per item

cor(data_multipleGroup[,"Y_g1"], rowMeans(data_multipleGroup[,-outcomeIdx])) # 0.21
cor(data_multipleGroup[,"Y_g2"], rowMeans(data_multipleGroup[,-outcomeIdx])) # 0.16

# wide to long
data <- data.frame(data_multipleGroup)
data <- tidyr::pivot_longer(data, cols = !c(Y_g1, Y_g2), names_to = "item", values_to = "value")
data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
data$Y <- ifelse(data$group == "g1", data$Y_g1, data$Y_g2)

##### without functions #####
# lambda_g1 <- cbind(TS = rep(1, J),
#                    specificEffects = diag(J))
# lambda_g2 <- cbind(TS = rep(1, J),
#                    specificEffects = diag(J))
# # without group
# lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
#                          rep(0, params$nOutcome),
#                          rep(0, params$nOutcome)),
#                    c(rep(0, (J + 1)*nGroups), rep(1, params$nOutcome), 0),
#                    c(rep(0, (J + 1)*nGroups), 0, rep(1, params$nOutcome)))
# colnames(lambdaSEM) <- c(paste0(rep(c("T", 
#                                       sapply(seq(J), 
#                                              function(iJ) {paste0("S", iJ)})), times = 2), 
#                                 "_g", rep(seq_len(nGroups), each = J+1)), "Y_g1", "Y_g2")
# covSEM <- diag(nFactors) * c(
#   rep(c(Sigma_T, diag(Sigma_S)), times = 2),  
#   params$sigmaErrY, params$sigmaErrY) # residual variance!
# betaMat <- matrix(rep(0, nFactors + nFactors * (nFactors-1)), 
#                   ncol = nFactors, byrow = T)
# # betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0, 0)
# betaMat[nFactors-1,] <- c(rep(c(betaT, 0, 0), c(1, J, J+params$M)), 0, 0)
# betaMat[nFactors,] <- c(rep(c(0, betaT, 0), c(J+params$M, 1, J)), 0, 0)
# betaMat
##### check data #####
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
groupIdx <- which(colnames(data) == "group")
cor(data$group, data[,-groupIdx]) # correlations with group  are all zero

cov(data$group, data$Y)
mean(data$group*data$Y)

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) # group is uncorrelated with outcome; other variables around ~0.2

groupIdx <- which(colnames(data) == "group" | 
                    colnames(data) == "Y_g1" |
                    colnames(data) == "Y_g2") 
cor(data[,"Y"], rowMeans(data[,-c(groupIdx, outcomeIdx)]))
# 0.29

# for the multiple group analysis, the data looks fine in both groups (see first test above)
fit_multipleGroups <- lavaan::sem(m0, data = data[,-c(1, 2)], group = "group")
fit_multipleGroups@Model@GLIST

# with identical factor models, modelling both groups at once works as well
fit_singleGroup <- lavaan::sem(m0, data = data)
fit_singleGroup@Model@GLIST


##### simulate from multiple group SEM with group--outcome correlation and configural MI violation #####
nGroups <- 2
miItems <- 1
miWeight <- list(c(0))

lambdaSEM <- defineLambda(params$dgp$nGroups, J, miItems = miItems, 
                          miWeight = miWeight, params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)

covSEM <- definePhi(nFactors, Sigma_T, Sigma_S[1, 1], J)

defineBeta <- function(nFactors,) {
  # beta matrix without any effect
  betaMat <- matrix(rep(0, nFactors + nFactors*(nFactors-1)), ncol = nFactors, byrow = T)
  
  # outcome variables are sorted in the following order: Y1_g1 Y2_g1 Y2_g2 Y2_g2
  params$dgp$nGroups
}


# betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0, 0)
betaMat[nFactors-1,] <- c(rep(c(betaT, 0, 0), c(1, J, J+params$M)), 0, 0)
betaMat[nFactors,] <- c(rep(c(0, betaT, 0), c(J+params$M, 1, J)), 0, 0)
betaMat

betaMat <- matrix(rep(0, 
                      ((J + params$M)*nGroups + params$nOutcome) + # diagonal
                        (((J + params$M)*nGroups + params$nOutcome)*((J + params$M)*nGroups + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(betaT, beta_S), times = 2), 0)

##### simulate from multiple group SEM with group--outcome correlation and metric MI violation #####
nGroups <- 2
miItems <- 1
miWeight <- list(c(-0.8))
lambdaSEM <- defineLambda(params$dgp$nGroups, J, miItems = miItems, 
                          miWeight = miWeight, params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)

covSEM <- definePhi(nFactors, Sigma_T, Sigma_S[1, 1], J)

##### simulate from multiple group SEM with group--outcome correlation and scalar MI violation #####

covSEM <- definePhi(nFactors, Sigma_T, Sigma_S[1, 1], J, scalarMI = T, miItems = 1)


# test correlation between group and outcome depending on outcome mean differences
#     between group for all MI violation conditions

##### error calculation for multiple items with ME violation #####

##### simulate nuance factor (without groups!)

################################################################################
# fit data
################################################################################
# N = {50, 100, 250, 500, 1000}
# simulated R² = {0.1, 0.3, 0.5, 0.8}
#   -> this will either be manipulated by manipulating beta coefficients or error variance of the outcome
# different dgp:
#   ... without MI violations, configural MIvio, metric MIvio, scalar MIvio
#   ... correlation between outcome and MI violation inducing group: r = {0, 0.2, 0.5}
#   ... with and without specific item effect
#   ... with and without nuance factor 


# prediction based on linear regression and based on elastic net regression
#   ... prediction based on items and based on composite
#   ... prediction based on predictor set with group, without group, with group predicted from items

# dependent measures
#   ... estBeta
#   ... R^2

######
# # calculate the composite score
# #     for the composite we do not use true W but unit weights for all items in every group
# Xcomp = data[, paste0("X", seq_len(J))] %*% W # this is the composite (out of all manifest variables); matrix of composites scores weighted by W (NxM)
# 
# data <- data.frame(data)
# # Creating Train-Test sets
# # X = (manifest) item values = T + S + e 
# itms.df = data.frame(data$Y, data[, paste0("X", seq_len(J))])
# itms.train=itms.df[1:N.train,]
# itms.test=itms.df[(N.train+1):N,]
# 
# # ("manifest") composite values = (1*T1 + 1*T2 + ... + 1*Tj) + (S1 + S2 + ... + Sj) + 
# comp.df = data.frame(data$Y,Xcomp)
# comp.train=comp.df[1:N.train,]
# comp.test=comp.df[(N.train+1):N,]
# 
# # combine items + construct
# comb.df = data.frame(data$Y, data[, paste0("X", seq_len(J))], Xcomp)
# comb.train = comb.df[1:N.train,]
# comb.test = comb.df[(N.train+1):N,]
# 
# # true prediction by using Construct and specific item effect both as predictors
# #   ->  this is not feasible anymore; specific item effects are simulated as latent 
# #       factors that are not identified and therefore cannot be estimated within the model
# 
# # Regressions
# lm.items = lm(data.Y~.,itms.train)
# lm.comp = lm(data.Y~.,comp.train)
# lm.comb = lm(data.Y~.,comb.train)
# 
# # observed coefficients excluding intercept
# items.coef = coef(lm.items)[2:(J + 1)]
# comp.coef = coef(lm.comp)[2:(params$M + 1)]
# comb.coef = coef(lm.comb)[2:(J + params$M + 1)]
# 
# # Compute the w/in sample bias for both models at the observation level
# #   ->  this is not feasible anymore; specific item effects are simulated as latent 
# #       factors that are not identified and therefore cannot be estimated within the model
# 
# #Calculate Holdout Prediction
# itms.pred = predict(lm.items,itms.test)
# comp.pred = predict(lm.comp,comp.test)
# comb.pred = predict(lm.comb,comb.test)
# 
# #Preparing Output
# mat = rep(NA, 28)
# 
# #Condition Information
# mat[1] = seed; mat[2] = N.train; mat[3] = J #note that N -holdout is returned not total N
# mat[4] = betaS;mat[5] = betaT; mat[6] = measVar.cond #beta_S, beta_T, and measVar conditions
# mat[7] = nGroups; mat[8] = miItems; mat[9] = miWeight; mat[10] = mSpecific
# #Bias
# mat[11] = NA
# mat[12] = NA
# mat[13] = NA
# 
# #Variance
# mat[14] = NA
# mat[15] = NA
# mat[16] = NA
# 
# #within samp error
# mat[17] = NA
# mat[18] = NA
# mat[19] = NA
# 
# #Holdout RMSE
# mat[20] = calcRMSE(itms.pred,itms.test[,"data.Y"])
# mat[21] = calcRMSE(comp.pred,comp.test[,"data.Y"])
# mat[22] = calcRMSE(comb.pred,comb.test[,"data.Y"])
# mat[23] = NA
# 
# # test R2
# mat[24] <- evalPerformance(itms.pred,itms.test[,"data.Y"])["Rsquared"]
# mat[25] <- evalPerformance(comp.pred,comp.test[,"data.Y"])["Rsquared"]
# mat[26] <- evalPerformance(comb.pred,comb.test[,"data.Y"])["Rsquared"]
# mat[27] <- NA
# 
# names(mat) <-  c("Seed","N","J","beta_S","beta_T","measVar",
#                  "nGroups", "miItems", "miWeight", "mSpecific",
#                  "item.bias","comp.bias","true.bias",
#                  "item.var","comp.var","true.var",
#                  "item.err","comp.err","true.err",
#                  "item.rmse","comp.rmse", "comb.rmse", "true.rmse",
#                  "item.R2", "comp.R2", "comb.R2", "true.R2",
#                  "sim.R2")
# # SEM is not identified; thus, we cannot estimate SEM to check simulation
# # To Do: simulate MI violations in SEM framework
# 
# m0 <- '# measurement model
#        T =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7 + 1*X8 + 1*X9
#        S1 =~ 1*X1
#        S2 =~ 1*X2
#        S3 =~ 1*X3
#        S4 =~ 1*X4
#        S5 =~ 1*X5
#        S6 =~ 1*X6
#        S7 =~ 1*X7
#        S8 =~ 1*X8
#        S9 =~ 1*X9
# 
#        # uncorrelated factor and nuances
#        T ~~ 0*S1
#        T ~~ 0*S2
#        T ~~ 0*S3
#        T ~~ 0*S4
#        T ~~ 0*S5
#        T ~~ 0*S6
#        T ~~ 0*S7
#        T ~~ 0*S8
#        T ~~ 0*S9
#        S1 ~~ 0*S2
#        S1 ~~ 0*S3
#        S1 ~~ 0*S4
#        S1 ~~ 0*S5
#        S1 ~~ 0*S6
#        S1 ~~ 0*S7
#        S1 ~~ 0*S8
#        S1 ~~ 0*S9
#        S2 ~~ 0*S3
#        S2 ~~ 0*S4
#        S2 ~~ 0*S5
#        S2 ~~ 0*S6
#        S2 ~~ 0*S7
#        S2 ~~ 0*S8
#        S2 ~~ 0*S9
#        S3 ~~ 0*S4
#        S3 ~~ 0*S5
#        S3 ~~ 0*S6
#        S3 ~~ 0*S7
#        S3 ~~ 0*S8
#        S3 ~~ 0*S9
#        S4 ~~ 0*S5
#        S4 ~~ 0*S6
#        S4 ~~ 0*S7
#        S4 ~~ 0*S8
#        S4 ~~ 0*S9
#        S5 ~~ 0*S6
#        S5 ~~ 0*S7
#        S5 ~~ 0*S8
#        S5 ~~ 0*S9
#        S6 ~~ 0*S7
#        S6 ~~ 0*S8
#        S6 ~~ 0*S9
#        S7 ~~ 0*S8
#        S7 ~~ 0*S9
#        S8 ~~ 0*S9
# 
#        # all latent representation of manifest outcome
#        fY =~ 1 * Y
#        Y ~~ 0*Y
#        
#        # structure model
#        fY ~ T + S1'
# 
# m0 <- '# measurement model
#        T =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7 + 1*X8 + 1*X9
#        S1 =~ 1*X1
#        S2 =~ 1*X2
#        S3 =~ 1*X3
#        S4 =~ 1*X4
#        S5 =~ 1*X5
#        S6 =~ 1*X6
#        S7 =~ 1*X7
#        S8 =~ 1*X8
#        S9 =~ 1*X9
# 
#        # all latent representation of manifest outcome
#        fY =~ 1 * Y
#        Y ~~ 0*Y
# 
#        # structure model
#        fY ~ T + S1'

##### simulate helper item for specific item effect factor #####
# # simulate helper item that does not load on the factor to have identied specific effect
# #   factors and check if the variance is simulated correctly
# lambdaSEM <- cbind(trueScore = c(rep(c(1, 0), J), rep(0, params$nOutcome)),
#                    specificEffects = rbind(apply(diag(J), 1, rep, each = 2), 0),
#                    Y = c(rep(0, J*2), rep(1, params$nOutcome)))
# nFactors <- ncol(lambdaSEM)
# colnames(lambdaSEM) <- c("T", 
#                          sapply(seq(J), function(iJ) {paste0("S", iJ)}),
#                          "Y")
# rownames(lambdaSEM) <- c(paste0("X", rep(seq_len(J), each = 2), "_", rep(c("true", "helper"), times = J)), "Y")
# lambdaSEM
# 
# # variances for T and S from the Jacobucci Paper
# #   variance of the outcome is unknwon but we only need residual variance of the
# #   outcome since y is an endogenuous variable in the structural model
# # factors are not correlated 
# #   T, S and e are sampled from different multivariate normal distributions, thus,
# #   factors need to be independent/uncorrelated
# #   relating the model to substantive research: nuance information is independent 
# #   from facet information (this is the dgm, i.e., extra information in nuances 
# #   that has nothing to do with what is already explained by facet)
# # variables from structure model need to be included in factor correlation matrix 
# # for the SEM simulation 
# covSEM <- diag(nFactors) * c(Sigma_T, diag(Sigma_S), params$sigmaErrY)
# covSEM
# 
# # # these is the variance of the composite with the true score entering multiple times
# # # varFormComp <- 40.5873
# # varFormComp <- 45
# # covSEM <- diag(nFactors) * c(varFormComp, diag(Sigma_S), params$sigmaErrY)
# 
# # path coefficient matrix (effect of factors on the target/criterion)
# # ! factors in order of appearance within model syntax
# # in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
# betaMat <- matrix(rep(0, 
#                       (J + params$M + params$nOutcome) + # diagonal
#                         ((J + params$M + params$nOutcome)*(J + params$M + params$nOutcome - 1))), # off diagonal
#                   ncol = nFactors, byrow = T)
# betaMat[nFactors,] <- c(betaT, beta_S, 0)
# 
# # error - Kovarianzmatrix der Residualvariablen mit p x p (= Anzahl der manifesten Variablen)
# # outcome as manifest variable but all latent! (therefore the error variance is zero?)
# thetaSEM <- diag(c(rep(
#   c(params$sigma[measVar.cond,"sigmaE"], 0), 
#   # c(params$sigma[measVar.cond,"sigmaE"], 1-params$sigma[measVar.cond,"sigmaS"]), 
#   times = J),  
#   rep(0, params$nOutcome)))
# 
# # set mean values for manifest variables
# # check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
# M <- c(rep(params$mS, J*2), rep(0, params$nOutcome))
# 
# # calculate model implied covariance matrix
# # for simulating from SEM
# sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
#   t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
# 
# # simulate date from multivariate normal distribution
# set.seed(seed)
# dataHelperItems <- mvtnorm::rmvnorm(N.train + params$holdout, mean = M, sigma = sigmaSEM)
# colnames(dataHelperItems) <- c(paste0("X", rep(seq_len(J), each = 2), "_", rep(c("true", "helper"), times = J)), "Y")
# 
# head(dataHelperItems)
# apply(dataHelperItems, 2, var) # this looks fine
# 
# mSpecific <- '# measurement model
#               S1 =~ X1_true + X1_helper
#               S2 =~ X2_true + X2_helper
#               S3 =~ X3_true + X3_helper
#               S4 =~ X4_true + X4_helper
#               S5 =~ X5_true + X5_helper
#               S6 =~ X6_true + X6_helper
#               S7 =~ X7_true + X7_helper
#               S8 =~ X8_true + X8_helper
#               S9 =~ X9_true + X9_helper
#               T =~ X1_true + X2_true + X3_true + X4_true + X5_true + X6_true + X7_true + X8_true + X9_true
#               
#               # all latent representation of manifest outcome
#               fY =~ 1 * Y
#               Y ~~ 0*Y
#               
#               T ~~ 0*S1
#               T ~~ 0*S2
#               T ~~ 0*S3
#               T ~~ 0*S4
#               T ~~ 0*S5
#               T ~~ 0*S6
#               T ~~ 0*S7
#               T ~~ 0*S8
#               T ~~ 0*S9
#               S1 ~~ 0*S2
#               S1 ~~ 0*S3
#               S1 ~~ 0*S4
#               S1 ~~ 0*S5
#               S1 ~~ 0*S6
#               S1 ~~ 0*S7
#               S1 ~~ 0*S8
#               S1 ~~ 0*S9
#               S2 ~~ 0*S3
#               S2 ~~ 0*S4
#               S2 ~~ 0*S5
#               S2 ~~ 0*S6
#               S2 ~~ 0*S7
#               S2 ~~ 0*S8
#               S2 ~~ 0*S9
#               S3 ~~ 0*S4
#               S3 ~~ 0*S5
#               S3 ~~ 0*S6
#               S3 ~~ 0*S7
#               S3 ~~ 0*S8
#               S3 ~~ 0*S9
#               S4 ~~ 0*S5
#               S4 ~~ 0*S6
#               S4 ~~ 0*S7
#               S4 ~~ 0*S8
#               S4 ~~ 0*S9
#               S5 ~~ 0*S6
#               S5 ~~ 0*S7
#               S5 ~~ 0*S8
#               S5 ~~ 0*S9
#               S6 ~~ 0*S7
#               S6 ~~ 0*S8
#               S6 ~~ 0*S9
#               S7 ~~ 0*S8
#               S7 ~~ 0*S9
#               S8 ~~ 0*S9
#               
#               # structure model
#               fY ~ T'
# 
# fit_helperItems <- lavaan::sem(mSpecific, data = dataHelperItems)
# fit_helperItems@Model@GLIST
# fit_helperItems@Model@GLIST$lambda