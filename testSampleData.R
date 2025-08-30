################################################################################
# simulate data from SEM 
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

# step by step approach
# 1. simulate specific item effect without measurement invariance violations
# 2. simulate from multiple group SEM without any MI violation
# 3. simulate from multiple group SEM with configural MI violation
# 4. simulate from multiple group SEM with metric measurement invariance violation
# 5. simulate from multiple group SEM with scalar measurement invariance violation
# 6. add correlation between MI violating group variable and the outcome
#     -> start simulating Y for each group
# 7. test MI violations for multiple Items

################################################################################
# set up
################################################################################
source('utils/simTools.R')
source('utils/setParam.R')
params$dgp$mS <- params$dgp$mItems # changed the variable name and do not want to adapt it here...

require(lavaan)
require(tidyr)
require(mvtnorm)
library(magic)

# check SEM data
m0 <- '# measurement model
       T =~ X1 + X2 + X3 + X4

       # all latent representation of manifest outcome 
       #    loading = 1, residual variance of manifest variable = 0
       fY =~ 1 * Y
       Y ~~ 0*Y
       
       # structure model
       fY ~ T'

################################################################################
# 0. simulate only measurement model 
################################################################################
# factor loading matrix without outcome; factors in columns, items in rows
# -> to simulate specific item effects, the outcome has to be simulated within the
#    SEM context
# params$dgp$tauFL[5]: p = 4, cAlpha = 0.9
lambda <- cbind(trueScore = rep(params$dgp$tauFL[5], params$dgp$P[1]), 
                specificEffects = t(diag(params$dgp$P[1])))
nFactors <- ncol(lambda)
covMat <- diag(nFactors) * c(params$dgp$varT[1], rep(params$dgp$varS[1], params$dgp$P[1]))

varE <- getErrorVar(c(params$dgp$tauFL[5], 1), c(params$dgp$varT[1], params$dgp$varS[1]), 1)
theta <- diag(params$dgp$P[1]) * varE
# for simulating only from measurement model
sigmaMat <- lambda %*% covMat %*% t(lambda) + theta

# variance in Y is way smaller in SEM simulation compared to Jacobucci simulation
#   betaS.cond == 1: var(Y) = 12.31 vs. var(data$Y) = 2.16154
#   betaS.cond == 4: var(Y) = 11.80889 vs. var(data$Y) = 1.178903
#   ->  this is due to the rather "formative" conception of the true score! 
#       the true score enters y multiple times in the original Jacobucci simulation
#   -> now that we simulate a reflective model in the original simulation the variance
#      matches with the sem formulation of the data 



# to do: 
# test multiple group SEM with grouping variable for all MI violation conditions
# test MI violations for multiple Items
# write functions for data generating model


################################################################################
# 1. simulate from SEM 
# (with lots of explanations; without multiple groups + without MI violations)
################################################################################
# factor loading matrix with outcome; but variance of the outcome is unknwon
#   therefore: factor loading matrix without outcome below

lambdaSEM <- defineLambda(flT = params$dgp$varFLgrid$flT[5], flS = params$dgp$varFLgrid$flS[5], 
                          nGroups = 1, nItems = params$dgp$varFLgrid$P[5], 
                          miItems = 0, miWeight = NA, nOutcome = params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)

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
covSEM <- diag(nFactors) * c(params$dgp$varFLgrid$varT[5], 
                             rep(params$dgp$varFLgrid$varS[5], params$dgp$varFLgrid$P[5]), 
                             params$dgp$sigmaEy)
covSEM

# path coefficient matrix (effect of factors on the target/criterion)
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- matrix(rep(0, 
                      (params$dgp$varFLgrid$P[5] + params$dgp$Q + params$dgp$nOutcome) + # diagonal
                        ((params$dgp$varFLgrid$P[5] + params$dgp$Q + params$dgp$nOutcome)*(params$dgp$varFLgrid$P[5] + params$dgp$Q + params$dgp$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(params$dgp$betaT, rep(params$dgp$betaS, params$dgp$varFLgrid$P[5]), 0)

# error - Kovarianzmatrix der Residualvariablen mit p x p (= Anzahl der manifesten Variablen)
# outcome as manifest variable but all latent! (therefore the error variance is zero?)
thetaSEM <- diag(c(rep(params$dgp$varFLgrid$varE[5], params$dgp$varFLgrid$P[5]), 
                   rep(0, params$dgp$nOutcome))) # without helper items

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
M <- c(rep(params$dgp$mS, params$dgp$varFLgrid$P[5]), 
       rep(0, params$dgp$nOutcome)) # without helper items

# calculate model implied covariance matrix
# for simulating from SEM
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
set.seed(7829)
data <- mvtnorm::rmvnorm(params$sim$nTrain + params$sim$testN, 
                         mean = M, sigma = sigmaSEM)
colnames(data) <- c(paste0("X", seq_len(params$dgp$varFLgrid$P[5])), "Y") # without helper items

# get correlation between Y and other variables
outcomeIdx <- which(colnames(data) == "Y")
cor(data[,"Y"], data[,-outcomeIdx]) 
# old: 0.228212 0.2486797 0.2208203 0.2284228 0.2128141 0.2395708 0.2245688 0.2379814 0.2142187
# new: 0.2710696 0.2213895 0.2268947 0.2417323

cor(data[,"Y"], rowMeans(data[,-outcomeIdx]))
# old: 0.3
# new: 0.3046089

# check SEM data
m0 <- '# measurement model
       T =~ X1 + X2 + X3 + X4

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

# score variance and Cronbach's alpha a slightly below
# get variance of the score; 4 + 2* (4*(4-1)/2) * 0.5 = 10
(scoreVar <- var(rowSums(data[,seq_len(params$dgp$varFLgrid$P[5])], na.rm = TRUE), na.rm = TRUE))

# get Cronbach's Alpgha; # simulated 0.8 so close enough?
getCronbachsAlpha(data[,seq_len(params$dgp$varFLgrid$P[5])]) 

################################################################################
# 2. simulate from multiple group SEM without any MI violation
# (with lots of explanations; without MI violations)
################################################################################
# block-diagonal binding of matrices 
# nGroups <- 2
lambda_g1 <- cbind(TS = rep(params$dgp$varFLgrid$flT[5], params$dgp$varFLgrid$P[5]),
                   specificEffects = diag(params$dgp$varFLgrid$P[5]))
lambda_g2 <- cbind(TS = rep(params$dgp$varFLgrid$flT[5], params$dgp$varFLgrid$P[5]),
                   specificEffects = diag(params$dgp$varFLgrid$P[5]))
# without group
lambdaSEM <- rbind(cbind(magic::adiag(lambda_g1, lambda_g2), 
                         rep(0, params$dgp$nOutcome)),
                   c(rep(0, (params$dgp$varFLgrid$P[5] + 1)*params$dgp$nGroups), rep(1, params$dgp$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(params$dgp$varFLgrid$P[5]), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[5]+1)), "Y")

# ! defineLambda function cannot model only one outcome despite multiple groups

# # with group
# lambdaSEMg <- rbind(cbind(adiag(lambda_g1, lambda_g2),
#                          # group = rep(c(0, 1), each = J),
#                          group = rep(c(0, 0), each = params$dgp$varFLgrid$P[5]),
#                          rep(0, params$dgp$nOutcome)),
#                    c(rep(c(0, 1), each = (params$dgp$varFLgrid$P[5] + 1)), 1, rep(0, params$dgp$nOutcome)),
#                    c(rep(0, (params$dgp$varFLgrid$P[5] + 1)*params$dgp$nGroups + 1), rep(1, params$dgp$nOutcome)))
# colnames(lambdaSEMg) <- c(paste0(rep(c("T",
#                                       sapply(seq(params$dgp$varFLgrid$P[5]),
#                                              function(iJ) {paste0("S", iJ)})), times = params$dgp$nGroups),
#                                 "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[5]+1)),
#                          "group", "Y")

nFactors <- ncol(lambdaSEM)

# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
# without group
covSEM <- diag(nFactors) * c(
  rep(c(params$dgp$varFLgrid$varT[5], rep(params$dgp$varFLgrid$varS[5], params$dgp$varFLgrid$P[5])), times = 2),  
  params$dgp$sigmaEy) # residual variance!

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
                      ((params$dgp$varFLgrid$P[5] + params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome) + # diagonal
                        (((params$dgp$varFLgrid$P[5] + params$dgp$Q)*params$dgp$nGroups + 
                            params$dgp$nOutcome)*((params$dgp$varFLgrid$P[5] + 
                                                     params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(params$dgp$betaT, rep(params$dgp$betaS, params$dgp$varFLgrid$P[5])), times = 2), 0)

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
thetaSEM <- diag(c(rep(params$dgp$varFLgrid$varE[5], 
                       params$dgp$varFLgrid$P[5] * params$dgp$nGroups), 
                   rep(0, params$dgp$nOutcome))) # without helper items

# # with group
# thetaSEM <- diag(c(rep(diag(Sigma_e), nGroups), 0, rep(0, params$nOutcome))) # without helper items
thetaSEM

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
M <- c(rep(params$dgp$mS, params$dgp$varFLgrid$P[5]*params$dgp$nGroups), rep(0, params$dgp$nOutcome)) # without helper items

# # with group
# M <- c(rep(params$mS, J*nGroups), 0, rep(0, params$nOutcome)) # without helper items

# calculate model implied covariance matrix
# for simulating from SEM
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(93901)
data_multipleGroup <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, 
                                       mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(params$dgp$varFLgrid$P[5]), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[5])), "Y")

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

# check data 
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
apply(data, 2, var)
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

################################################################################
# 3. simulate from multiple group SEM with configural measurement invariance violation 
# (remove the factor loadings of one of the items)
################################################################################
# when removing the loading of 1 item on the true score in one of the groups, this
#   item decreases in variance
# introduce MI violation in last item to avoid using it as marker item in SEM
# nGroups <- 2
# miItems <- 1
# miWeight <- 0
lambda_g1 <- cbind(TS = rep(params$dgp$varFLgrid$flT[5], params$dgp$varFLgrid$P[5]),
                   specificEffects = diag(params$dgp$varFLgrid$P[5]))
lambda_g2 <- cbind(TS = rep(c(params$dgp$varFLgrid$flT[5], params$dgp$miWeight[1]), 
                            c(params$dgp$varFLgrid$P[5] - params$dgp$miItems[1], params$dgp$miItems[1])),
                   specificEffects = diag(params$dgp$varFLgrid$P[5]))

lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                         rep(0, params$dgp$nOutcome)),
                   c(rep(0, (params$dgp$varFLgrid$P[5] + 1)*params$dgp$nGroups), rep(1, params$dgp$nOutcome)))
colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                      sapply(seq(params$dgp$varFLgrid$P[5]), 
                                             function(iJ) {paste0("S", iJ)})), times = 2), 
                                "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[5]+1)), "Y")


nFactors <- ncol(lambdaSEM)


# variances for T and S from the Jacobucci Paper
#   variance of the outcome is unknwon but we only need residual variance of the
#   outcome since y is an endogenuous variable in the structural model
covSEM <- diag(nFactors) * c(
  rep(c(params$dgp$varFLgrid$varT[5], rep(params$dgp$varFLgrid$varS[5], params$dgp$varFLgrid$P[5])), times = 2),  
  params$dgp$sigmaEy) # residual variance!

# path coefficients
betaMat <- matrix(rep(0, 
                      ((params$dgp$varFLgrid$P[5] + params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome) + # diagonal
                        (((params$dgp$varFLgrid$P[5] + params$dgp$Q)*params$dgp$nGroups + 
                            params$dgp$nOutcome)*((params$dgp$varFLgrid$P[5] + params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(rep(c(params$dgp$betaT, rep(params$dgp$betaS, params$dgp$varFLgrid$P[5])), times = params$dgp$nGroups), 0)

thetaSEM <- diag(c(rep(params$dgp$varFLgrid$varE[5], params$dgp$varFLgrid$P[5]),  # group 1
                   rep(c(params$dgp$varFLgrid$varE[5], 
                         params$dgp$varFLgrid$varE[5]+params$dgp$varFLgrid$varT[5]), 
                       c(params$dgp$varFLgrid$P[5]-params$dgp$miItems[1], params$dgp$miItems[1])), # group 2
                   rep(0, params$dgp$nOutcome))) # outcome

sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(28930)
data_multipleGroup <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, 
                                       mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(params$dgp$varFLgrid$P[5]), 
                                                      function(iJ) {paste0("X", iJ)})), times = 2), 
                                         "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[5])), "Y")

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

# check data 
# without specific item effect on outcome
# check variances across groups now
apply(data, 2, var)

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
apply(data, 2, var)
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

################################################################################
# 4. simulate from multiple group SEM with metric measurement invariance violation 
################################################################################
# nGroups <- 2
# miItems <- 1
# miWeight <- -0.8

# # this works: 
# iCond <- 41 # works too: 37, 38, 39, 40, 41, 42
# miWeight <- params$dgp$varFLgrid$flT[iCond] +
#   (params$dgp$varFLgrid$b1[iCond]*(params$dgp$varFLgrid$varT[iCond]/params$dgp$varItems)) * params$dgp$varFLgrid$flT[iCond]

# # this works: 43, 44, 45; for varT = 0.5 this is the same as above
# iCond <- 53, 54, # zu viel Varianz in x4
# iCond <- 54
# miWeight <- params$dgp$varFLgrid$flT[iCond] + (params$dgp$varFLgrid$b1[iCond] * params$dgp$varFLgrid$flT[iCond])
# miWeight <- params$dgp$varFLgrid$flT[iCond] +
#   ((params$dgp$varFLgrid$b1[iCond] * params$dgp$varFLgrid$flT[iCond]) * (params$dgp$varItems - params$dgp$varFLgrid$varT[iCond]))
# miWeight <- params$dgp$varFLgrid$flT[iCond] +
#   (params$dgp$varFLgrid$b1[iCond]*(params$dgp$varFLgrid$varT[iCond]/params$dgp$varItems)) * params$dgp$varFLgrid$flT[iCond]

# iCond <- 53
# miWeight <- params$dgp$varFLgrid$flT[iCond] + 
#   (params$dgp$varFLgrid$b1[iCond]*(params$dgp$varItems - params$dgp$varFLgrid$varT[iCond])) * params$dgp$varFLgrid$flT[iCond]

condVec <- which(# params$dgp$grid$P == 4 &
                 #params$dgp$grid$cAlpha == 0.8 &
                 params$dgp$grid$miItems > 1 &
                 #params$dgp$grid$varT == 0.5 &
                 params$dgp$grid$MIvio == "metric")

for (iCond in condVec) {
  miWeight <- params$dgp$grid$flT[iCond] + params$dgp$grid$b1[iCond]
  lambda_g1 <- cbind(TS = rep(params$dgp$grid$flT[iCond], params$dgp$grid$P[iCond]),
                     # specificEffects = diag(params$dgp$varFLgrid$P[iCond]))
                     specificEffects = diag(params$dgp$grid$P[iCond]) * 
                       rep(c(1, params$dgp$grid$flS[iCond]), 
                           c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])))
  
  lambda_g2 <- cbind(TS = rep(c(params$dgp$grid$flT[iCond], miWeight), 
                              c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])),
                     specificEffects = diag(params$dgp$grid$P[iCond]) * 
                       rep(c(1, params$dgp$grid$flS[iCond]), 
                           c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])))
  
  lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                           rep(0, params$dgp$nOutcome)),
                     c(rep(0, (params$dgp$grid$P[iCond] + 1)*params$dgp$nGroups), 
                       rep(1, params$dgp$nOutcome)))
  colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                        sapply(seq(params$dgp$grid$P[iCond]), 
                                               function(iJ) {paste0("S", iJ)})), times = 2), 
                                  "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond]+1)), "Y")
  rownames(lambdaSEM) <- c(paste0(rep(c(sapply(seq(params$dgp$grid$P[iCond]), 
                                               function(iJ) {paste0("X", iJ)})), times = 2), 
                                  "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond])), "Y")
  
  print(lambdaSEM)
  
  nFactors <- ncol(lambdaSEM)
  
  # variances for T and S from the Jacobucci Paper
  # #   variance of the outcome is unknwon but we only need residual variance of the
  # #   outcome since y is an endogenuous variable in the structural model
  # covSEM <- diag(nFactors) * c(
  #   rep(c(params$dgp$varFLgrid$varT[iCond], rep(params$dgp$varFLgrid$varS[iCond], params$dgp$varFLgrid$P[iCond])), times = 2),  
  #   params$dgp$sigmaEy) # residual variance!
  
  covSEM <- diag(nFactors) * c(
    rep(c(params$dgp$grid$varT[iCond], params$dgp$grid$varStable[iCond], params$dgp$grid$varS[iCond]), 
        c(1, params$dgp$grid$P[iCond] - params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])), # group 1
    rep(c(params$dgp$grid$varT[iCond], params$dgp$grid$varStable[iCond], params$dgp$grid$varS[iCond]), 
        c(1, params$dgp$grid$P[iCond] - params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])), # group 2
    params$dgp$sigmaEy)
  print(covSEM)
  # covSEM[10, 10] <- 0
  # covSEM[5, 5] <- 0
  
  # path coefficients
  betaMat <- matrix(rep(0, 
                        ((params$dgp$grid$P[iCond] + params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome) + # diagonal
                          (((params$dgp$grid$P[iCond] + params$dgp$Q)*params$dgp$nGroups + 
                              params$dgp$nOutcome)*((params$dgp$grid$P[iCond] + 
                                                       params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome - 1))), # off diagonal
                    ncol = nFactors, byrow = T)
  betaMat[nFactors,] <- c(rep(c(params$dgp$betaT, rep(params$dgp$grid$betaS[iCond], params$dgp$grid$P[iCond])), times = params$dgp$nGroups), 0)
  
  
  # error variance
  # # here: when I remove the specific item variance and instead model the MI violation,
  # #       I should not need to adjust the error variance for the respective item, right?
  # adjustedError <- getErrorVar(loadingVec = c(params$dgp$miWeight[2], params$dgp$varFLgrid$flS[5]),
  #                              varianceVec = c(params$dgp$varFLgrid$varT[5], params$dgp$varFLgrid$varS[5]),
  #                              targetVar = params$dgp$varItems)
  # 
  # thetaSEM <- diag(c(rep(params$dgp$varFLgrid$varE[5], params$dgp$varFLgrid$P[5]),  # group 1
  #                    rep(c(params$dgp$varFLgrid$varE[5], adjustedError), 
  #                        c(params$dgp$varFLgrid$P[5]-params$dgp$miItems[1], params$dgp$miItems[1])), # group 2
  #                    rep(0, params$dgp$nOutcome))) # outcome
  
  thetaSEM <- diag(c(rep(params$dgp$grid$varE[iCond], params$dgp$grid$P[iCond] * 
                           params$dgp$nGroups), rep(0, params$dgp$nOutcome))) # outcome
  
  # introducing negative factor loadings changes correlation with all items of the factor
  # variance is reduced, too
  sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
    t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
  # sigmaSEM
  
  # set mean values for manifest variables
  # check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
  M <- c(rep(params$dgp$mS, params$dgp$grid$P[iCond]*params$dgp$nGroups), rep(0, params$dgp$nOutcome)) # without helper items
  # defineMeanVec(mItems = params$dgp$mItems, mOut = params$dgp$mOut, nItems = params$dgp$varFLgrid$P[iCond],
  #               nGroups = params$dgp$nGroups, nOutcome = params$dgp$nOutcome)
  
  
  # simulate date from multivariate normal distribution
  # if we have this block matrix structure we only need N/nGroups observations
  #   as variables for each group are simulated separately and otherwise we would
  #   end up with N * nGroups observations
  set.seed(372893)
  data_multipleGroup <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, 
                                         mean = M, sigma = sigmaSEM)
  # without group
  colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(params$dgp$grid$P[iCond]), 
                                                        function(iJ) {paste0("X", iJ)})), times = 2), 
                                           "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond])), "Y")
  
  # wide to long
  data <- data.frame(data_multipleGroup)
  data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
  data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
  data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
  
  ##### check data #####
  # # without specific item effect on outcome
  # # check variances across groups now
  # apply(data, 2, var)
  
  # check correlation of confounder/group with outcome and other variables
  data$group <- ifelse(data$group == "g1", 0, 1)
  varVec <- apply(data, 2, var)
  
  print(varVec)
  
}


# # checḱ variances of all items
# #   -> variances in all items (independent of group) are 1
# apply(data_multipleGroup, 2, var)
# 
# # get correlation between Y and other variables
# outcomeIdx <- which(colnames(data_multipleGroup) == "Y")
# cor(data_multipleGroup[,"Y"], data_multipleGroup[,-outcomeIdx]) 
# # ~ around 0.2 per item; except X9_g2 which correlates around -0.15
# 
# cor(data_multipleGroup[,"Y"], rowMeans(data_multipleGroup[,-outcomeIdx]))
# # 0.45

# # wide to long
# data <- data.frame(data_multipleGroup)
# data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
# data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
# data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
# 
# # check data 
# # # without specific item effect on outcome
# # # check variances across groups now
# # apply(data, 2, var)
# 
# # check correlation of confounder/group with outcome and other variables
# data$group <- ifelse(data$group == "g1", 0, 1)
# apply(data, 2, var)
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

################################################################################
# 5. simulate from multiple group SEM with scalar measurement invariance violation 
################################################################################
# nGroups <- 2
# miItems <- 1
# miWeight <- 1

condVec <- which(# params$dgp$grid$P == 4 &
  #params$dgp$grid$cAlpha == 0.8 &
  params$dgp$grid$miItems >= 1 &
    #params$dgp$grid$varT == 0.5 &
    params$dgp$grid$MIvio == "scalar")

iCond <- condVec[1] 
params$dgp$grid[iCond,]
for (iCond in condVec) {
  lambda_g1 <- cbind(TS = rep(params$dgp$grid$flT[iCond], params$dgp$grid$P[iCond]),
                     specificEffects = diag(params$dgp$grid$P[iCond]) * 
                       rep(c(1, params$dgp$grid$flS[iCond]), 
                           c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])))
  lambda_g2 <- cbind(TS = rep(c(params$dgp$grid$flT[iCond], params$dgp$grid$flT[iCond]), 
                              c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])),
                     specificEffects = diag(params$dgp$grid$P[iCond]) * 
                       rep(c(1, params$dgp$grid$flS[iCond]), 
                           c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])))
  
  
  lambdaSEM <- rbind(cbind(adiag(lambda_g1, lambda_g2), 
                           rep(0, params$dgp$nOutcome)),
                     c(rep(0, (params$dgp$grid$P[iCond] + 1)*params$dgp$nGroups), rep(1, params$dgp$nOutcome)))
  colnames(lambdaSEM) <- c(paste0(rep(c("T", 
                                        sapply(seq(params$dgp$grid$P[iCond]), 
                                               function(iJ) {paste0("S", iJ)})), times = 2), 
                                  "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond]+1)), "Y")
  rownames(lambdaSEM) <- c(paste0(rep(c(sapply(seq(params$dgp$grid$P[iCond]), 
                                               function(iJ) {paste0("X", iJ)})), times = 2), 
                                  "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond])), "Y")
  
  lambdaSEM
  
  nFactors <- ncol(lambdaSEM)
  
  # variances for T and S from the Jacobucci Paper
  #   variance of the outcome is unknwon but we only need residual variance of the
  #   outcome since y is an endogenuous variable in the structural model
  covSEM <- diag(nFactors) * c(
    rep(c(params$dgp$grid$varT[iCond], params$dgp$grid$varStable[iCond], params$dgp$grid$varS[iCond]), 
        c(1, params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])), # group 1
    rep(c(params$dgp$grid$varT[iCond], params$dgp$grid$varStable[iCond], params$dgp$grid$varS[iCond]), 
        c(1, params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])), # group 2
    params$dgp$sigmaEy) # residual variance!
  
  # MIflag = T
  # definePhi(nFactors+1,
  #           varT = params$dgp$grid$varT[iCond],
  #           varStable = params$dgp$grid$varStable[iCond],
  #           varS = params$dgp$grid$varS[iCond],
  #           nItems = params$dgp$grid$P[iCond], MIvio = MIflag,
  #           miItems = params$dgp$grid$miItems[iCond], #here: this will become part of the grid!
  #           # miItems = 0, # check warning message (;
  #           varNames = c(colnames(lambdaSEM), "Y_g2"))
  
  # path coefficients
  betaMat <- matrix(rep(0, 
                        ((params$dgp$grid$P[iCond] + params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome) + # diagonal
                          (((params$dgp$grid$P[iCond] + params$dgp$Q)*params$dgp$nGroups + 
                              params$dgp$nOutcome)*((params$dgp$grid$P[iCond] + 
                                                       params$dgp$Q)*params$dgp$nGroups + params$dgp$nOutcome - 1))), # off diagonal
                    ncol = nFactors, byrow = T)
  betaMat[nFactors,] <- c(rep(c(params$dgp$betaT, rep(params$dgp$grid$betaS[iCond], params$dgp$grid$P[iCond])), times = params$dgp$nGroups), 0)
  
  thetaSEM <- diag(c(rep(params$dgp$grid$varE[iCond], params$dgp$grid$P[iCond] * 
                           params$dgp$nGroups), rep(0, params$dgp$nOutcome))) # outcome
  
  # set mean values for manifest variables
  # check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
  mSpecific <- params$dgp$grid$b0[iCond]
  M <- c(rep(params$dgp$mItems, params$dgp$grid$P[iCond]), # group 1 
         rep(c(params$dgp$mItems, mSpecific), c(params$dgp$grid$P[iCond]-params$dgp$grid$miItems[iCond], params$dgp$grid$miItems[iCond])), # group 2 
         rep(0, params$dgp$nOutcome)) # without helper items
  
  # defineMeanVec(mItems = params$dgp$mItems, mOut = params$dgp$mOut, nItems = params$dgp$varFLgrid$P[iCond],
  #               nGroups = params$dgp$nGroups, nOutcome = params$dgp$nOutcome, mXg2 = params$dgp$varFLgrid$b0[iCond], 
  #               miItems = params$dgp$miItems[1])
  
  # introducing negative factor loadings changes correlation with all items of the factor
  # variance is reduced, too
  sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
    t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
  sigmaSEM
  
  # simulate date from multivariate normal distribution
  # if we have this block matrix structure we only need N/nGroups observations
  #   as variables for each group are simulated separately and otherwise we would
  #   end up with N * nGroups observations
  set.seed(2892)
  data_multipleGroup <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, 
                                         mean = M, sigma = sigmaSEM)
  # without group
  colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(params$dgp$grid$P[iCond]), 
                                                        function(iJ) {paste0("X", iJ)})), times = 2), 
                                           "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$grid$P[iCond])), "Y")
  
  # wide to long
  data <- data.frame(data_multipleGroup)
  data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
  data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
  data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
  
  # check data
  # without specific item effect on outcome
  # check variances across groups now
  # -> across group variances of all items are indeed 1
  # apply(data, 2, var)
  
  # check correlation of confounder/group with outcome and other variables
  data$group <- ifelse(data$group == "g1", 0, 1)
  checkVarVec <- apply(data, 2, var)
  print(checkVarVec)  
}


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

# # wide to long
# data <- data.frame(data_multipleGroup)
# data <- tidyr::pivot_longer(data, cols = !Y, names_to = "item", values_to = "value")
# data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
# data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
# 
# # check data
# # without specific item effect on outcome
# # check variances across groups now
# # -> across group variances of all items are indeed 1
# # apply(data, 2, var)
# 
# # check correlation of confounder/group with outcome and other variables
# data$group <- ifelse(data$group == "g1", 0, 1)
# apply(data, 2, var)
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

################################################################################
##### correlate confounding group with outcome #####
################################################################################
# in every MI violating condition, group is not correlated with the outcome
# correlation between group and items:
#   - configural and metric MI violations do not create correlations between items and group
#   - scalar MI violation creates correlation between MI violating item and the group

# a) simulate from SEM without any MI violation
# b) simulate from multiple group SEM with group--outcome correlation and scalar MI violation
# c) simulate from multiple group SEM with group--outcome correlation and metric MI violation

# iCond <- 5
iCond <- 18 + 5
params$dgp$varFLgrid[iCond, ]

# factor loadings
# lambdaSEM <- defineLambda(nGroups, nItems = J, miItems = 0, miWeight = list(c(NA)), params$dgp$nOutcome)
lambdaSEM <- defineLambda(flT = params$dgp$varFLgrid$flT[iCond], 
                          flS = params$dgp$varFLgrid$flS[iCond], 
                          nGroups = params$dgp$nGroups, 
                          #nGroups = 1, 
                          nItems = params$dgp$varFLgrid$P[iCond], 
                          miItems = 0, miWeight = list(c(NA)), nOutcome= params$dgp$nOutcome)
lambdaSEM
nFactors <- ncol(lambdaSEM)

# factor variances
MIflag <- if (params$dgp$varFLgrid$MIvio[iCond] == "specific") F else T
covSEM <- definePhi(nFactors, params$dgp$varFLgrid$varT[iCond], params$dgp$varFLgrid$varS[iCond], 
                    nItems = params$dgp$varFLgrid$P[iCond], MIvio = MIflag, 
                    miItems = params$dgp$miItems[1], varNames = colnames(lambdaSEM))
covSEM

# path coefficient matrix (effect of factors on the target/criterion)
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- defineBeta(nFactors, params$dgp$betaT, rep(0, params$dgp$varFLgrid$P[iCond]), 
                      params$dgp$varFLgrid$P[iCond], colnames(lambdaSEM))

# residual variance matrix of the items 
thetaSEM <- defineTheta(varE = params$dgp$varFLgrid$varE[iCond], nItems = params$dgp$varFLgrid$P[iCond],
                        params$dgp$nGroups, params$dgp$nOutcome)
thetaSEM

# # check for one variable if error matches to meet the targetVar
# facVec <- names(which(lambdaSEM["X4_g2",] != 0))
# getErrorVar(loadingVec = lambdaSEM["X4_g2",facVec], 
#             varianceVec = diag(covSEM[facVec, facVec]),
#             targetVar = 1)
# # nVariables <- dim(lambdaSEM)[1]
# rownames(lambdaSEM)[grep("^X", rownames(lambdaSEM))]
# itemNameIdxVec <- grep("^X", rownames(lambdaSEM))
# 
# # check for all variables if error matches to meet the targetVar
# sapply(itemNameIdxVec, function(iVar) {
#   iVarName <- rownames(lambdaSEM)[iVar]  
#   facVec <- names(which(lambdaSEM[iVarName,] != 0))
#   getErrorVar(loadingVec = lambdaSEM[iVarName,facVec], 
#               varianceVec = diag(covSEM[facVec, facVec]),
#               targetVar = 1)
# })

# calculate model implied covariance matrix
# for simulating from SEM
sigmaSEM <- getSigma(lambdaSEM, betaMat, covSEM, thetaSEM)
sigmaSEM

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
Y1g1Idx <- which(colnames(sigmaSEM) == "Y1_g1")
Y1g2Idx <- which(colnames(sigmaSEM) == "Y1_g2")

# we assume balance design with 2 groups! 
varWithinY <- 0.5 * sigmaSEM[Y1g1Idx, Y1g1Idx] + 0.5 * sigmaSEM[Y1g2Idx, Y1g2Idx]
targetRho <- 0.5
meanYg2 <- getMeanYg2(targetRho, params$dgp$varG, varWithinY)

M <- c(rep(params$dgp$mS, params$dgp$varFLgrid$P[iCond]*params$dgp$nGroups), 
       rep(0, params$dgp$nOutcome), rep(meanYg2, params$dgp$nOutcome)) 
M



# simulate date from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(2893)
data_multipleGroup <- mvtnorm::rmvnorm((params$sim$nTrain + params$sim$testN)/params$dgp$nGroups, mean = M, sigma = sigmaSEM)
# without group
colnames(data_multipleGroup) <- c(paste0(rep(c(sapply(seq(params$dgp$varFLgrid$P[iCond]), 
                                                      function(iJ) {paste0("X", iJ)})), times = params$dgp$nGroups), 
                                         "_g", rep(seq_len(params$dgp$nGroups), each = params$dgp$varFLgrid$P[iCond])), "Y_g1", "Y_g2")

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
# check data 
# without specific item effect on outcome
# check variances across groups now

# check correlation of confounder/group with outcome and other variables
data$group <- ifelse(data$group == "g1", 0, 1)
apply(data, 2, var)
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




# test correlation between group and outcome depending on outcome mean differences
#     between group for all MI violation conditions

##### error calculation for multiple items with ME violation #####

##### simulate nuance factor (without groups!)

