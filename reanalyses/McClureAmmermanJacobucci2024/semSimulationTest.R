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
#   ->  this is due to the rather formative conception of the true score! 
#       the true score enters y multiple time in the original Jacobucci simulation
#   -> now that we simulate a reflective model in the original simulation the variance
#      matches with the sem formulation of the data 

setwd("../reanalysis/McClureAmmermanJacobucci2023/Simulation/")
source('simTools.R')
source('setParameters.R')
source('Linear-Regression-Bias-Variance_functions.R')

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
betaS.cond <- 4
beta_S <- setBetaS(betaS.cond, params$sigma[measVar.cond,"sigmaS"], J, params)
betaS <- params$betaS[betaS.cond]/params$sigma[measVar.cond,"sigmaS"]
params$nOutcome <- 1

N = params$holdout + N.train
W = popWeight(params$M,J) #creates a weight matrix that partitions J in M~equal unit-weighted composites

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

##### simulate from SEM #####
# factor loading matrix with outcome; but variance of the outcome is unknwon
#   therefore: factor loading matrix without outcome below
lambdaSEM <- cbind(trueScore = c(rep(1, J), rep(0, params$nOutcome)),
                   specificEffects = t(cbind(diag(J), 0)),
                   Y = c(rep(0, J), rep(1, params$nOutcome)))
nFactors <- ncol(lambdaSEM)
colnames(lambdaSEM) <- c("T", 
                         sapply(seq(J), function(iJ) {paste0("S", iJ)}),
                         "Y")

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

# # these is the variance of the composite with the true score entering multiple times
# # varFormComp <- 40.5873
# varFormComp <- 45
# covSEM <- diag(nFactors) * c(varFormComp, diag(Sigma_S), params$sigmaErrY)

# path coefficient matrix
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- matrix(rep(0, 
                      (J + params$M + params$nOutcome) + # diagonal
                        ((J + params$M + params$nOutcome)*(J + params$M + params$nOutcome - 1))), # off diagonal
                  ncol = nFactors, byrow = T)
betaMat[nFactors,] <- c(betaT, beta_S, 0)

# error - Kovarianzmatrix der Residualvariablen mit p x p (= Anzahl der manifesten Variablen)
# outcome as manifest variable but all latent! (therefore the error variance is zero?)
thetaSEM <- diag(c(diag(Sigma_e), rep(0, params$nOutcome)))

# calculate model implied covariance matrix
# for simulating from SEM
sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% covSEM %*%
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM

# set mean values for manifest variables
# check this: mean of the outcome (y) is zero now... this might be different in the Jacobucci simulation?
M <- c(rep(params$mS, J), rep(0, params$nOutcome))

# simulate date from multivariate normal distribution
set.seed(seed)
data <- mvtnorm::rmvnorm(N.train + params$holdout, mean = M, sigma = sigmaSEM)
colnames(data) <- c(paste0("X", seq_len(J)), "Y")

head(data)
apply(data, 2, var) # this looks fine
######
# check SEM data
######
m0 <- '# measurement model
       T =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9

       # all latent representation of manifest outcome
       fY =~ 1 * Y
       Y ~~ 0*Y
       # structure model
       fY ~ T'

m0 <- '# measurement model
       T =~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9
       S1 =~ X1

       # all latent representation of manifest outcome
       fY =~ 1 * Y
       Y ~~ 0*Y
        
       S1 =~ 1*X1
       X1 ~~ 0.25*X1

       # structure model
       fY ~ T'


# fit <- sem(mTrue, data = dataList[["lambdaN0.2_error0.2_betaNY0.6_betaFY0.6"]])
apply(data, 2, var)
fit <- lavaan::sem(m0, data = data)
fit@Model@GLIST

#dataSpecific <- data 
apply(dataSpecific, 2, var)
fit <- lavaan::sem(m0, data = dataSpecific)
fit@Model@GLIST
######
# calculate the composite score
#     for the composite we do not use true W but unit weights for all items in every group
Xcomp = data[, paste0("X", seq_len(J))] %*% W # this is the composite (out of all manifest variables); matrix of composites scores weighted by W (NxM)

data <- data.frame(data)
# Creating Train-Test sets
# X = (manifest) item values = T + S + e 
itms.df = data.frame(data$Y, data[, paste0("X", seq_len(J))])
itms.train=itms.df[1:N.train,]
itms.test=itms.df[(N.train+1):N,]

# ("manifest") composite values = (1*T1 + 1*T2 + ... + 1*Tj) + (S1 + S2 + ... + Sj) + 
comp.df = data.frame(data$Y,Xcomp)
comp.train=comp.df[1:N.train,]
comp.test=comp.df[(N.train+1):N,]

# combine items + construct
comb.df = data.frame(data$Y, data[, paste0("X", seq_len(J))], Xcomp)
comb.train = comb.df[1:N.train,]
comb.test = comb.df[(N.train+1):N,]

# true prediction by using Construct and specific item effect both as predictors
#   ->  this is not feasible anymore; specific item effects are simulated as latent 
#       factors that are not identified and therefore cannot be estimated within the model

# Regressions
lm.items = lm(data.Y~.,itms.train)
lm.comp = lm(data.Y~.,comp.train)
lm.comb = lm(data.Y~.,comb.train)

# observed coefficients excluding intercept
items.coef = coef(lm.items)[2:(J + 1)]
comp.coef = coef(lm.comp)[2:(params$M + 1)]
comb.coef = coef(lm.comb)[2:(J + params$M + 1)]

# Compute the w/in sample bias for both models at the observation level
#   ->  this is not feasible anymore; specific item effects are simulated as latent 
#       factors that are not identified and therefore cannot be estimated within the model

#Calculate Holdout Prediction
itms.pred = predict(lm.items,itms.test)
comp.pred = predict(lm.comp,comp.test)
comb.pred = predict(lm.comb,comb.test)

#Preparing Output
mat = rep(NA, 28)

#Condition Information
mat[1] = seed; mat[2] = N.train; mat[3] = J #note that N -holdout is returned not total N
mat[4] = betaS;mat[5] = betaT; mat[6] = measVar.cond #beta_S, beta_T, and measVar conditions
mat[7] = nGroups; mat[8] = miItems; mat[9] = miWeight; mat[10] = mSpecific
#Bias
mat[11] = NA
mat[12] = NA
mat[13] = NA

#Variance
mat[14] = NA
mat[15] = NA
mat[16] = NA

#within samp error
mat[17] = NA
mat[18] = NA
mat[19] = NA

#Holdout RMSE
mat[20] = calcRMSE(itms.pred,itms.test[,"data.Y"])
mat[21] = calcRMSE(comp.pred,comp.test[,"data.Y"])
mat[22] = calcRMSE(comb.pred,comb.test[,"data.Y"])
mat[23] = NA

# test R2
mat[24] <- evalPerformance(itms.pred,itms.test[,"data.Y"])["Rsquared"]
mat[25] <- evalPerformance(comp.pred,comp.test[,"data.Y"])["Rsquared"]
mat[26] <- evalPerformance(comb.pred,comb.test[,"data.Y"])["Rsquared"]
mat[27] <- NA

names(mat) <-  c("Seed","N","J","beta_S","beta_T","measVar",
                 "nGroups", "miItems", "miWeight", "mSpecific",
                 "item.bias","comp.bias","true.bias",
                 "item.var","comp.var","true.var",
                 "item.err","comp.err","true.err",
                 "item.rmse","comp.rmse", "comb.rmse", "true.rmse",
                 "item.R2", "comp.R2", "comb.R2", "true.R2",
                 "sim.R2")
# SEM is not identified; thus, we cannot estimate SEM to check simulation
# To Do: simulate MI violations in SEM framework

m0 <- '# measurement model
       T =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7 + 1*X8 + 1*X9
       S1 =~ 1*X1
       S2 =~ 1*X2
       S3 =~ 1*X3
       S4 =~ 1*X4
       S5 =~ 1*X5
       S6 =~ 1*X6
       S7 =~ 1*X7
       S8 =~ 1*X8
       S9 =~ 1*X9

       # uncorrelated factor and nuances
       T ~~ 0*S1
       T ~~ 0*S2
       T ~~ 0*S3
       T ~~ 0*S4
       T ~~ 0*S5
       T ~~ 0*S6
       T ~~ 0*S7
       T ~~ 0*S8
       T ~~ 0*S9
       S1 ~~ 0*S2
       S1 ~~ 0*S3
       S1 ~~ 0*S4
       S1 ~~ 0*S5
       S1 ~~ 0*S6
       S1 ~~ 0*S7
       S1 ~~ 0*S8
       S1 ~~ 0*S9
       S2 ~~ 0*S3
       S2 ~~ 0*S4
       S2 ~~ 0*S5
       S2 ~~ 0*S6
       S2 ~~ 0*S7
       S2 ~~ 0*S8
       S2 ~~ 0*S9
       S3 ~~ 0*S4
       S3 ~~ 0*S5
       S3 ~~ 0*S6
       S3 ~~ 0*S7
       S3 ~~ 0*S8
       S3 ~~ 0*S9
       S4 ~~ 0*S5
       S4 ~~ 0*S6
       S4 ~~ 0*S7
       S4 ~~ 0*S8
       S4 ~~ 0*S9
       S5 ~~ 0*S6
       S5 ~~ 0*S7
       S5 ~~ 0*S8
       S5 ~~ 0*S9
       S6 ~~ 0*S7
       S6 ~~ 0*S8
       S6 ~~ 0*S9
       S7 ~~ 0*S8
       S7 ~~ 0*S9
       S8 ~~ 0*S9

       # all latent representation of manifest outcome
       fY =~ 1 * Y
       Y ~~ 0*Y
       
       # structure model
       fY ~ T + S1'

m0 <- '# measurement model
       T =~ 1*X1 + 1*X2 + 1*X3 + 1*X4 + 1*X5 + 1*X6 + 1*X7 + 1*X8 + 1*X9
       S1 =~ 1*X1
       S2 =~ 1*X2
       S3 =~ 1*X3
       S4 =~ 1*X4
       S5 =~ 1*X5
       S6 =~ 1*X6
       S7 =~ 1*X7
       S8 =~ 1*X8
       S9 =~ 1*X9

       # all latent representation of manifest outcome
       fY =~ 1 * Y
       Y ~~ 0*Y

       # structure model
       fY ~ T + S1'