# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R") # general functions for simulation management
source("utils/simTools.R") # general functions for data simulation

# functions to generate data
source("utils/sampleData.R")

# check if all packages are available; install packages that are not available
# default is to call functions from packages explicitly (with ::) instead of loading
#   packages into search space/attaching everything
packageVec <- c("parallel", "mvtnorm", "magic", "lavaan", "semTools", "tidyr")
sapply(packageVec, function(iPackage) {
  if (!require(iPackage, character.only = TRUE)) {
    install.packages(iPackage)
  }
})

require(mvtnorm)
library(magic)
library(lavaan)
library(parallel)

timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")

# generate folder for data files
dataFolder <- "data"
createFolder(dataFolder)

# generate folder for log files
logFolder <- "log"
createFolder(logFolder)

# set number of cores to use in parallel computing
nCoresSampling <- 30 # simulate data

# grid to simulate data with mapply later
#   simulation via mapply to easily simulate subsets of parameter combinations
#   create factor structure for predictor variables to increase reliability
# grid with ...
# ... P = params$dgp$P # number of items per factor
# ... Q = params$dgp$Q # number of factors
# ... cAlpha = params$dgp$cAlpha # Cronbach's Alpha
# ... varT = 
# ... miItems = params$dgp$miItems # number of MI violating items
# ... corGO = params$dgp$rGroupOut # correlation group -- outcome
# ... R2 = 
# ... varStable =
# ... flT
# ... flS
# ... varE
# ... MIvio = params$dgp$MIvio # measurement invariance violations
# ... b0
# ... b1
# ... betaT
# ... varS
# ... betaS

grid <- params$dgp$grid
# ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
str(grid) 

# # sanity checks for data
# params$sim$nSamples <- 1

##### test it #####
# # test simulation condition that crushes SEM estimation
# condVec <- which(params$dgp$grid$P == 8 &
#                    params$dgp$grid$cAlpha == 0.6 & # 0.6
#                    params$dgp$grid$miItems == 4 &
#                    params$dgp$grid$corGO == 0 &
#                    params$dgp$grid$R2 == 0.8 &
#                    params$dgp$grid$varT == 0.25 & # 0.25
#                    params$dgp$grid$MIvio == "nuanceFactor")
# iCond <- condVec[1]
# P <- params$dgp$grid$P[iCond]
# Q <- params$dgp$grid$Q[iCond]
# cAlpha <- params$dgp$grid$cAlpha[iCond]
# varT <- params$dgp$grid$varT[iCond]
# miItems <- params$dgp$grid$miItems[iCond]
# corGO <- params$dgp$grid$corGO[iCond]
# R2 <- params$dgp$grid$R2[iCond]
# varStable <- params$dgp$grid$varStable[iCond]
# flT <- params$dgp$grid$flT[iCond]
# flS <- params$dgp$grid$flS[iCond]
# varE <- params$dgp$grid$varE[iCond]
# MIvio <- params$dgp$grid$MIvio[iCond]
# b0 <- params$dgp$grid$b0[iCond]
# b1 <- params$dgp$grid$b1[iCond]
# betaT <- params$dgp$grid$betaT[iCond]
# varS <- params$dgp$grid$varS[iCond]
# betaS <- params$dgp$grid$betaS[iCond]
# covNuanceFactor <- params$dgp$grid$covNuanceFactor[iCond]
# sampleSeed <- grid$sampleSeed[iCond]
# 
# params$dgp$grid[iCond,]

#####
options(warn = -1)

# sample data in parallel
createData <- function(P, Q, cAlpha, varT, miItems, corGO, R2, varStable, flT, flS,
                       varE, MIvio, b0, b1, betaT, varS, betaS, covNuanceFactor, sampleSeed){
  
  # checkFlag: check deviation of parameter estimates from true parameter values
  #            and get (estimated) standardized factor loadings 
  checkFlag <- FALSE # set this to true to get checkSample Data 
  
  ##### measurement model #####
  # ... always simulate nGroups = 2 in every condition but keep parameters identical! 
  
  # generate lambda matrix
  lambdaSEM <- defineLambda(flT = flT, flS = flS, nGroups = params$dgp$nGroups, 
                            nItems = P, miItems = miItems, 
                            # for MIvio %in% c("specific", "scalar") miWeight is automatically flT
                            miWeight = flT + b1, nOutcome = params$dgp$nOutcome, 
                            nFactors = Q)
  
  nLatentSEM <- ncol(lambdaSEM)
  
  # generate covariance matrix for the factors
  #   variance of the outcome is unknwon but we only need residual variance of the
  #   outcome since y is an endogenuous variable in the structural model
  # MIflag <- if (MIvio != "specific") T else F
  MIflag <- if (!(MIvio %in% c("specific", paste0("NF_", params$dgp$corNuanceFactor)))) T else F
  covSEM <- definePhi(nLatentSEM, varT = varT, varStable = varStable, varS = varS,
                      nItems = P, MIvio = MIvio, miItems = miItems, 
                      covNuanceFactor = covNuanceFactor, 
                      varNames = colnames(lambdaSEM))
  
  # generate residual variance matrix
  # I remove the specific item variance and instead model the MI violation with exactly
  #   the same variance as for the specific item variance
  #   thus, I do not need to adjust the residual variance of the MI violating item, right?
  thetaSEM <- defineTheta(varE = varE, nItems = P, 
                          nGroups = params$dgp$nGroups, nOutcome = params$dgp$nOutcome)
  
  ##### structure model #####
  # generate beta matrix
  betaSEM <- defineBeta(nLatentSEM, betaT = betaT, betaS = betaS, miItems = miItems, 
                        nItems = P, varNames = colnames(lambdaSEM), nFactors = 1)
  
  ##### complete SEM #####
  # calculate model implied covariance matrix for ...
  #     simulating from SEM by sampling from multivariate normal distribution 
  sigmaSEM <- getSigma(lambdaSEM, betaSEM, covSEM, thetaSEM)
  
  # generate mean vector for data sampling
  # simulate group-outcome correlations by mean differences in the outcome
  #   -> calculate difference in outcome for target correlation assuming balanced design with 2 groups! 
  if (corGO != 0) {
    Y1g1Idx <- which(colnames(sigmaSEM) == "Y1_g1")
    Y1g2Idx <- which(colnames(sigmaSEM) == "Y1_g2")
    varWithinY <- 0.5 * sigmaSEM[Y1g1Idx, Y1g1Idx] + 0.5 * sigmaSEM[Y1g2Idx, Y1g2Idx]
    meanYg2 <- getMeanYg2(targetRho = corGO, params$dgp$varG, varWithinY)
  } else {
    # depending on correlation between group and outcome; if cor = 0: mYg2 = mOut
    meanYg2 <- params$dgp$mOut
  }
  
  M <- defineMeanVec(mItems = params$dgp$mItems, mOut = params$dgp$mOut, 
                     nItems = P, nGroups = params$dgp$nGroups, 
                     nOutcome = params$dgp$nOutcome, mXg2 = b0, mYg2 = meanYg2,
                     miItems = miItems)
  
  ##############################################################################
  # start random sampling (in parallel)
  ##############################################################################
  environment(sampleData) <- environment()  

  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "sampleDataStatus_", timeStampFolder, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(sampleSeed)
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  # function that does random sampling
  # only sample data for the largest N! draw random subsets for smaller Ns when fitting the data
  sampleData() # run function to actually create dataset
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
}

# actually sample data (or save checks - checkFlag = T)
out <- do.call(mapply, c(FUN = createData, grid))

options(warn = 0)
