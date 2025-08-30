################################################################################
# fit data
################################################################################

# Model fitting
library(parallel) # fit models in parallel across samples
library(glmnet)   # ENET
library(caret)    # GBM

# load parameters & custom functions 
source("utils/setParam.R") # parameter values
# source("utils/simTools.R") 
source("utils/utilityFUNS.R")
source("utils/fitTools.R")

# functions to fit different models
source("utils/fitOLS.R")
source("utils/fitPenalizedReg.R")

# functions to calculate and save outcome measures from models
# source("utils/saveResults.R")

# build folder structure
dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 30

grid <- params$dgp$grid

# only load biggest dataset! draw random subsamples and fit smaller subsets
# this is done to avoid excessive memory requirements to save the data!
# N = {50, 100, 250, 500, 1000}

# prediction based on linear regression and based on elastic net regression
#   ... prediction based on items and based on composite
#   ... prediction based on predictor set with group, without group, with group predicted from items

results <- lapply(seq_len(nrow(grid)), function(iSim) {
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "fitDataStatus",
                                     timeStamp, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(grid[iSim, "sampleSeed"])
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  # file name for all samples in one big rda file or folder name for single sample rda files
  conditionFolder <- paste0("simData_MIvio", grid$MIvio[iSim], "_miItems", grid$miItems[iSim], 
                     "_P", grid$P[iSim], "_Q", grid$Q[iSim], "_cAlpha", grid$cAlpha[iSim], 
                     "_varT", grid$varT[iSim], "_corGO", grid$corGO[iSim], "_R2", grid$R2[iSim])
  
  
  tStart <- Sys.time() # time monitoring
  
  # iterate over different training samples in parallel
  estRes <- parLapply(cl, seq_len(params$sim$nTrain), function(iSample) {
    
    # monitor progress of the simulation study (irrespective of model)
    # write progress into log file to keep track of simulation progress
    if (iSample == 1) {
      print(paste0(conditionFolder,
                   "; first Sample at ", format(Sys.time(), "%d%m%y_%H%M%S")))  
    }
    
    # load training and holdout test data for single sample rda files
    # load training sample
    load(paste0(dataFolder, "/", conditionFolder, "/sample_", iSample, ".rda"))
    
    # load test data
    # setParam$dgp$nTest == 1
    load(paste0(dataFolder, "/", conditionFolder, "/sample_test1.rda"))
    
    # lapply over different conditions of sample sizes
    estN <- lapply(seq_along(params$dgp$trainN), function(iN) { 
      
      # random sampling from training sample (reproducible?)
      sampleIdx <- sample(seq_len(max(params$dgp$trainN)), size = params$dgp$trainN[iN], replace = FALSE)
      dataSubset <- data[sampleIdx, ]
      
      # # checks
      # length(sampleIdx)
      # length(unique(sampleIdx))
      
      
      # set foldid 
      # ... for every N separately
      # ... but identical foldids for every model within the same condition
      foldid <- sample(rep(seq_len(params$fit$nfolds), length.out =  params$dgp$trainN[iN]), replace = FALSE)
      
      # run fit functions with different settings 
      # ... fit data on training sample
      # ... fit different datasets within every fit function {onlyItems, onlyScores, Items + Score}
      # ... test fit model on holdout data within fit function
      
      # linear regression 
      fit_LM <- fitOLS(trainData = dataSubset, testData = ?, foldid = foldid) 
      
      # ridge regression
      fit_Ridge <- fitRegReg(trainData = dataSubset, testData = ?, 
                             foldid = foldid, regMode = "ridge") 
      
      # elastic net regression
      fit_ENET <- fitRegReg(trainData = dataSubset, testData = ?, 
                            foldid = foldid, regMode = "enet") 
      
      # elastic net regression only item penalty
      fit_ENETitems <- fitRegReg(trainData = dataSubset, testData = ?, 
                                 foldid = foldid, regMode = "itemENET") 
      
    })
    
    # combine results for every subsampled N 
    
    
  }) # end of parallel fitting for samples
  tEnd <- Sys.time() # time monitoring
  difftime(tEnd, tStart)
})

######
#  true prediction by using Construct and specific item effect both as predictors
#    ->  this is not feasible anymore; specific item effects are simulated as latent 
#        factors that are not identified and therefore cannot be estimated within the model
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
