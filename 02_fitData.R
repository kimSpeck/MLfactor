################################################################################
# fit data
################################################################################
# Model fitting
library(parallel) # fit models in parallel across samples
library(glmnet)   # ENET
library(caret)    # GBM
library(tidyr)    # GBM

# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")
source("utils/fitTools.R")

# functions to fit different models
source("utils/fitOLS.R")
source("utils/fitPenalizedReg.R")

# build folder structure
dataFolder <- "data"

resFolder <- "results"
createFolder(resFolder)

logFolder <- "log"
createFolder(logFolder)

timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
# nCoresSampling <- detectCores() - 1 
nCoresSampling <- 20

# grid <- params$dgp$grid
# grid <- params$dgp$grid[which(params$dgp$grid$MIvio %in% paste0("NF_", params$dgp$corNuanceFactor)), ]
grid <- params$dgp$grid[which(params$dgp$grid$MIvio %in% paste0("scalarReal", c(0.4, 0.6))), ]

# only load biggest dataset! draw random subsamples and fit smaller subsets
# this is done to avoid excessive memory requirements to save the data!

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
  
  # # test it
  # iSim <- 73
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
    trainSample <- loadRData(paste0(dataFolder, "/", conditionFolder, "/sample_", iSample, ".rda"))
    
    # load test data
    testSample <- loadRData(paste0(dataFolder, "/", conditionFolder, "/sample_test1.rda"))
    
    # lapply over different conditions of sample sizes
    estN <- lapply(seq_along(params$dgp$trainN), function(iN) { 
      
      # random sampling from training sample (reproducible?)
      # this needs to be done separately for the two groups!
      #   otherwise the manipulation of the MI violation is not correct anymore
      # 1. generate two sets of indices for group = 0 and for group = 1
      # 2. sample the params$dgp$trainN[iN]/2 indices from both sets of indices
      # 3. generate the trainData subset of size params$dgp$trainN[iN] with 50%
      #     group 0 and group 1 respectively by using the union of both index sets
      gZeroIDX <- which(trainSample$group == 0)
      gOneIDX <- which(trainSample$group == 1)
      sampleIdx <- c(sample(gZeroIDX, size = (params$dgp$trainN[iN]/params$dgp$nGroups), 
                            replace = FALSE),
                     sample(gOneIDX, size = (params$dgp$trainN[iN]/params$dgp$nGroups), 
                            replace = FALSE))
      trainSubset <- trainSample[sampleIdx, ]
      
      # # checks
      # length(sampleIdx)
      # length(unique(sampleIdx))
      
      # set foldid 
      # ... for every N separately
      # ... but identical foldids for every model within the same condition
      foldid <- sample(rep(seq_len(params$fit$nfolds), length.out =  params$dgp$trainN[iN]), 
                       replace = FALSE)
      
      # run fit functions with different settings 
      # ... fit data on training sample
      # ... fit different datasets within every fit function {onlyItems, onlyScores, Items + Score}
      # ... test fit model on holdout data within fit function
      modelList <- vector(mode = "list", length = 8)
      names(modelList) <- c("fit_LM", "fit_Ridge", "fit_ENET", "fit_ENETitems")
      
      # linear regression 
      modelList[["fit_LM"]] <- fitOLS(trainData = trainSubset, 
                       testData = testSample,  
                       foldid = foldid, group = params$fit$fitGroup) 
    
      # ridge regression
      modelList[["fit_Ridge"]] <- fitRegReg(trainData = trainSubset,
                             testData = testSample,
                             foldid = foldid, regMode = "ridge", group = params$fit$fitGroup) 
      
      # elastic net regression
      modelList[["fit_ENET"]] <- fitRegReg(trainData = trainSubset, 
                            testData = testSample,
                            foldid = foldid, regMode = "enet", group = params$fit$fitGroup) 
     
      # elastic net regression only item penalty
      modelList[["fit_ENETitems"]] <- fitRegReg(trainData = trainSubset, 
                                 testData = testSample,
                                 foldid = foldid, regMode = "itemENET", group = params$fit$fitGroup) 
      
      if (params$fit$fitGroup) {
        list(performVec = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                                model = rep(c("LM", "Ridge", "ENET", "itemENET"), 
                                            c(6, 4, 4, 2)),
                                predSet = c("items", "score", "comb", 
                                            "items_group", "score_group", "comb_group", 
                                            "items", "comb", 
                                            "items_group", "comb_group", 
                                            "items", "comb",
                                            "items_group", "comb_group",
                                            "comb", "comb_group"),
                                do.call(rbind, lapply(modelList, `[[`, "performVec"))),
             # predError = do.call(rbind, lapply(modelList, `[[`, "predError")),
             regCoef = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                             model = rep(c("LM", "Ridge", "ENET", "itemENET"), 
                                         c(6, 4, 4, 2)),
                             predSet = c("items", "score", "comb", 
                                         "items_group", "score_group", "comb_group", 
                                         "items", "comb", 
                                         "items_group", "comb_group", 
                                         "items", "comb",
                                         "items_group", "comb_group",
                                         "comb", "comb_group"),
                             do.call(rbind, lapply(modelList, `[[`, "regCoef"))),
             hyperParams = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                                 model = rep(c("Ridge", "ENET", "itemENET"), 
                                             c(4, 4, 2)),
                                 predSet = c("items", "comb", 
                                             "items_group", "comb_group", 
                                             "items", "comb",
                                             "items_group", "comb_group",
                                             "comb", "comb_group"),
                                 do.call(rbind, lapply(modelList, `[[`, "hyperParams"))))
      } else { # result object without predictor sets including group as a covariate
        list(performVec = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                                model = rep(c("LM", "Ridge", "ENET", "itemENET"), c(3, 2, 2, 1)),
                                predSet = rep(c("items", "score", "comb", "items", "comb", "items", "comb", "comb")),
                                do.call(rbind, lapply(modelList, `[[`, "performVec"))),
             # predError = do.call(rbind, lapply(modelList, `[[`, "predError")),
             regCoef = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                             model = rep(c("LM", "Ridge", "ENET", "itemENET"), c(3, 2, 2, 1)),
                             predSet = rep(c("items", "score", "comb", "items", "comb", "items", "comb", "comb")),
                             do.call(rbind, lapply(modelList, `[[`, "regCoef"))),
             hyperParams = cbind(sample = iSample, N = rep(params$dgp$trainN[iN]), 
                                 model = rep(c("Ridge", "ENET", "itemENET"), c(2, 2, 1)),
                                 predSet = c("items", "comb", "items", "comb", "comb"),
                                 do.call(rbind, lapply(modelList, `[[`, "hyperParams"))))
      }
      
      
    })
    
    # combine results for every subsampled N 
    # save condition information!
    # save sample number! 
    resAllN <- vector(mode = "list", 3)
    names(resAllN) <- c("performVec", "regCoef", "hyperParams")
    
    resAllN[["performVec"]] <- do.call(rbind, lapply(estN, `[[`, "performVec"))
    # do.call(rbind, lapply(estN, `[[`, "predError"))
    resAllN[["regCoef"]] <- do.call(rbind, lapply(estN, `[[`, "regCoef"))
    resAllN[["hyperParams"]] <- do.call(rbind, lapply(estN, `[[`, "hyperParams"))
    
    resAllN # this list with performVec for all models and all Ns is saved for every sample 
    
  }) # end of parallel fitting for samples
  
  tEnd <- Sys.time() # time monitoring
  difftime(tEnd, tStart)
  
  # manage data
  # save separate each simulated condition in rds files  
  #   -> pro: skip fitted condGrid conditions if algorithm does not run till the end...
  resFileName <- paste0(resFolder, "/", "resFile_MIvio", grid$MIvio[iSim], 
                        "_miItems", grid$miItems[iSim], "_P", grid$P[iSim], 
                        "_Q", grid$Q[iSim], "_cAlpha", grid$cAlpha[iSim], 
                        "_varT", grid$varT[iSim], "_corGO", grid$corGO[iSim], 
                        "_R2", grid$R2[iSim], ".rds")
  saveRDS(estRes, file = resFileName)
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
  gc()
  
})
