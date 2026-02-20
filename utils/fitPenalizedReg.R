# fit penalized linear regression with cross-validation to evaluate identical performance 
#     measures for different models 
# set arguments in function to control regularization: 
#   ... ridge regression (alpha = 0)
#   ... elastic net regression 
#   ... elastic net regression with regularization restricted to items
#       thus, scores have to be included in the model and only deviations from 
#             score-based predictions are represented in items 

fitRegReg <- function(trainData, testData, foldid, regMode) {
  
  # # test it
  # trainData <- dataSubset
  # testData <- data[which(!(seq_len(dim(data)[1]) %in% sampleIdx)),]
  
  # checks
  if (length(foldid) != dim(trainData)[1]) {
    stop("Fold-ID vector and data have not the same length!")
  }
  
  if (!(regMode %in% c("ridge", "enet", "itemENET"))) {
    stop("Your specified regMode does not exist.")
  }
  
  # fit only items, only scores, fit items + scores
  # ... implement this by different vectors of variable names
  itemNames <- colnames(trainData)[stringr::str_detect(colnames(trainData), "^X[:digit:]+")]
  
  # fit every combination of predictors for regMode = {ridge, enet}
  if (regMode != "itemENET") { 
    # one predictor model does not work with glmnet
    #     for one predictor, shrinkage 
    #     glmnet: "Error: x should be a matrix with 2 or more columns"
    trainPredSets <- list(itemNames, c(itemNames, "score"))  
    
    # create penaltyFactor vector for these conditions (with only defaultValue)
    #   inside the lapply()-function that iterates over predictor sets instead of 
    #   saving a list of penaltyFactor vectors with all predictor sets
  } else { # fit only item + scores for regMode = itemENET
    trainPredSets <- list(c(itemNames, "score"))
    
    # automatically set the penalty factor that controls regularization of only items
    # penalty factor = vector of length nPredictors in data (excluding intercept; intercept never penalized)
    # default penalization: 1; unpenalized features: 0
    #   Unpenalized features are still estimated jointly, so their coefficients can 
    #   change with lambda because the penalized part of the model changes; they’re just not directly shrunk.
    penaltyFactor <- c(rep(params$fit$defaultPenalty, length(itemNames)), 
                       rep(params$fit$noPenalty,  params$dgp$Q))
    names(penaltyFactor) <- c(itemNames, "score")
  }
  
  # # test it 
  # iPredSet <- trainPredSets[[2]]
  
  predSetRes <- lapply(trainPredSets, function(iPredSet) {
    # here: fix automatically generating model syntax for caret::train
    modelSyntax <- paste0("Y ~ (", paste(iPredSet, sep = "", collapse = "+"), ")") # all items
    # print(modelSyntax)
    
    # here: produce list with penalty factor = default penalty with respective predictor set length
    if (regMode != "itemENET") { # for itemENET penaltyFactor vector already exists
      penaltyFactor <- rep(params$fit$defaultPenalty, length(iPredSet))
      names(penaltyFactor) <- c(iPredSet) # match names to trainPredSet
    }
    
    # tuning grid depending on regMode
    if (regMode == "ridge") { # fix alpha to zero to fit ridge regression
      tmpTuneGrid <- expand.grid(.alpha = 0,
                                 .lambda = params$fit$lambda)
    } else { # cv over alpha and lambda
      tmpTuneGrid <- expand.grid(.alpha = params$fit$alpha,
                                 .lambda = params$fit$lambda)
    }
    
    # fixed fold IDs in caret
    # arguments in trainControl: 
    # ... index = a list with elements for each resampling iteration. Each list element 
    #     is the sample rows used for training at that iteration.
    # ... indexOut = a list (the same length as index) that dictates which sample 
    #     are held-out for each resample. If NULL, then the unique set 
    #     of samples not contained in index is used.
    # caret expects TRAINING indices in index, TEST indices in indexOut
    indexTrainCV <- setNames(lapply(seq_len(params$fit$nfolds), function(iFold) {
      which(foldid != iFold)
      }), paste0("Fold", seq_len(params$fit$nfolds)))
    indexOutTestCV <- setNames(lapply(seq_len(params$fit$nfolds), function(iFold) {
      which(foldid == iFold)
      }), paste0("Fold", seq_len(params$fit$nfolds)))
    
    # cross validation for training model
    # here: automatically supply model syntax! How do I need to write the model syntax for caret?
    fit <- caret::train(as.formula(modelSyntax),
                        trainData,
                        method = "glmnet", # elastic net regression
                        penalty.factor = penaltyFactor, # penalty factor
                        metric = "RMSE",
                        # standardize and center due to different scale the score lives on!
                        #   otherwise regularization would affect items and score differently 
                        # standardization is only applied to predictors; not to the outcome! 
                        #   thus, we do not get standardized coefficienst automatically!
                        preProcess = c("center", "scale"), 
                        #model training using ten-fold cross-validation
                        trControl = trainControl(method = "cv", 
                                                 number = params$fit$nfolds,
                                                 # fixed fold IDs for training and testing
                                                 index = indexTrainCV,
                                                 indexOut = indexOutTestCV,
                                                 # "best" = min criterion; oneSE rule
                                                 selectionFunction = params$fit$caretCrit,
                                                 search = "grid"), 
                        tuneGrid = tmpTuneGrid, #setting the tuning grid
                        verbose = FALSE)
    
    # R² across all test sets in the CV procedure (for the final hyperparameters)
    #   those are the same values as the resamples()-function gives back for the models
    #   therefore this is the result, on which the limits of the DDM should be evaluated
    performCVtest <- matrix(unlist(fit$results[rownames(fit$results) == rownames(fit$bestTune), 
                                       c("RMSE", "Rsquared", "MAE")]), , ncol = 3, nrow = 1)
    
    # standardize and center full training and full test for refitting model on 
    #     full data and test performance!
    fullTrain <- cbind(trainData[,"Y"], 
                       apply(trainData[, which(colnames(trainData) %in% iPredSet)], 
                             2, scale, center = T, scale = T))
    fullTest <- cbind(testData[,"Y"], 
                      apply(testData[, which(colnames(testData) %in% iPredSet)], 
                            2, scale, center = T, scale = T)) 
    
    # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
    predTrain <- predict(fit, fullTrain)
    performTrain <- matrix(evalPerformance(predTrain, fullTrain$Y), ncol = 3, nrow = 1)
    
    # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
    predTest <- predict(fit, fullTest)  
    performTest <- matrix(evalPerformance(predTest, fullTest$Y), ncol = 3, nrow = 1)
    
    # here: calculate prediction bias and prediction variance
    predBias <- NA_real_
    predVar <- NA_real_
    
    # build results object which saves all coefficients and standard errors of the estimated coefficients
    #   -> this should have identical dimensions across ...
    #       ... all simulated conditions (P = {4, 8}) 
    #       ... and all datasets/feature sets used for fitting the model
    # coefMat <- matrix(rep(NA_real_, (1 + max(params$dgp$P) + max(params$dgp$Q)) * 2), nrow = 2) # with SE
    coefMat <- matrix(rep(NA_real_, (1 + max(params$dgp$P) + max(params$dgp$Q))), nrow = 1) # without SE
    colnames(coefMat) <- c("(Intercept)", 
                           paste0("X", seq_len(max(params$dgp$P))), 
                           rep("score", max(params$dgp$Q)))
    
    # extract coefficients of the final model
    betaVec <- predict(fit$finalModel, type = "coefficients", s = fit$bestTune$lambda) 
    
    # fill in regression coefficients
    # due to perfect collinearity of some predictor variables, regression coefficient 
    #   might become NA in some simulated conditions; therefore check is.na()
    nameIDX <- match(rownames(betaVec), colnames(coefMat))
    keep <- !is.na(nameIDX)
    coefMat[1, nameIDX[keep]] <- betaVec[keep]
    
    # do not save SEs of the ridge or elastic net regression coefficient estimates
    # ... for ridge regression approximate SEs with fixed λ would be possible to calculate
    # ... for elastic net SEs only via bootstrapping which would arguably take too long
    # -> SEs not important enough for extra computational effort
    
    # object with results to return to predSetRes 
    list(
      # {RMSE, Rsquared, MAE} x {CVtest, train, HOtest}
      performVec = cbind(performCVtest, performTrain, performTest), 
      # prediction bias and prediction variance 
      # predError = matrix(c(predBias, predVar), ncol = 2, nrow = 1),
      # predictor choice by extracting regression coefficients (for all possible coefficients)
      # matrix with regression coefficients in row 1 and (no) SEs in row 2 
      regCoef = coefMat,
      # save hyperparameter values to check if cross validation chose ridge 
      hyperParams = matrix(c(alpha = fit$bestTune$alpha, lambda = fit$bestTune$lambda), 
                           ncol = 2, nrow = 1)) 

  })
  
  # build results object to return
  # resPenReg <- vector(mode = "list", 4)
  # names(resPenReg) <- c("performVec", "predError", "regCoef", "hyperParams")
  resPenReg <- vector(mode = "list", 3)
  names(resPenReg) <- c("performVec", "regCoef", "hyperParams")
  
  # object with results to return to fitData Script 
  resPenReg[["performVec"]] <- do.call(rbind, lapply(predSetRes, `[[`, "performVec"))
  # name columns of performance Matrix
  colnames(resPenReg[["performVec"]]) <- c(unite(expand.grid(c("RMSE", "Rsquared", "MAE"), 
                                                          c("CVtest", "train", "HOtest")), 
                                              "names"))[["names"]]
  
  # resPenReg[["predError"]] <- do.call(rbind, lapply(predSetRes, `[[`, "predError"))
  resPenReg[["regCoef"]] <- do.call(rbind, lapply(predSetRes, `[[`, "regCoef"))
  resPenReg[["hyperParams"]] <- do.call(rbind, lapply(predSetRes, `[[`, "hyperParams"))
  
  return(resPenReg)
  
}