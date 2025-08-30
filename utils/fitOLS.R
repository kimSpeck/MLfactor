# fit linear regression with cross-validation to evaluate identical performance 
#     measures for different models and avoid overfitting

fitOLS <- function(trainData, testData, foldid) {
  
  # # test it
  # trainData <- dataSubset
  # testData <- data[which(!(seq_len(dim(data)[1]) %in% sampleIdx)),]
  
  # checks
  if (length(foldid) != dim(trainData)[1]) {
    stop("Fold-ID vector and data have not the same length!")
  }
  
  if (Q > 1) {
    stop("Result object is not implemented for multiple factors.")
  }
  
  # fit only items, only scores, fit items + scores
  # ... implement this by different vectors of variable names
  itemNames <- paste0("X", seq_len(P))
  trainPredSets <- list(itemNames, "score", c(itemNames, "score"))
  
  # # test it 
  # iPredSet <- trainPredSets[[1]]
  
  # iterate over vectors of variable names to fit items, scores and both combined
  predSetRes <- lapply(trainPredSets, function(iPredSet) {
    lm_modelSyntax <- paste0("Y ~ (1 +", paste(iPredSet, sep = "", collapse = "+"), ")") # all items
    # print(lm_modelSyntax)
    
    # result objects for CV
    yhatCV <- rep(NA_real_,  dim(trainData)[1])
    performPerFold <- vector("list", params$fit$nfolds)

    # cross-validation with set foldids
    for (iFold in seq_len(params$fit$nfolds)) {
      # separate data into (nFold - 1) folds as training data and the test fold 
      testIDX  <- which(foldid == iFold) # indizes in this specific fold
      trainIDX <- setdiff(seq_len(dim(trainData)[1]), testIDX) # all other indices
      
      # fit regression on the (nFold - 1) folds training data
      tmpFitCV <- lm(lm_modelSyntax, data = trainData[trainIDX, ])
      
      # use fitted model to predict outcome in test fold
      pred <- as.numeric(predict(tmpFitCV, newdata = trainData[testIDX, ]))
      yTest  <- trainData$Y[testIDX] # true outcome values in test fold
      
      # performance metrics on test data for every fold separately
      RMSE <- sqrt(mean((yTest - pred)^2))
      SST <- sum((yTest - mean(yTest))^2)
      SSE <- sum((yTest - pred)^2)
      R2 <- if (SST > 0) 1 - (SSE/SST) else NA_real_
      MAE <- mean(abs(pred - yTest))
      
      # # this is overfitting in case of cross-validation because it essentially does rescaling!
      # # this is why I do not use evalPerformance function here!
      # R2cor <- cor(pred, yTest, use = "pairwise.complete.obs")**2
      
      yhatCV[testIDX] <- pred
      performPerFold[[iFold]] <- data.frame(
        fold = iFold, # nTrain = length(trainIDX), nTest = length(testIDX),
        RMSE = RMSE, Rsquared = R2, MAE = MAE
      )
    }
    
    performPerFold <- do.call(rbind, performPerFold)
    performCVtest <- matrix(colMeans(performPerFold[, c("RMSE", "Rsquared", "MAE")]), ncol = 3, nrow = 1)
    
    # fit regression on full training data (no tuning parameters, thus new model)
    # attention: one predictor column is an exact linear combination of others (or is all zeros). 
    # In that case the corresponding regression coefficient is not identifiable, so 
    # lm drops it and reports it as NA with the message “(1 not defined because of singularities)
    fit <- lm(lm_modelSyntax, data = trainData)
    
    # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model training
    predTrain <- predict(fit, trainData)
    performTrain <- matrix(evalPerformance(predTrain, trainData$Y), ncol = 3, nrow = 1)
    
    # # demonstrate real quick that R2 and R2cor are identical now
    # SST <- sum((trainData$Y - mean(trainData$Y))^2)
    # SSE <- sum((trainData$Y - predTrain)^2)
    # (R2 <- if (sst > 0) 1 - (SSE/SST) else NA_real_)
    # (R2cor <- cor(predTrain, trainData$Y, use = "pairwise.complete.obs")**2)
    
    # get Root Mean Squared Error (RMSE), explained variance (R2), and Mean Absolute Error (MAE) for model testing
    predTest <- predict(fit, testData)  
    performTest <- matrix(evalPerformance(predTest, testData$Y), ncol = 3, nrow = 1)
    
    # here: calculate prediction bias and prediction variance
    predBias <- NA_real_
    predVar <- NA_real_
    
    # build results object which saves all coefficients and standard errors of the estimated coefficients
    #   -> this should have identical dimensions across ...
    #       ... all simulated conditions (P = {4, 8}) 
    #       ... and all datasets/feature sets used for fitting the model
    coefMat <- matrix(rep(NA_real_, (1 + max(params$dgp$P) + max(params$dgp$Q)) * 2), nrow = 2)
    colnames(coefMat) <- c("(Intercept)", 
                           paste0("X", seq_len(max(params$dgp$P))), 
                           rep("score", max(params$dgp$Q)))
    
    # fill in regression coefficients
    # due to perfect collinearity of some predictor variables, regression coefficient 
    #   might become NA in some simulated conditions; therefore check is.na()
    nameIDX <- match(names(coef(fit)), colnames(coefMat))
    keep <- !is.na(nameIDX)
    coefMat[1, nameIDX[keep]] <- coef(fit)[keep]
    
    # fill in SE of the regression coefficient estimates
    vcovFit   <- tryCatch(vcov(fit), error = function(e) NULL)
    seFit   <- if (!is.null(vcovFit)) sqrt(diag(vcovFit)) else rep(NA_real_, length(coef(fit)))
    if (length(seFit) == length(coef(fit))) coefMat[2, nameIDX[keep]] <- seFit[keep]
    
    # object with results to return to predSetRes 
    list(
      # {RMSE, Rsquared, MAE} x {CVtest, train, HOtest}
      performVec = cbind(performCVtest, performTrain, performTest), 
      # prediction bias and prediction variance 
      predError = matrix(c(predBias, predVar), ncol = 2, nrow = 1),
      # predictor choice by extracting regression coefficients (for all possible coefficients)
      # matrix with regression coefficients in row 1 and SEs in row 2 
      regCoef = coefMat)
  })
  
  # build results object to return
  resOLS <- vector(mode = "list", 3)
  names(resOLS) <- c("performVec", "predError", "regCoef")
  
  # object with results to return to fitData Script 
  resOLS[["performVec"]] <- do.call(rbind, lapply(predSetRes, `[[`, "performVec"))
  # name columns of performance Matrix
  colnames(resOLS[["performVec"]]) <- c(unite(expand.grid(c("RMSE", "Rsquared", "MAE"), 
                                                        c("CVtest", "train", "HOtest")), 
                                            "names"))[["names"]]
  
  resOLS[["predError"]] <- do.call(rbind, lapply(predSetRes, `[[`, "predError"))
  resOLS[["regCoef"]] <- do.call(rbind, lapply(predSetRes, `[[`, "regCoef"))
  
  return(resOLS)
}