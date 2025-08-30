# fit linear regression with cross-validation to evaluate identical performance 
#     measures for different models and avoid overfitting

fitOLS <- function(trainData, testData, foldid) {
  
  # # check it
  # testData <- data[which(!(seq_len(dim(data)[1]) %in% sampleIdx)),]
  
  # checks
  if (length(foldid) != dim(trainData)[1]) {
    stop("Fold-ID vector and data have not the same length!")
  }
  
  # fit only items, only scores, fit items + scores
  # ... implement this by different vectors of variable names
  itemNames <- paste0("X", seq_len(P))
  trainPredSets <- list(itemNames, "score", c(itemNames, "score"))
  
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
    
    
    
    # object with results to return to predSetRes 
    list(
      # {RMSE, Rsquared, MAE} x {CVtest, train, HOtest}
      performVec = cbind(performCVtest, performTrain, performTest), 
      # prediction bias and prediction variance 
      predError = matrix(c(predBias, predVar), ncol = 2, nrow = 1),
      # predictor choice by extracting regression coefficients (for all possible coefficients)
      # here: due to perfect collinearity of some predictor variables, regression coefficient 
      #   might become NA in some simulated conditions
      regCoef = coef(fit) # regression coefficients
      summary(fit)$coefficients # SE?
      )
  })
  
  # object with results to return to fitData Script 
  # name columns of performance Matrix
  colNames <- c(unite(expand.grid(c("RMSE", "Rsquared", "MAE"), c("CVtest", "train", "HOtest")), "names"))[["names"]]
  
}