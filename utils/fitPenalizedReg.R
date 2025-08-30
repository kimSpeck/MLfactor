# fit penalized linear regression with cross-validation to evaluate identical performance 
#     measures for different models 
# set arguments in function to control regularization: 
#   ... ridge regression (alpha = 0)
#   ... elastic net regression 
#   ... elastic net regression with regularization restricted to items
#       thus, scores have to be included in the model and only deviations from 
#             score-based predictions are represented in items 

fitRegReg <- function(trainData, testData, foldid, regMode) {
  # checks
  if (length(foldid) != dim(trainData)[1]) {
    stop("Fold-ID vector and data have not the same length!")
  }
  
  if (!(regMode %in% c("ridge", "enet", "itemENET"))) {
    stop("Your specified regMode does not exist.")
  }
  
  # fit only items, only scores, fit items + scores
  # ... implement this by different vectors of variable names
  itemNames <- paste0("X", seq_len(P))
  # fit every combination of predictors for regMode = {ridge, enet}
  if (regMode != "itemENET") { 
    trainPredSets <- list(itemNames, "score", c(itemNames, "score"))  
  } else { # fit only item + scores for regMode = itemENET
    trainPredSets <- list(c(itemNames, "score"))
    
    # automatically set the penalty factor that controls regularization of only items
    # penalty factor = vector of length nPredictors in data (excluding intercept; intercept never penalized)
    # default penalization: 1; unpenalized features: 0
    #   Unpenalized features are still estimated jointly, so their coefficients can 
    #   change with lambda because the penalized part of the model changes; they’re just not directly shrunk.
    pf <- vector()
  }
  
  
  # option to use penalty only on items
  # - penalty.factor is a numeric vector of length p (number of columns in X, excluding the intercept).
  # - Default is 1 for all features.
  # - Set penalty.factor[j] = 0 to leave feature j unpenalized (it’s estimated but not shrunk).
  # - The intercept is never penalized (when intercept=TRUE).
  
}