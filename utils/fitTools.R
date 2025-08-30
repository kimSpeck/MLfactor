evalPerformance <- function(pred, obs) {
  ###
  # evaluate performance (R2, RMSE, MAE) for train or test data
  # see postResample function of the caret package
  # https://github.com/topepo/caret/blob/5f4bd2069bf486ae92240979f9d65b5c138ca8d4/pkg/caret/R/postResample.R#L126C3-L144C49
  # adjusted caret package function to return 0 instead of NA if there are no 
  #     predictors chosen in elastic net or lasso regression
  ## input:
  # pred    - [vector] outcome as predicted based on model 
  # obs     - [vector] observed outcome 
  ## output:
  #         - [vector] c(RMSE, Rsquared, MAE)
  ###
  
  # # this is the original caret code; it produces NA if no predictor is chosen
  # #     with (silent) warning message standard deviation is zero
  # #   -> NAs bias the estimation of Rsquared in simulated conditions with low R2
  # resamplCor <- try(cor(pred, obs, use = "pairwise.complete.obs"), silent = TRUE)
  
  if (!is.factor(obs) && is.numeric(obs)) {
    # empirical R2
    # use tryCatch to set correlation and therefore R2 to zero if no predictor is 
    #     chosen in elastic net or lasso regression
    resamplCor <- tryCatch(cor(pred, obs, use = "pairwise.complete.obs"),
                           warning = function(w) 0, silent = TRUE) 
    mse <- mean((pred - obs)^2)
    mae <- mean(abs(pred - obs))
  } else {
    resamplCor <- NA
    mse <- NA
    mae <- NA
  }
  
  return(c(RMSE = sqrt(mse), Rsquared = resamplCor^2, MAE = mae))
}