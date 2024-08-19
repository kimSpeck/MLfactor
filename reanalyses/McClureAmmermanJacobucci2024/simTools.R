# function to create folder for output, plots, logfiles
createFolder <- function(folderPath) {
  if (!file.exists(folderPath)){
    dir.create(folderPath)
  }
}

setBetaS <- function(betaS.cond, SigmaS, J, params) {
  # item specific effects; S - Y relationships
  if(betaS.cond == 1){ # all betaS = 0; no item effect 
    beta_S = rep(params$betaS[1], J)
  }else if(betaS.cond == 2){ # all betaS = 0.1/sigmaS; small item specific effect on every item
    val = params$betaS[2]/SigmaS #this needs modification for M > 1
    beta_S = rep(val, J)
  }else if(betaS.cond == 3){ # betaS = 0.3/sigmaS; moderate item specific effect on 1/3 of the items
    val = params$betaS[3]/SigmaS
    beta_S = c(rep(val, J/3), rep(0, (J-(J/3)))) #assumes j%%3==0
  }else if(betaS.cond == 4){ # betaS = 0.5/sigmaS; large item specific effect in one item
    val = params$betaS[4]/SigmaS
    beta_S = c(val, rep(0,J-1))
  }else if(betaS.cond == 5){
    val1 = params$betaS[2]/SigmaS
    val2= params$betaS[3]/SigmaS
    val3= params$betaS[4]/SigmaS
    beta_S = c(rep(val1, 3), rep(val2, 2), val3, rep(0,J-6))
  }
  return(beta_S)
}

setMIitems <- function(miItems.cond, J) {
  if (miItems.cond == 2) {
    miItems <- J/3
  } else if (miItems.cond == 3) {
    miItems <- floor(J/2)
  } else {
    miItems <- miItems.cond
  }
  return(miItems)
}

calcRMSE = function(pred,obs){
  return(sqrt(mean((pred-obs)^2)))
}


evalPerformance <- function(pred, obs) {
  ###
  # evaluate performance (R2, RMSE, MAE) for train or test data
  # see postResample function of the caret package
  # https://github.com/topepo/caret/blob/5f4bd2069bf486ae92240979f9d65b5c138ca8d4/pkg/caret/R/postResample.R#L126C3-L144C49
  # adjusted caret package function to return 0 instead of NA if there are no 
  #     predictors chosen in elastic net or lasso regression
  # to do: only do calculation if obs != factor && obs == numeric; otherwise we did 
  #         classification and need an pseudo R²?
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