# Functions for Bias-Variance Script calculations

#Use this to calculate case-wise bias
calcCaseBias2 = function(T.scores,S.scores,reg.coef,T.coef,S.coef,W){#assumes W is correctly specified
  N = nrow(T.scores)
  caseBias=rep(NA,N)
  
  # applies W to beta_T to get a J x M beta_T vector
  beta_T.star = W%*%T.coef #used when reg.coef = items
  #applies W to beta_T to get a JxM beta_X vector -> used when reg.coef=composite
  beta_X.star = W%*%reg.coef
  
  for(i in 1:N){
    if(length(reg.coef)==length(S.coef)){#if item scores
      caseBias[i] = T.scores[i,]%*%t(W) %*%(beta_T.star-reg.coef) + 
        S.scores[i,]%*%(S.coef-reg.coef)
      
    }else if(length(reg.coef)==length(T.coef)){

      caseBias[i] = T.scores[i,]%*%t(W)%*%(beta_T.star-beta_X.star)+
        S.scores[i,]%*%(S.coef-beta_X.star)
    }
  }
  return(caseBias)
}

sampAvgBias2 = function(caseBias){
  return(sum(caseBias^2)/length(caseBias))
}

calcSampVar = function(p,N,sigmaEps){
  return((p/N)*sigmaEps)
}

avgPredErr = function(sigmaEps,avgbias.sq,var){
  return(sigmaEps+avgbias.sq+var)
}

popWeight = function(M,J){#populates a simple weight matrix for the dimensionality of T
  if(M==1){W=matrix(1,nrow=J,ncol=M)}
  else{
    W = matrix(0,nrow=J,ncol=M)
    if(J%%M==0){
      Jm = J/M
      start = 1
      for(m in 1:M){
        #print(start); print(start+Jm)
        W[start:(start+(Jm-1)),m] = rep(1,Jm)
        start = start+Jm
      }
    }else{
      Jm = floor(J/M)
      JM = J - Jm*(M-1)
      start = 1
      for(m in 1:(M-1)){
        #print(start); print(start+Jm)
        W[start:(start+(Jm-1)),m] = rep(1,Jm)
        start = start+Jm
      }
      W[start:J,M] = rep(1,JM)
    }
  }
  return(W)
}