#Bias - Variance Linear Regression Simulation
# open code from https://osf.io/s4fv9/
# replicate the following conditions that are printed in the paper 
#   N = {50, 100, 250, 500, 1000}
#   J = {9, 21}
#   betaS.cond = {Condition1: betaS = 0, Condition4: one betaS = 0.5/sigmaS, other betaS = 0}
#   betaT = 0.5 
#   variance components = {Condition1: sigmaT = 0.5, sigmaS = 0.25, sigmaE = 0.25}
#
# add measurement invariance condition:
#   two groups 
#   change the weight of one variable on the composite
#     one group with w9 = 1 and the other group with w9 = 0
#     --> wie lässt sich diese Manipulation mit Pratiwi et al. (2023) in Zusammenhang bringen?
# ! set the weight for sigmaS to zero to eliminate the item specific effect

# compare: (for variance components condition 1 as this is the plot in the paper)
#   no item effect (= item specific effects condition 1) vs. 
#   large effect for 1 item (= item specific effects condition 4) vs. 
#   measurement invariance violation of 1 item (= the MI equivalent of item specific effects condition 4)

# implement types of measurement invariance
# violating ...
#   ... configural MI by setting miWeight = 0 for specific items
#       streng genommen wissen wir über configural MI nichts, weil 1-faktorielles Modell!
#   ... metric MI by setting miWeight != {0, 1} for specific items
#   ... scalar MI by setting the mean value of the specific effect to another value? (check this)

# to do:
#   - in finaler Simualtion auf jeden Fall nicht nur 1 sample mit jeweiliger sample size;
#     simuliere auf jeden Fall ~ 1000 samples (auch Varianz in outcomes interessant?)
################################################################################
# implement MI by using row wise kronecker products to vary weighting of items 
#   depending on subgroups in the data

# test stuff
# # easy example of row wise kronecker product
# con <- matrix(c(1:6), byrow = 1, ncol = 1)
# weights <- matrix(c(1, 1, 
#                     1, 1, 
#                     1, 1, 
#                     1, 0, 
#                     1, 0, 
#                     1, 0), byrow = T, ncol = 2)
# grpnet::row.kronecker(con, weights)

################################################################################
require(MASS) #MASS for simulation
# require(caret) # enet
require(tidyr)
require(ggplot2)
require(grpnet)
library(parallel)
# source('../reanalysis/McClureAmmermanJacobucci2023/Simulation/Linear-Regression-Bias-Variance_functions.R')
setwd("reanalyses/McClureAmmermanJacobucci2024")
#setwd("Dokumente/Promotion/MLfactor/reanalysis/McClureAmmermanJacobucci2023/Simulation/")
source('original_Linear-Regression-Bias-Variance_functions.R')
source('simTools.R')
source('setParameters.R')


plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

# N = number of observations/training sample size 
# M = number of composites
# J = number of items/variables per composite

# to do: 
# reduce grid according to conditions
# test  betaS = 0 = betaS.cond == 0 vs. 
#       betaS = 0.5/sigmaS = betaS.cond == 4 vs. 
#       MI with betaS = 0 = betaS.cond == 0 vs.
#       MI with betaS = 0.5/sigmaS = betaS.cond == 4 (in real data MI violations and specific effects occur)
dim(params$gridFull) # 185000 conditions is too much

replicateFig2 <- params$gridFull[which(
  # only sigmaT = 0.5, sigmaS = 0.25, sigmaE = 0.25
  params$gridFull$measVar.cond == 1 &
    # true score effect = 0.5 (instead of {0, 0.25, 0.5, 0.75})
    params$gridFull$betaT == 0.5 &
    (params$gridFull$MI == TRUE & # no MI violations == MI is true
       # only all specific item effects are 0 (i.e., betaS.cond = 1) or 1 specific item effect is 0.5/sigmaS (betaS.cond = 4)
       (params$gridFull$betaS.cond == 1 | params$gridFull$betaS.cond == 4))),]

# 1. manipulate one item; different measurement invariance violations
#   -> this produces 01_plotMImanipulation
#   with condition 1 (no specific item effect) in the left column and condition 4
#     (large item specific effect on one item) in the second column from the left
#   other conditions are measurement invariance violations in one item instead of
#     the large effect on one item
#   for the scalar measurement invariance violation different means for the groups
#     means variance inflation! (specific item effect with M = 5 -> var(S) = 6.5
#     instead of 0.25)
simMIconditions <- params$gridFull[which(
  # only sigmaT = 0.5, sigmaS = 0.25, sigmaE = 0.25
  params$gridFull$measVar.cond == 1 &
    # true score effect = 0.5 (instead of {0, 0.25, 0.5, 0.75})
    params$gridFull$betaT == 0.5 &
    (params$gridFull$MI == FALSE & # MI violations == MI is violated (FALSE)
       # # mean difference between groups in scalar MIvio
       # # test different mean differences in scalar Mvio and therefore varying degree of variance
       # params$gridFull$mSpecific == 5 & 
       # # MIvio only in one item
       # # test different different numbers of MIvio items
       # params$gridFull$miItems.cond == 1 &
       # no specific item effects additionally to MIvio
       (params$gridFull$betaS.cond == 1 | params$gridFull$betaS.cond == 4))), ]

# every condition for multiple scalar MIvio condition (in terms of MIdifference, i.e., variance)
# every condition for different number of items with MI violation
# every condition for MIvio with and without large specific item effect
dim(simMIconditions) # 3600 MI violation conditions
dim(replicateFig2) # 100 replication conditions

params$subGrid <- rbind(replicateFig2, simMIconditions)
# enumerate conditions to keep track of the simulation progress
params$subGrid$cnt <- seq_along(params$subGrid$N.train) 

# 2. add conditions with different deviations in mean value of the specific item effect
#   -> this produces 02_plotMImeanSpecificEffect and 02_plotScalarMImean
#   increasing deviation of mean values in scalar measurement invariance violation
#     increase the superiority of items of composites!
# -> is this due to increase in variance?
#     with mean differences from 0 of {1, 2, 5, 10} the variance of the specific
#     item effect increases from 0.25 to {0.25, 1, 6.25, 25} for the first item!
#     -> therefore mean difference of 1 keeps variance in measurement invariant
#         item constant
# to do:
# -> what happens for the linear regression when all variables are standardized?
#     benefit for items disappears! benefit of items with measurement invariance
#       violations is only due to higher variance in variable
#     thus, for elastic net models, scalar MI is no potential explanation for superiority
#       of items as predictor variables are usually standardized
# -> what happend in model with different means but identical variance in the predictor variable?
# params$subGrid <- params$gridFull[which(params$gridFull$measVar.cond == 1 &
#                                           params$gridFull$betaT == 0.5 &
#                                           (params$gridFull$MI == FALSE &
#                                              params$gridFull$betaS.cond == 1 &
#                                              params$gridFull$miItems.cond == 1)), ]

# 3. add conditions with varying number of items that are subject to measurement invariance violations
# !!! with more than one variable with measurement invariance violations even for
#     the standardized version items are better than composites (with 2 or 3 MI violating variables)
# params$subGrid <- params$gridFull[which(params$gridFull$measVar.cond == 1 &
#                                           params$gridFull$betaT == 0.5 &
#                                           (params$gridFull$MI == FALSE &
#                                              params$gridFull$mSpecific == 1 &
#                                              # params$gridFull$mSpecific == 5 &
#                                              params$gridFull$betaS.cond == 1)), ]

# 4. add conditions with specific effect on outcome
# this simulated conditions provides a different view and might answer another question
# beforehand we tried to answer the question whether MI violations might be able to
#     show item superiority without any item specific effect present. 
#     we demonstrated that items become superior in terms of prediction over the composite if ... 
#     ... scalar measurement invariance violations increased the variance of the specific item effect
#     ... multiple items with scalar measurement invariance violations but with fixed variance of the specific item effect
# now we want to test if measurement invariance violations might cover up the superiority
#     of the items over composites as well (the aim is to understand and become able to explain or predict
#     which data might profit from item based or composite based prediction)
# params$subGrid <- params$gridFull[which(params$gridFull$measVar.cond == 1 &
#                                           params$gridFull$betaT == 0.5 &
#                                           ((params$gridFull$MI == TRUE &
#                                               (params$gridFull$betaS.cond == 1 |
#                                                  params$gridFull$betaS.cond == 4)) |
#                                              (params$gridFull$MI == FALSE &
#                                                 params$gridFull$miItems.cond == 1 &
#                                                 params$gridFull$mSpecific == 1 &
#                                                 (params$gridFull$betaS.cond == 1 |
#                                                    params$gridFull$betaS.cond == 4)))), ]

dim(params$subGrid)

##### test it #####
# N.train <- 1000
# J <- 9
# betaT <- 0.5
# seed <- params$subGrid$seed[1]
# measVar.cond <- 1
# sigmaT <- 0.5
# sigmaS <- 0.25
# sigmaS <- 0.25
# # 1. messinvariante Items ohne item specific effect
# miItems.cond <- 0
# betaS.cond <- 1
# MI <- T
# miWeight <- NA
# 2. messinvariante Items mit large effect of 1 Item
# miItems.cond <- 0
# betaS.cond <- 4 # betaS = {2.000000,  1.111111,  1.515152, 10.000000,1.000000}
# MI <- T
# miWeight <- 0
# # 3. mess-non-invariante Items
# betaS.cond <- 1
# MI <- F
# miItems.cond <- 1
# nGroups <- 2
# miWeight <- 0
# # 4. mess-non-invariante Items und large effect of 1 Item
# betaS.cond <- 4 # betaS = {2.000000,  1.111111,  1.515152, 10.000000,1.000000}
# MI <- F

# simFit <- function(N.train, J, betaT, betaS.cond, measVar.cond, MI, 
#                    nGroups, miItems.cond, miWeight, mSpecific, 
#                    sigmaErrY = 1, reflective = F, 
#                    # random seed by default; give seed for reproducibility purposes
#                    # think about this in parallel version with multiple repetitions
#                    seed = sample(3789941, size = 1)){

nCoresSampling <- 30
timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")
logFolder <- "logs"
createFolder(logFolder)

# fit multiple repetitions of every simulated condition in parallel
manageFitting <- function(N.train, J, betaT, betaS.cond, measVar.cond, MI, 
                          nGroups, miItems.cond, miWeight, mSpecific, 
                          # default values are overwritten by params$subGrid
                          sigmaErrY = 1, reflective = F, 
                          # random seed by default; give seed for reproducibility purposes
                          # seed is applied for each condition to generate RNG
                          conditionSeed = sample(3789941, size = 1), cnt){
  
  # copy entire manageFitting environment to the simFit function to be able to work
  #   with all parameters
  environment(simFit) <- environment()
  
  # Initiate cluster; type = "FORK" only on Linux/MacOS: contains all environment variables automatically
  cl <- makeCluster(nCoresSampling, type = "FORK",
                    outfile = paste0(logFolder, "/", "fitDataStatus",
                                     timeStampFolder, ".txt"))
  
  # set seed that works for parallel processing
  set.seed(conditionSeed)
  s <- .Random.seed
  clusterSetRNGStream(cl = cl, iseed = s)
  
  # condition start printen
  tmpTime <- format(Sys.time(), "%d%m%y_%H%M%S")
  print(paste0("running condition: ", cnt, " at ", tmpTime))
  
  res <- simFit() # run function to actually create dataset
  
  # close cluster to return resources (memory) back to OS
  stopCluster(cl)
  
  # write results of the condition into file
  fileName <- paste0("simResults_", timeStampFolder, ".csv")
  write.table(res, file = fileName, append = TRUE, 
              row.names = F, col.names = !file.exists(fileName), sep = ",")
  
  # condition stop printen
  tmpTime <- format(Sys.time(), "%d%m%y_%H%M%S")
  print(paste0("finishing condition: ", cnt, " at ", tmpTime))
  
  return(res)
}


simFit <- function(){
  
  sampleRes <- parSapply(cl, seq_len(params$nRep), function(iRep) {
    # set parameters
    N = params$holdout + N.train
    W = popWeight(params$M, J) #creates a weight matrix that partitions J in M~equal unit-weighted composites
    
    # note that beta isn't a standardized coefficient
    # beta_T = true score effect
    beta_T = rep(betaT, params$M)
    
    # set betaS to the correct value! The object to work with will be beta_S
    beta_S <- setBetaS(betaS.cond, params$sigma[measVar.cond,"sigmaS"], J, params)
    betaS <- params$betaS[betaS.cond]/params$sigma[measVar.cond,"sigmaS"]
    
    # number of items with measurement invariance violation depending on number of
    #   variables (J)
    # for J = 9: 0, 1, 3, 4
    # for J = 21: 0, 1, 7, 10
    miItems <- setMIitems(miItems.cond, J)
    
    # grouping variable to describe MI violation with 0 as "standard" and 1 "MI violation"
    ## for configural or metric invariance violations
    miGroups <- c(rep(c(0, 1), each = N.train/nGroups), 
                  rep(c(0, 1), each = params$holdout/nGroups))
    
    # weights for manifest variables on composite
    if (MI) { # measurement invariant items
      trueW <- W # unit weighting of all items (= basic sum score)
    } else {
      # first item with measurement invariance violation!
      # this is where configural and metric measurement invariance violations are implemented
      trueW <- matrix(c(
        # training sample
        rep(W, N.train/nGroups),               # weights as in MI condition
        rep(c(rep(miWeight, miItems),   # mi violation weights in mi violating items
              rep(1, J - miItems)), N.train/nGroups),
        # holdout sample
        rep(W, params$holdout/nGroups),        # weights as in MI condition
        rep(c(rep(miWeight, miItems),   # mi violation weights in mi violating items
              rep(1, J - miItems)), params$holdout/nGroups)), 
        byrow = T, ncol = J)
    }
    
    # set variances 
    if (params$M == 1) { # one component case = one value
      Sigma_T <- params$sigma[measVar.cond,"sigmaT"]
    } else { # multiple component case = matrix with values
      Sigma_T = matrix(params$sigma[measVar.cond,"sigmaT"], 
                       nrow=params$M , ncol=params$M)
      diag(Sigma_T) <- params$sigma[measVar.cond,"sigmaT"]
    }
    
    Sigma_S=diag(params$sigma[measVar.cond,"sigmaS"],nrow = J,ncol = J)
    Sigma_e=diag(params$sigma[measVar.cond,"sigmaE"],nrow = J,ncol = J)
    
    ##############################################################################
    # sample values: 
    # set seed now!
    # set.seed(seed)
    
    # simulate/sample true score, specific score and error values
    #     note that data sampling is already done for training and holdout data!
    #     check this! no difference between this approach and sampling twice from same dgm!
    #     Construct [N x M] = [N x 1]; dim(t(W)) = [1 x J] = [1 x 9]; 
    Construct = mvrnorm(N, rep(0,params$M), Sigma = Sigma_T) # sample true score data; NxM
    # ! true score is the same for all items
    
    # sample specific item effects (with or without violations of scalar measurement invariance)
    # to avoid higher variance in the specific item effect in conditions with scalar measurement 
    #     invariance violations as a confounding factor we have to model the deviations from 
    #     the mean as dummy regression
    # 
    if (MI == F & miWeight == 1) { # scalar measurement invariance violation
      # this is the scalar measurement invariance violation producing all the variance of 
      #   the specific item effect by group difference; not by variability between persons
      # # matrix(rep(rep(1:J), each = N.train/nGroups), ncol = J) # test
      SmiItems = rbind(matrix(rep(rep(params$mS, times = miItems), each = N.train/nGroups), ncol = miItems),
                       matrix(rep(rep(mSpecific, miItems),   # mi violation weights in mi violating items
                                  each = N.train/nGroups), ncol = miItems),
                       matrix(rep(rep(params$mS, times = miItems), each = params$holdout/nGroups), ncol = miItems),
                       matrix(rep(rep(mSpecific, miItems),   # mi violation weights in mi violating items
                                  each = params$holdout/nGroups), ncol = miItems))
      
      # to do: müsste das nicht umgekehrt sein? zuerst schwankende items, dann Gruppenitem?
      # S = cbind(mvrnorm(N, rep(params$mS, J - miItems), 
      #                   Sigma = Sigma_S[(miItems+1):J, (miItems+1):J]), 
      #           SmiItems)
      S = cbind(SmiItems, mvrnorm(N, rep(params$mS, J - miItems), 
                                  Sigma = Sigma_S[(miItems+1):J, (miItems+1):J]))
      

      
    } else {
      S = mvrnorm(N, rep(params$mS,J), Sigma = Sigma_S) # sample specific item data S; NxJ  
    }
    
    e = mvrnorm(N, rep(0,J), Sigma = Sigma_e) # sample error data e; NxJ
    
    # calculate observed item values (based on T, S and e)
    #   true score + specific effect + error
    #     Construct [N x M] = [N x 1]; dim(t(W)) = [1 x J] = [1 x 9]; 
    #     Construct %*% t(W) = [N x J]
    if (MI) { # measurement invariant items
      # manifeste Variablen als **weighted** True Composite Score + Item specific effect + error 
      X = Construct %*% t(trueW) + S + e # generates observed scores (NxJ) = manifest variables
      # standardizing at this point would mean leakage and S with more variance still provided to 
      #     calculation of Y!!!
    } else {
      # include measurment invariance by manipulating the weight matrix for generating X (= observed scores)
      # use row wise Kronecker Product to change the weights depending on group
      # need to know: 1:N.train = training data; (N.train+1):N = holdout sample
      # -> weights depending on group in training and in holdout sample
      # generate weight matrix for row wise kronecker product
      # [N x 1], [N x J]
      X = grpnet::row.kronecker(Construct, trueW) + S + e # generates observed scores (NxJ) = manifest variables
      # standardizing at this point would mean leakage and S with more variance still provided to 
      #     calculation of Y!!!
    }
    
    # calculate the composite score (i.e., sum score of manifest variables with W as unit vector)
    #     for the composite we do not use true W but unit weights for all items in every group
    # [N x J] [J x 1] = [N x 1]
    Xcomp = X %*% W # this is the composite (out of all manifest variables); matrix of composites scores weighted by W (NxM)
    
    # Generate Y [N x 1]
    err_y = rnorm(N, 0, sqrt(sigmaErrY))
    
    # to simulate y we use the true weights! 
    if (reflective) {
      # based on a reflective model instead of a formative model
      Y = Construct %*% beta_T + S %*% beta_S + err_y  
      R2.true = var(Construct %*% beta_T + S %*% beta_S) / var(Y) 
    } else {
      if (MI) {
        # [N x 1] [1 x J] [J x 1] [1 x 1] + [N x J] [J x 1] + [1 x 1]
        # [2000 x 1] [1 x 9] [9 x 1] [1 x 1] + [2000 x 9] [9 x 1] + [1 x 1]
        # t(trueW) %*% trueW = 9; therefore the true score enters 9 times in y
        #   -> this is rather a formative instead of a reflexice model
        # var(Construct %*% t(trueW) %*% trueW) = 41.20175
        # var(Construct %*% t(trueW) %*% trueW %*% beta_T) = 10.3
        # therefore, the specific item effect only contributes a really small amount
        #   of variance to the criterion/target
        # this is why the R² is so unbelievably high in every simulated condition
        Y = Construct %*% t(trueW) %*% trueW %*% beta_T + S %*% beta_S + err_y
        R2.true = var(Construct %*% t(trueW) %*% trueW %*% beta_T + S %*% beta_S) / var(Y)
        
      } else {
        # wenn hier W verwendet wird geht bei MI weight = 0 der Mittelwert des latenten Konstruktes (M(T) = 0) ein
        Y = grpnet::row.kronecker(Construct, trueW) %*% W %*% beta_T + S %*% beta_S + err_y
        R2.true = var(grpnet::row.kronecker(Construct, trueW) %*% W %*% beta_T + S %*% beta_S)/var(Y)
      }
    }
    
    ##############################################################################
    # Creating Train-Test sets
    # X = (manifest) item values = T + S + e 
    
    itms.df = data.frame(Y, X, miGroups)
    itms.train = itms.df[1:N.train, -dim(itms.df)[2]]
    itms.test = itms.df[(N.train+1):N,-dim(itms.df)[2]]
    itmsCov.train = itms.df[1:N.train, ]
    itmsCov.test = itms.df[(N.train+1):N, ]
    
    # ("manifest") composite values = (1*T1 + 1*T2 + ... + 1*Tj) + (S1 + S2 + ... + Sj) + 
    comp.df = data.frame(Y,Xcomp, miGroups)
    comp.train = comp.df[1:N.train,-dim(comp.df)[2]]
    comp.test = comp.df[(N.train+1):N,-dim(comp.df)[2]]
    compCov.train = comp.df[1:N.train,]
    compCov.test = comp.df[(N.train+1):N,]
    
    # true prediction by using Construct and specific item effect both as predictors
    true.df = data.frame(Y,Construct,S)
    true.train = true.df[1:N.train,]
    true.test = true.df[(N.train+1):N,]
    
    # Regressions
    lm.items = lm(Y~.,itms.train)
    lm.comp = lm(Y~.,comp.train)
    lm.itemsCov = lm(Y~.,itmsCov.train)
    lm.compCov = lm(Y~.,compCov.train)
    lm.true = lm(Y~.,true.train)
    
    #Calculate Holdout Prediction
    itms.pred = predict(lm.items,itms.test)
    comp.pred = predict(lm.comp,comp.test)
    itmsCov.pred = predict(lm.itemsCov,itmsCov.test)
    compCov.pred = predict(lm.compCov,compCov.test)
    true.pred = predict(lm.true,true.test)
    
    #Preparing Output
    resVec = as.matrix(rep(NA, 20), nrow = 1)
    dim(resVec) <- c(1, 20)
    
    #Condition Information
    resVec[1:14] <- c(iRep, N.train, J, betaT, betaS, measVar.cond, 
                      MI, nGroups, miItems.cond, miWeight, mSpecific, 
                      sigmaErrY, reflective, round(R2.true, 5))
    
    #Holdout RMSE
    resVec[15:19] <- round(c(calcRMSE(itms.pred,itms.test[,"Y"]), 
                       calcRMSE(comp.pred,comp.test[,"Y"]),
                       calcRMSE(itmsCov.pred,itmsCov.test[,"Y"]), 
                       calcRMSE(compCov.pred,compCov.test[,"Y"]),
                       calcRMSE(true.pred,true.test[,"Y"])), 5)
    
    # test R2
    resVec[20:24] <- round(c(evalPerformance(itms.pred,itms.test[,"Y"])["Rsquared"],
                       evalPerformance(comp.pred,comp.test[,"Y"])["Rsquared"],
                       evalPerformance(itmsCov.pred,itmsCov.test[,"Y"])["Rsquared"],
                       evalPerformance(compCov.pred,compCov.test[,"Y"])["Rsquared"],
                       evalPerformance(true.pred,true.test[,"Y"])["Rsquared"]), 5)
    
    resVec
  }, simplify = F)
  
  # res <- do.call(Map,c(c,sampleRes))
  res <- do.call(rbind, sampleRes)
  colnames(res) <-  c("Rep", "N","J","betaT","betaS","measVar",
                   "MI", "nGroups", "miItems", "miWeight", "mSpecific",
                   "sigmaEy", "formVSref", "simR2",
                   "item_rmse","comp_rmse", "itemCov_rmse","compCov_rmse",  "true_rmse",
                   "item_R2", "comp_R2", "itemCov_R2", "compCov_R2", "true_R2")
  # print(res)
  return(res)
}

# run simulation for all conditions in params$subGrid

# # to simulate from kind of formative measurement model (R^2 dependent on number of items)
# # with the error variance from Jacobucci on the outcome (sigma^2_e = 1): 
# #   ... R^2 is around 0.9 for J = 9
# #   ... R^2 is arount 0.99 for J = 21
# # due to incredibly high R² the specific item effect has almost no effect on true R²
# out <- do.call(mapply, c(FUN = simFun, params$subGrid))

# params$nRep <- 10
# replicateFig2$cnt <- seq_along(replicateFig2$N.train)
# out <- do.call(mapply, c(FUN = manageFitting, replicateFig2[1:10,], SIMPLIFY = F))
out <- do.call(mapply, c(FUN = manageFitting, params$subGrid, SIMPLIFY = F))

save(out, file = "simResults.rda")
################################################################################
# analyse results
################################################################################
load("simResults.rda")
out <- data.frame(t(out))
names(out) <-  c("Rep", "N","J","betaT","betaS","measVar",
                 "MI", "nGroups", "miItems", "miWeight", "mSpecific",
                 "sigmaEy", "formVSref", "simR2",
                 "item_rmse","comp_rmse", "itemCov_rmse","compCov_rmse",  "true_rmse",
                 "item_R2", "comp_R2", "itemCov_R2", "compCov_R2", "true_R2")

# the simulated R2 is really high! 
range(out$sim.R2)
# 0.8631196 1.0009733 # "formative"
# 0.0996968 0.5620854 # reflective with error in target
# 0.899537 1.002019 # reflective without error in target

rmCols <- c("Seed", "true.bias", "true.var", "true.err", "true.rmse", "true.R2",
            "sim.R2", "item.err","comp.err")
subset <- out
subset[, rmCols] <- list(NULL)

subset <- tidyr::pivot_longer(subset, 
                              !c(N, J, beta_S, beta_T, measVar, nGroups, miItems, miWeight, mSpecific),
                              names_to = "measures", values_to = "values")

subset <- separate(subset, measures, c("model", "measures"), sep = "\\.", remove = F)

str(subset)

subset$condition <- NA
subset$condition[which(subset$beta_S != 0)] <- "largeItemEffect"
subset$condition[which(subset$beta_S == 0 &
                         subset$nGroups == 1)] <- "noItemEffect"
subset$condition[which(subset$miWeight == 0)] <- "configuralMIvio"
subset$condition[which(subset$miWeight == -0.8)] <- "metricMIvio"
subset$condition[which(subset$miWeight == 1)] <- "scalarMIvio"

cols2fac <- c("N", "J", "model", "miItems", "mSpecific", "miWeight")
subset[cols2fac] <- lapply(subset[cols2fac], factor)
subset$measures <- factor(subset$measures, 
                          levels = c("R2", "rmse", "bias", "var"))
subset$condition <- factor(subset$condition, 
                           levels = c("noItemEffect", "largeItemEffect", 
                                      "configuralMIvio", "metricMIvio", "scalarMIvio"))

################################################################################
# plot results
################################################################################
# these are results graphics for different parameter manipulation grids to show
#   different effects on the superiority of items over composites

################################################################################
# # 1. manipulate one item; different measurement invariance violations
################################################################################
# plot RMSE, bias², var, R²
unique(subset$model)
subset$model <- factor(subset$model,
                       levels = c("comp", "comb", "item"))
(plotMImani <- ggplot(subset[subset$model != "comb",],
                      aes(y = values, x = N,
                          group = interaction(J, model),
                          shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    #scale_color_manual(values = c("#CC0000", "#FF3300", "#006600")) +
    scale_color_manual(values = c("#CC0000", "#006600")) +
    facet_grid(measures ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# # ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations.eps"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_reflective.eps"),
#                 plot = plotMImani,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# # ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_reflective.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_reflective_highR2.png"),
#                 plot = plotMImani,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# wie passt Abwesenheit von Unterschieden in prediction auf Basis von items und
#   composite mit Pratiwi Simulation zusammen?

# skalare Messinvarianzverletzung durch Mittelwertsunterschiede in specific item
#   effects zwischen Gruppen; aber specific item effects ohne Beziehung zu Outcome
# Wie wird durch Verletzung der skalaren MI item prediction > composite prediction?
# Wie ist das übertragbar/vereinbar mit der Vorstellung eines inhaltlichen Nuancen-Faktors?
#   wenn Mittelwertsunterschiede zwischen Items ohne Effekt auf Outcome prediction
#     vergrößern, brauchen wir dann items für Nuancen?
#   ist Nuance nichts anderes als die Variable, die die Mittelwertsunterschiede erlären kann?
#   Nuance als specific item effect als latente Variable?

################################################################################
# 2. add conditions with different deviations in mean value of the specific item effect
################################################################################
unique(subset$mSpecific)
unique(subset$model)
unique(subset$miItems)
subset$model <- factor(subset$model,
                       levels = c("comp", "comb", "item"))

(plotMIspecific <- ggplot(subset[which(subset$model != "comb" &
                                         subset$mSpecific != 10), ],
                          aes(y = values, x = N,
                              group = interaction(J, model, mSpecific),
                              shape = J, linetype = model, colour = mSpecific)) +
    geom_point() +
    geom_line() +
    #scale_linetype_manual(values = c("dashed", "dotted", "solid")) + 
    scale_color_manual(values = c("#0033FF", "#0033FF", "#006600", "#FF6600", "#CC0000")) +
    facet_grid(measures ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# # ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect.eps"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_reflective.eps"),
#                 plot = plotMIspecific,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# # ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_reflective.png"),
#                 plot = plotMIspecific,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# only condition with constant variance in specific item effect
(plotMIspecific_constant <- ggplot(subset[which(subset$model != "comb" &
                                                  (subset$mSpecific == 1 | 
                                                     subset$mSpecific == 0)), ],
                                   aes(y = values, x = N,
                                       group = interaction(J, model, mSpecific),
                                       shape = J, linetype = model, colour = mSpecific)) +
    geom_point() +
    geom_line() +
    #scale_linetype_manual(values = c("dashed", "dotted", "solid")) + 
    scale_color_manual(values = c("#0033FF", "#0033FF", "#006600", "#FF6600", "#CC0000")) +
    facet_grid(measures ~ condition, scales = "free_y", labeller = label_both, switch = "y"))


# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_fixedVarS.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_fixedVarS_reflective.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_fixedVarS_reflective_highR2.png"),
#                 plot = plotMIspecific_constant,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


# the additional manipulation of different means in scalar MI violation only affects the
#   scalar MI  violation condition -> only plot this subset of the data
# mSpecific = 1 means that the variance of the specific item effect is fixed to 0.25
onlyScalarMI <- subset[which(subset$condition == "scalarMIvio" &
                               subset$model != "comb"),]
(plotMIscalarSpecific <- ggplot(onlyScalarMI,
                                aes(y = values, x = N,
                                    group = interaction(J, model, mSpecific),
                                    shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    # scale_color_manual(values = c("#CC0000", "#FF3300", "#006600")) +
    scale_color_manual(values = c("#CC0000", "#006600")) +
    facet_grid(measures ~ mSpecific, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean.eps"),
#                 plot = plotMIscalarSpecific,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean_reflective.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean_reflective_highR2.png"),
#                 plot = plotMIscalarSpecific,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")



# with any kind of measurement invariance violation:
#   for small sample sizes variance in item based prediction is higher (for small sample sizes)
#   compared to composite based prediction; 
#   since item based prediction is worse for item based prediction supperior predictive 
#   performance of items cannot be effect of selective reporting
# the stronger mean deviations for specific item effects, the higher the bias in the composite
#   more bias for large number of items in the composite although number of items deviating from 
#   the mean of 0 is constant (only 1 item!) (this is counter intuitive!!!)
# the stronger mean deviations for specific item effects, the higher advantages in 
#   the predictive performance (R2)
#     - for specific effect 0 vs. 2 there is seems to be no advantage in item based 
#       prediction
#     - for specific effect 0 vs. {5, 10} there are advantages in predictive performance 
#       for items compared with the composite

################################################################################
# 3. add conditions with varying number of items that are subject to measurement invariance violations
################################################################################
unique(subset$miItems)
unique(subset$mSpecific)

# the actual number of items in miItems.cond = {2, 3} differs depending on the number of
#   items in the composite {J/3, floor(J/2)}
#   -> backtracking actual item numbers (depending on J) to item number conditions (independent of J)
subset$miItems.cond <- plyr::mapvalues(subset$miItems,
                                       from = c(0, 1, 3, 7, 4, 10),
                                       to = c(0, 1, 2, 2, 3, 3))

str(subset)
(plotMInrItems <- ggplot(subset[subset$model != "comb", ],
                         aes(y = values, x = N,
                             group = interaction(J, model, miItems.cond),
                             shape = J, linetype = model, colour = miItems.cond)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#0033FF", "#006600", "#FF6600", "#CC0000")) +
    facet_grid(measures ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems.eps"),
#                 plot = plotMInrItems,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_reflective.png"),
#                 plot = plotMInrItems,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


onlyScalarMI <- subset[which(subset$condition == "scalarMIvio"),]
(plotMIscalarNrItems <- ggplot(onlyScalarMI[onlyScalarMI$model != "comb",],
                               aes(y = values, x = N,
                                   group = interaction(J, model, miItems.cond),
                                   shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#006600")) +
    facet_grid(measures ~ miItems.cond, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotScalarMInrItems.eps"),
#                 plot = plotMIscalarNrItems,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotScalarMInrItems.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotScalarMInrItems_reflective.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotScalarMInrItems_reflective_highR2.png"),
#                 plot = plotMIscalarNrItems,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# the number of manipulated items does not seem to make any difference for configural 
#   and metric measurement invariance violations, i.e., the predictive performance of
#   the item based models do not surpass the composite based condition (unit weighted!)
# only for scalar measurement invariance violations the number of items with deviating
#   means of the specific effect on the items does make a difference in terms of 
#   predictive performance (R2, RMSE and bias)

################################################################################
# 4. add conditions with specific effect on outcome
################################################################################
unique(subset$miItems) # number of MI violating items is constant across all conditions
unique(subset$mSpecific) # mean value of scalar MI violating items only has one value
unique(subset$miWeight) # there are all types of MI violation (none, configural, metric, scalar) 
unique(subset$beta_S)
unique(subset$model)
subset <- subset[which(subset$model != "comb"),]

subset$beta_S <- factor(subset$beta_S)
subset <- tidyr::unite(subset, "model_betaS", c(model, beta_S), sep = " x ", remove = F)
unique(subset$model_betaS)
subset$model_betaS <- factor(subset$model_betaS, 
                             levels = c("item x 0", "item x 2", "comp x 0", "comp x 2"))
ggplot(subset,
       aes(y = values, x = N,
           group = interaction(J, model, model_betaS),
           shape = J, linetype = model, colour = model_betaS)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("#006600", "#33CC33", "#CC0000", "#FF3300")) +
  facet_grid(measures ~ condition, scales = "free_y", labeller = label_both, switch = "y")

(plotMIxBetaS <- ggplot(subset[which(subset$measures == "R2"),],
                        aes(y = values, x = N,
                            group = interaction(J, model),
                            shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#006600")) +
    facet_grid(beta_S ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS.eps"),
#                 plot = plotMIxBetaS,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS_reflective.png"),
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS_reflective_highR2.png"),
#                 plot = plotMIxBetaS,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# zusätzlich einen spezifischen Effekt dazu zu simulieren katapultiert bei scalar
#   MI violations die predictive performance unabhängig von Anzahl der Items weit über 
#   die predictive performance der composite
# bei configural und metric measurement invariance ist bei 21 items die predictive 
#   performance trotzdem nicht viel besser

# combined condition (items + factor) macht keinen Unterschied (unit weighting) 

################################################################################
# add regularized regressions?
# enet.items <-  train(Y ~ .,
#                      itms.train,
#                      method = "glmnet", #elastic net regularized regression
#                      metric = "RMSE",
#                      preProc = "scale",
#                      trControl = trainControl(method="cv", number=10), tuneLength = 21) #model training using ten-fold cross-validation
# 
# # enet.comp <-  train(Y ~ .,
# #                     comp.train,
# #                     method = "glmnet", #elastic net regularized regression
# #                     metric = "RMSE",
# #                     preProc = "scale",
# #                     trControl = trainControl(method="cv", number=10), tuneLength = 21) #model training using ten-fold cross-validation
# # 
# enet.true <-  train(Y ~ .,
#                     true.train,
#                     method = "glmnet", #elastic net regularized regression
#                     metric = "RMSE",
#                     preProc = "scale",
#                     trControl = trainControl(method="cv", number=10), tuneLength = 21) #model training using ten-fold cross-validation
# 
# pred.itemsTrain <- predict(enet.items, itms.train) 
# pred.itemsTest <- predict(enet.items, itms.test)
# 
# # pred.compTrain <- predict(enet.comp, comp.train) 
# # pred.compTest <- predict(enet.comp, comp.test)
# 
# pred.trueTrain <- predict(enet.true, true.train) 
# pred.trueTest <- predict(enet.true, true.test)