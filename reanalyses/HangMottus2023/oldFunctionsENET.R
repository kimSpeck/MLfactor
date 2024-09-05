################################################################################
# fit model only without leakage
# thus, stick to Hang & Mottus cv! (i.e., fixed alpha and varying lambda)
################################################################################
enetNOleak <- function(x1, y, pTrain = .75, inter = FALSE){
  # # test it
  # sex <- ifelse(adata$sex == 1, 0, 1) # as dummy
  # x1 <- cbind(items, sex)
  # y <- data$substance 
  # pTrain = .75
  
  # split data in train & test
  dat <- splitData(x1, y, pTrain)
  
  # median imputation AFTER splitting
  dat <- medianImpute(dat$xtrain, dat$xtest, dat$ytrain, dat$ytest)
  
  # calculate facet and domain values for test and train separately
  # unweighted sum scores to calculate facets and domains
  facetsTrain <- getSumScores(dat$xtrain, facet.names)
  facetsTest <- getSumScores(dat$xtest, facet.names)
  
  domainsTrain <- getSumScores(dat$xtrain, domain.names)
  domainsTest <- getSumScores(dat$xtest, domain.names)
  
  # z-standardising data AFTER splitting
  # scale by default does centering and standardising
  # do this after imputation of missing values and separately for training and test data
  dat$xtrain = sapply(dat$xtrain, scale) # z-Standardisierung
  dat$xtest = sapply(dat$xtest,scale) # z-Standardisierung
  
  facetsTrain = sapply(facetsTrain, scale) # z-Standardisierung
  facetsTest = sapply(facetsTest, scale) # z-Standardisierung
  
  domainsTrain = sapply(domainsTrain, scale) # z-Standardisierung
  domainsTest = sapply(domainsTest, scale) # z-Standardisierung
  
  # only covariates to add to facet and domain data (already standardized)
  covarTrain <- dat$xtrain[,!(colnames(dat$xtrain) %in% item.names)]
  covarTest <- dat$xtest[,!(colnames(dat$xtest) %in% item.names)] 
  
  # add covariates to facet and domain data
  # this works as well if there are no covariates in the data
  facetsTrain <- cbind(facetsTrain, covarTrain)
  domainsTrain <- cbind(domainsTrain, covarTrain)
  
  facetsTest <- cbind(facetsTest, covarTest)
  domainsTest <- cbind(domainsTest, covarTest)
  
  # fit elastic net regression with standard tuning of lambda and fixed alpha (close to ridge)
  itemsR2 <- fitENET.glmnet(dat$xtrain, dat$ytrain, dat$xtest, dat$ytest, inter)
  facetsR2 <- fitENET.glmnet(facetsTrain, dat$ytrain, facetsTest, dat$ytest, inter)
  domainsR2 <- fitENET.glmnet(domainsTrain, dat$ytrain, domainsTest, dat$ytest, inter)
  
  # square correlations of predicted and true criterion values to arrive at R² values
  cbind(domainsR2, facetsR2, itemsR2)
}

################################################################################
# fit model
# without leakage and with fixed cross-validation of alpha and lambda
# R² as outcome; not only correlation of predicted and actual criterion values
################################################################################

# fit data with caret (without leakage)
# -> missing data imputation for train and test separately
# -> z-Standardizing for train and test separately
enetCaret <- function(x1, y, pTrain = .75, inter = FALSE) {
  
  # split data in train & test
  dat <- splitData(x1, y, pTrain)
  
  # median imputation AFTER splitting
  dat <- medianImpute(dat$xtrain, dat$xtest, dat$ytrain, dat$ytest)
  
  # calculate facet and domain values for test and train separately
  # unweighted sum scores to calculate facets and domains
  facetsTrain <- getSumScores(dat$xtrain, facet.names)
  facetsTest <- getSumScores(dat$xtest, facet.names)
  
  domainsTrain <- getSumScores(dat$xtrain, domain.names)
  domainsTest <- getSumScores(dat$xtest, domain.names)
  
  # z-standardising data AFTER splitting
  # scale by default does centering and standardising
  # do this after imputation of missing values and separately for training and test data
  dat$xtrain = sapply(dat$xtrain, scale) # z-Standardisierung
  dat$xtest = sapply(dat$xtest,scale) # z-Standardisierung
  
  facetsTrain = sapply(facetsTrain, scale) # z-Standardisierung
  facetsTest = sapply(facetsTest, scale) # z-Standardisierung
  
  domainsTrain = sapply(domainsTrain, scale) # z-Standardisierung
  domainsTest = sapply(domainsTest, scale) # z-Standardisierung
  
  # only covariates to add to facet and domain data (already standardized)
  covarTrain <- dat$xtrain[,!(colnames(dat$xtrain) %in% item.names)]
  covarTest <- dat$xtest[,!(colnames(dat$xtest) %in% item.names)] 
  
  # add covariates to facet and domain data
  # this works as well if there are no covariates in the data
  facetsTrain <- cbind(facetsTrain, covarTrain)
  domainsTrain <- cbind(domainsTrain, covarTrain)
  
  facetsTest <- cbind(facetsTest, covarTest)
  domainsTest <- cbind(domainsTest, covarTest)
  
  # adjust variable names 
  #   data.frame as input to caret needs identical variable names in train and test
  #   data (which need to be fixed, otherwise covarTrain != covarTest)
  covNames <- colnames(xtrain)[!(colnames(xtrain) %in% item.names)]
  colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
  colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
  
  # fit glmnet with caret
  itemsR2 <- fitENET.caret(xtrain, ytrain, xtest, ytest, 
                           tuneGrid = tuneGrid, inter = inter)
  facetsR2 <- fitENET.caret(facetsTrain, ytrain, facetsTest, ytest, 
                            tuneGrid = tuneGrid, inter = inter)
  domainsR2 <- fitENET.caret(domainsTrain, ytrain, domainsTest, ytest, 
                             tuneGrid = tuneGrid, inter = inter)
  
  cbind(domainsR2, facetsR2, itemsR2)
}

# Warning message:
# In nominalTrainWorkflow(x = x, y = y, wts = weights, info = trainInfo,  :
#                         There were missing values in resampled performance measures.
#   This error happens when the model didn't converge in some cross-validation folds 
#   the predictions get zero variance. As a result, the metrics like RMSE or Rsquared 
#   can't be calculated so they become NAs
# -> these tuning parameters do not get chosen. but the model itself works fine
# taking a closer look on the cv models indicates that large values of lambda result 
#     in this issue; for small lambda values, every alpha converges

################################################################################
# run analysis with factor scores instead of sum scores 
# factor scores either for
################################################################################
# multiple group CFA cannot be estimated for test data due to small sample sizes 
#   -> 240 items and only 132 observations in one group 
factorNames <- list(first = facet.names, 
                    second = domain.names)

facetModel <- genModel(factorNames, p)

enetCaret.FS <- function(x1, y, p=.75, multipleGroup = FALSE) {
  # sample by default does sampling without replacement
  # split into train and validation data set by index
  s = sample(nrow(x1), p*nrow(x1)) 
  
  # if(is.null(x2)) x2 = x1 # copy x1 to get the test data if selection by index
  
  xtrain <- x1[s,]
  ytrain <- y[s]
  
  xtest <- x1[-s,]
  ytest <- y[-s]
  
  # imputation with median = median replace (impute missing values with median values)! 
  # do this separately for training and test data!
  # all missing values are removed
  xtrain = data.frame(xtrain) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  xtest = data.frame(xtest) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  # scale by default does centering and standardising
  # do this after imputation of missing values and separately for training and test data
  xtrain = sapply(xtrain, scale) # z-Standardisierung
  xtest = sapply(xtest,scale) # z-Standardisierung
  
  # only covariates to add to facet and domain data
  covarTrain <- xtrain[,!(colnames(xtrain) %in% item.names)]
  covarTest <- xtest[,!(colnames(xtest) %in% item.names)] 
  
  # calculate facet and domain values for test and train separately
  # fit CFA with second order factors based on covariance matrix 
  
  if (multipleGroup) {
    fitTrain <- lavaan::cfa(facetModel, 
                            data = xtrain, 
                            se = "none",
                            group = "sex")  
    fitTest <- lavaan::cfa(facetModel, 
                           data = xtest, 
                           se = "none",
                           group = "sex")
  } else {
    fitTrain <- lavaan::cfa(facetModel, 
                            data = xtrain, 
                            se = "none")
    fitTest <- lavaan::cfa(facetModel, 
                           data = xtest, 
                           se = "none")
  }
  
  # predict factor scores based in cfa model of the data
  fsTrain <- lavPredict(fitTrain, method = "Bartlett") 
  fsTest <- lavPredict(fitTest, method = "Bartlett") 
  
  # extract only facet or only domain scores
  facetsTrain <- fsTrain[,colnames(fsTrain) %in% facet.names]
  facetsTest <- fsTest[,colnames(fsTest) %in% facet.names]
  domainsTrain <- fsTrain[,colnames(fsTrain) %in% domain.names]
  domainsTest <- fsTest[,colnames(fsTest) %in% domain.names]
  
  # add covariates to facet and domain data
  # this works as well if there are no covariates in the data
  facetsTrain <- cbind(facetsTrain, covarTrain)
  domainsTrain <- cbind(domainsTrain, covarTrain)
  
  facetsTest <- cbind(facetsTest, covarTest)
  domainsTest <- cbind(domainsTest, covarTest)
  
  # adjust variable names 
  #   data.frame as input to caret needs identical variable names in train and test
  #   data (which need to be fixed, otherwise covarTrain != covarTest)
  covNames <- colnames(xtrain)[!(colnames(xtrain) %in% item.names)]
  colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
  colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
  
  # fit glmnet with caret
  itemsR2 <- fitENET(cbind(xtrain, ytrain = ytrain), 
                     xtest, ytest, tuneGrid = tuneGrid)
  facetsR2 <- fitENET(cbind(facetsTrain, ytrain = ytrain), 
                      facetsTest, ytest, tuneGrid = tuneGrid)
  domainsR2 <- fitENET(cbind(domainsTrain, ytrain = ytrain), 
                       domainsTest, ytest, tuneGrid = tuneGrid)
  
  cbind(domainsR2, facetsR2, itemsR2)
}

################################################################################
# run analysis with factor scores instead of sum scores 
# factor scores either for
################################################################################
# multiple group CFA cannot be estimated for test data due to small sample sizes 
#   -> 240 items and only 132 observations in one group 
factorNames <- list(first = facet.names, 
                    second = domain.names)

facetModel <- genModel(factorNames, p)

enetCaret.FS <- function(x1, y, p=.75, multipleGroup = FALSE) {
  # sample by default does sampling without replacement
  # split into train and validation data set by index
  s = sample(nrow(x1), p*nrow(x1)) 
  
  # if(is.null(x2)) x2 = x1 # copy x1 to get the test data if selection by index
  
  xtrain <- x1[s,]
  ytrain <- y[s]
  
  xtest <- x1[-s,]
  ytest <- y[-s]
  
  # imputation with median = median replace (impute missing values with median values)! 
  # do this separately for training and test data!
  # all missing values are removed
  xtrain = data.frame(xtrain) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  xtest = data.frame(xtest) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  # scale by default does centering and standardising
  # do this after imputation of missing values and separately for training and test data
  xtrain = sapply(xtrain, scale) # z-Standardisierung
  xtest = sapply(xtest,scale) # z-Standardisierung
  
  # only covariates to add to facet and domain data
  covarTrain <- xtrain[,!(colnames(xtrain) %in% item.names)]
  covarTest <- xtest[,!(colnames(xtest) %in% item.names)] 
  
  # calculate facet and domain values for test and train separately
  # fit CFA with second order factors based on covariance matrix 
  
  if (multipleGroup) {
    fitTrain <- lavaan::cfa(facetModel, 
                            data = xtrain, 
                            se = "none",
                            group = "sex")  
    fitTest <- lavaan::cfa(facetModel, 
                           data = xtest, 
                           se = "none",
                           group = "sex")
  } else {
    fitTrain <- lavaan::cfa(facetModel, 
                            data = xtrain, 
                            se = "none")
    fitTest <- lavaan::cfa(facetModel, 
                           data = xtest, 
                           se = "none")
  }
  
  # predict factor scores based in cfa model of the data
  fsTrain <- lavPredict(fitTrain, method = "Bartlett") 
  fsTest <- lavPredict(fitTest, method = "Bartlett") 
  
  # extract only facet or only domain scores
  facetsTrain <- fsTrain[,colnames(fsTrain) %in% facet.names]
  facetsTest <- fsTest[,colnames(fsTest) %in% facet.names]
  domainsTrain <- fsTrain[,colnames(fsTrain) %in% domain.names]
  domainsTest <- fsTest[,colnames(fsTest) %in% domain.names]
  
  # add covariates to facet and domain data
  # this works as well if there are no covariates in the data
  facetsTrain <- cbind(facetsTrain, covarTrain)
  domainsTrain <- cbind(domainsTrain, covarTrain)
  
  facetsTest <- cbind(facetsTest, covarTest)
  domainsTest <- cbind(domainsTest, covarTest)
  
  # adjust variable names 
  #   data.frame as input to caret needs identical variable names in train and test
  #   data (which need to be fixed, otherwise covarTrain != covarTest)
  covNames <- colnames(xtrain)[!(colnames(xtrain) %in% item.names)]
  colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
  colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
  
  # fit glmnet with caret
  itemsR2 <- fitENET(cbind(xtrain, ytrain = ytrain), 
                     xtest, ytest, tuneGrid = tuneGrid)
  facetsR2 <- fitENET(cbind(facetsTrain, ytrain = ytrain), 
                      facetsTest, ytest, tuneGrid = tuneGrid)
  domainsR2 <- fitENET(cbind(domainsTrain, ytrain = ytrain), 
                       domainsTest, ytest, tuneGrid = tuneGrid)
  
  cbind(domainsR2, facetsR2, itemsR2)
}

sex <- ifelse(adata$sex == 1, 0, 1) # as dummy
items <- cbind(items, sex)

# criterion to predict with ENET
criterionVec <- list(data$p_factor, data$intern, data$psych, data$substance)
nReps <- 100

nCoresSampling <- 10
timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
outfileCL <- paste0("log", "/", "fitDataStatus", timeStamp,".txt")
cl <- parallel::makeCluster(nCoresSampling, type = "FORK",
                            outfile = outfileCL)

# set seed that works for parallel processing
set.seed(3940904)
s <- .Random.seed
clusterSetRNGStream(cl = cl, iseed = s)


predResFS <- lapply(seq_along(criterionVec), function(iCrit) {
  
  tmp <- parLapply(cl, seq_len(nReps), function(iRep) {
    # tmp <- lapply(seq_len(nReps), function(iRep) {
    
    # run analysis with factor scores from CFA model (separately for train and test)
    #   with leakage fix and tuning fix
    enetCaret.FS(items, criterionVec[[iCrit]], p=.75)  # leakage and tuning fixed
  })
  
  topRes <- NCmisc::top(CPU = F,Table = T, procs = nProcs)$Table
  idxUSER <- which(topRes$"USER" == myUser & 
                     topRes$"COMMAND" == "rsession")
  sink(file = outfileCL, append = T, type = c("output", "message"))
  print(system('free -m', intern = T))
  # print(NCmisc::top(Table = T))
  print(topRes[idxUSER, c("PID", "USER", "%CPU", "%MEM")])
  sink(file = NULL)
  
  do.call(rbind, tmp)
})
names(predResFS) <- c("pFactor", "intern", "psych", "substance")

stopCluster(cl)

# results from factor scores are basically the same as for sum scores
# we cannot estimate multiple group CFA due to small sample sizes in one of the 
#     sex groups which means that there are more parameters than informations to 
#     fit the model
# we might try to fit the model for only one domain as multiple group and then 
#     predict the outcome based on this domain with multiple group factor scores
#     or factor scores that ignore measurement invariance violations
#   -> however, this does not have high priority 

################################################################################
# run analysis
################################################################################
# only add sex as a covariate; age did not improve prediction in any condition
#   age is not correlated with any of the outcomes... 


# add covariates that might result in measurement invariance
sex <- ifelse(adata$sex == 1, 0, 1) # as dummy
dataVec <- list(NULL, "sex") # thus reduce data set further in every iteration

# ! if there is no covariate we cannot add inter = TRUE since there are no interactions
# therefore remove these conditions (otherwise inter = FALSE condition would be run twice)
expand.grid(model = c("glmnet", "caret"), 
            factorScores = c(FALSE, TRUE),
            inter = c(FALSE, TRUE))

# criterion to predict with ENET
criterionVec <- list(data$p_factor, data$intern, data$psych, data$substance)
nReps <- 100

nCoresSampling <- 20
timeStamp <- format(Sys.time(), "%d%m%y_%H%M%S")
outfileCL <- paste0(logFolder, "/", "fitDataStatus", timeStamp, ".txt")
cl <- parallel::makeCluster(nCoresSampling, type = "FORK",
                            outfile = outfileCL)

nProcs = nCoresSampling
myUser = "kls"

# set seed that works for parallel processing
set.seed(3940904)
s <- .Random.seed
clusterSetRNGStream(cl = cl, iseed = s)

myUser <- "kls"
nProcs <- nCoresSampling

resR2 <- lapply(seq_along(dataVec), function(iData) {
  # add sex to the data
  items <- cbind(items, sex)
  
  topRes <- NCmisc::top(CPU = F,Table = T, procs = nProcs)$Table
  idxUSER <- which(topRes$"USER" == myUser & 
                     topRes$"COMMAND" == "rsession")
  sink(file = outfileCL, append = T, type = c("output", "message"))
  print(system('free -m', intern = T))
  # print(NCmisc::top(Table = T))
  print(topRes[idxUSER, c("PID", "USER", "%CPU", "%MEM")])
  sink(file = NULL)
  
  # remove extra predictors for each iteration (depending on condition)
  items[,c(dataVec[[iData]])] <- list(NULL)
  
  TMP <- lapply(seq_along(criterionVec), function(iCrit) {
    
    tmp <- parLapply(cl, seq_len(nReps), function(iRep) {
      #lapply(seq_len(nReps), function(iRep) {
      
      myENET(items, y = criterionVec[[iCrit]], pTrain = 0.75)
    })
    do.call(rbind, tmp)
  })
  names(TMP) <- c("pFactor", "intern", "psych", "substance")
  
  # remove all extra variables before the next iteration for different data set starts
  items[,c(dataVec[[2]])] <- list(NULL)
  
  TMP
})

names(resR2) <- c("sexANDage", "onlySex", "onlyAge", "none")

stopCluster(cl)