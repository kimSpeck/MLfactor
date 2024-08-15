# attention: new psych-package version (2.4.6.26) estimates p-factor values differently
#           instead of correlations of ~ 0.5 for p-factor as criterion in the ENET
#           we only arrive at ~ 0.4 with the new version of the psych package
#   -> unklar, weshalb unterschiedliche Ergebnisse je nach psych-Version ausgegeben werden
# do not replicate Hang & Mottus (2023) analysis but improve the analysis
# 1. there was leakage in the data due to the imputation of missing values by the
#     median in the full data set instead of for training and test samples separately
# 2. by adding more predictors to the data, the performance in terms of R² decreases
#     which might be an indicator for suboptimal hyper paramater tuning in the ENET
#     -> there is no tuning for alpha
#     -> there is no custom vector of lambda values provided to the algorithm
# 3. the fixed alpha value for the ENET is alpha = 0.05 which is really close to 
#     the ridge regression (alpha = 0)
#     therefore, variable selection might be rather liberal
# 4. as hyperparameter selection criterion they use min instead of the one-standard-error
#     rule; therefore model overfitting might take place

# Are the items still superior to Facets and Domains if leakage and tuning are fixed?
#   -> yes. Items are still superior to Facets and Domains.
#   however, for some of the outcomes, using sex as a covariate tends to increase 
#     the test R² for domains > facettes > items
#   this only happend if the outcome was correlated with sex

# more testing:
# 0. standardize facets and domains; sum scores with less variance than items 
#   -> var(facets) & var(domains) < 1 < var(items); predictive superiority durch 
#      Varianzeinschränkung?
# 1. calculate factor scores with the Bartlett method. Use Bartlett factor scores
#     for prediction
# 2. use interaction with sex as additional predictors on item, facet and domain
#     level
# 3. investigate which predictors are selected by the ENET regression if we use
#     Items, Facets and domains as predictors of the outcomes.
#   -> does the algorithm prefer facets and domains over items due to increased
#      reliability?
#   -> why does algorithm prefer items as unreliable predictors over facets and 
#      domains?
# 4. use random forest as another algorithm! run random forest for items, facets 
#     and domains. Does superiority for items disappear? Why? non-linear associations?

# nItems = 240
# nFacets = 30
# nDomains = 5

setwd("../reanalysis/HangMottus2023")
source("analysisTools.R")

logFolder <- "log"
createFolder(logFolder)

plotFolder <- "plots"
createFolder(plotFolder)

#setup----
library(tidyverse)
library(lavaan)
library(psych)
library(careless)
library(EFAtools)
# model fitting
library(caret)
library(glmnet)
library(randomForest)
library(parallel)
# 
library(bootnet)
library(qgraph)
#library(sjPlot)
library(matrixStats)

adata_full = read.csv("Dataset/AAAData (1).csv")

################################################################################
#clean data----
################################################################################
## drop data with more than 40 missing personality item responses
adata_full = adata_full %>% filter( (adata_full %>% select(210:457) %>% is.na %>% rowSums) <= 40 )

## quality check 
adata_full = adata_full %>% select(-contains("SS_"))

# substance abuse scores ASSIST questionnaire (never used is coded as 0)
# substance abuse scores build the third factor in the bifactor model later that 
#     will be used as a criterion/target predicted by personality items/features 
adata_full$tobacco <- select(adata_full, c("assist_q2_tobacco", "assist_q3_tobacco", "assist_q4_tobacco", "assist_q5_tobacco", "assist_q6_tobacco", "assist_q7_tobacco")) %>% rowMeans
adata_full$alcohol <- select(adata_full, c("assist_q2_alcohol", "assist_q3_alcohol", "assist_q4_alcohol", "assist_q5_alcohol", "assist_q6_alcohol", "assist_q7_alcohol")) %>% rowMeans
adata_full$cannabis <- select(adata_full, c("assist_q2_cannabis", "assist_q3_cannabis", "assist_q4_cannabis", "assist_q5_cannabis", "assist_q6_cannabis", "assist_q7_cannabis")) %>% rowMeans

# recode missing values as never used = 0
adata_full[,451:453] = adata_full %>% select(451:453) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), 0))) 

## restrict age to 18 plus 
adata = adata_full %>% filter(agegroup == 2 & age > 17) %>% select(-contains("SS_"))

## reverse personality items, calculate alpha
p = select(adata, 210:449) 
item.names <- colnames(p)

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
  map(~psych::alpha(p %>% select(grep(.x, names(p)))))

map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

################################################################################
#EFA----
################################################################################
# factors from EFA will become the criterion values to predict with items, facets and domains
#subset DSM data
data <- adata[,c(27:46, 451:453)]
colnames(data)

#test suitablity for factor analysis
KMO(data)
BARTLETT(data)

#select number of factors
PARALLEL(data, cor_method = "spearman") # PCA-determined eigenvalues:  3
VSS(data)
eigenvalues<-eigen(cor(data))$values
plot(eigenvalues, type = 'b', pch = 16,
     main = "Scree Plot", xlab="", ylab="Eigenvalues")
axis(1, at = 1:10, labels = 1:10)

#build bifactor model with 3 specific factors
efa_dsm <- psych::omega(data, fm="minres", nfactors = 3, poly = F)

#extract factors scores for further analysis
efa_scores <- as.data.frame(psych::factor.scores(data,
                                                 efa_dsm$schmid$sl)$scores)

data$p_factor <- efa_scores$g # general factor from bifactor model
data$intern <- efa_scores$`F1*` 
data$psych <- efa_scores$`F2*` 
data$substance <- efa_scores$`F3*` # tobacco, alcohol, cannabis


alpha <- seq(0, 1, length.out = 20)
lambda <- 10^seq(-1, 1, length = 100)
# tuneGrid <- data.frame(cbind(alpha = alpha, lambda = lambda))
tuneGrid <- expand.grid(alpha = alpha, lambda = lambda)

items <- p
item.names <- colnames(items)
facet.names <- paste(rep(c("N","E","O","A","C"), each = 6), rep(1:6, times = 5), sep= "")
domain.names <- c("N","E","O","A","C")

sum(is.na(adata$sex)) # there is exactly one missing values
idxNA <- which(is.na(adata$sex)) # remove missing value from domains, facets and items (and sex + age)

items <- items[-idxNA,]
adata <- adata[-idxNA,]
data <- data[-idxNA,]

################################################################################
# fit full ENET model with items, facets and domains
# what will be selected? Does it depend on reliability?
################################################################################
# this procedures does not make any sense;
#   -> if items are chosen more frequently this happens due to specific variance 
#       in items that explains variance of the outcome beyond the domains or facets
#   -> thus, reliability might not be the critical feature to chose among predictors
#       rather despite lower reliability items are selected more frequently due to more 
#       variance
# Frage danach, ob Domänen mehr als ihre Teile/Items sind. Nein. Domänen sagen nichts
#     über die Einzelitems hinaus aus.
# -> das ist ja der Witz an einem Faktorenmodell. Der Faktor schlägt sich in Items 
#     nieder, aber sagt nichts über Items hinaus aus... 

fullENET <- function(x1, y, pTrain = 0.75) {
  # # test it
  # x1 <- items
  # y <- data$substance 
  # y <- data$intern 
  # pTrain = .75
  
  # split data in train & test
  dat <- splitData(x1, y, pTrain)

  # median imputation AFTER splitting
  dat <- medianImpute(dat$xtrain, dat$xtest, dat$ytrain, dat$ytest)

  # unweighted sum scores to calculate facets and domains
  facetsTrain <- getSumScores(dat$xtrain, facet.names)
  facetsTest <- getSumScores(dat$xtest, facet.names)
  
  domainsTrain <- getSumScores(dat$xtrain, domain.names)
  domainsTest <- getSumScores(dat$xtest, domain.names)
  
  # join items, facets and domains into one data set!
  trainData <- cbind(dat$xtrain, facetsTrain, domainsTrain)
  testData <- cbind(dat$xtest, facetsTest, domainsTest)
  
  # z-standardising data AFTER splitting
  # scale by default does centering and standardising
  # do this after imputation of missing values and separately for training and test data
  trainData = sapply(trainData, scale) # z-Standardisierung
  testData = sapply(testData, scale) # z-Standardisierung
  
  if (model == "glmnet") {
    # fit elastic net regression with standard tuning of lambda and fixed alpha (close to ridge)
    model.glmnet <- cv.glmnet(trainData, dat$ytrain, nfolds = 10, alpha = .05)
    
    lambdaMinIdx <- model.glmnet$index["min",]
    estBeta <- model.glmnet[["glmnet.fit"]][["beta"]][,lambdaMinIdx]
    selectedVars <- names(estBeta)[estBeta != 0]
    
    # fit <- glmnet(x = trainData,
    #               y = dat$ytrain, 
    #               alpha = 0.05,
    #               lambda = model.glmnet$lambda.min,
    #               family = "gaussian")
    # 
    # estBeta <- as.matrix(fit$beta)
    
    # predict criterion values from the model for the test data
    predTMP <- predict(model, testData, s = "lambda.min")
    predCor <- cor(predTMP, dat$ytest)
    
  } else if (model == "caret") {
    # # adjust variable names   
    # #   data.frame as input to caret needs identical variable names in train and test
    # #   data (which need to be fixed, otherwise covarTrain != covarTest)
    # covNames <- colnames(dat$xtrain)[!(colnames(dat$xtrain) %in% item.names)]
    # colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
    # colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
    
    # fit model (including cross-validation)
    fullData <- cbind(trainData, ytrain = dat$ytrain)
    model <- train(ytrain ~ .,
                   fullData,
                   method = "glmnet", #elastic net regularized regression
                   metric = "RMSE",
                   #preProc = "scale",
                   trControl = trainControl(method = "cv", 
                                            number = 10,
                                            selectionFunction = "best"), 
                   tuneGrid = tuneGrid,
                   verbose = FALSE) #model training using ten-fold cross-validation 
    
    estBeta <- coef(model$finalModel, model$bestTune$lambda)
    selectedVars <- rownames(estBeta)[which(estBeta != 0)]
    
    predTMP <- predict(model, testData, s = "lambda.min")
    predCor <- cor(predTMP, dat$ytest)
  }

  nItems <- sum(selectedVars %in% item.names)
  nFacets <- sum(selectedVars %in% facet.names)
  nDomains <- sum(selectedVars %in% domain.names)
  
  return(list(R2 = predCor**2, estBeta = estBeta, selectedVars = selectedVars,
              nItems = nItems, nFacets = nFacets, nDomains = nDomains))

}

################################################################################
# fit random forest for items, facets and domains (separately)
################################################################################
# are items still superior? Why? 
# run random forest regression for items, for facets and for domains separately

# to do: 
# fix the tuning grid

# default for regression is p/3
maxMtry <- (ncol(items)/3)
# maxMtry <- sqrt(ncol(items)) # for classification

tuneGridRF <- expand.grid(# number of variables randomly sampled as candidates at each split
                        mtry = c(1, 10, 30, 45, 60, 80), 
                        # number of trees to grow (how big is ensemble to average across) 
                        ntree = c(100, 300, 600, 1000, 2000), 
                        # minimum size of terminal nodes
                        nodesize = c(5, 25, 50, 100))

tundGridRF <- data.frame(mtry = c(1, 10, 30, 45, 60, 80))

randomForestPrediction <- function(x1, y, pTrain = 0.75) {
  # # test it
  # x1 <- items
  # y <- data$substance 
  # y <- data$intern 
  # pTrain = .75
  
  # split data in train & test
  dat <- splitData(x1, y, pTrain)
  
  # median imputation AFTER splitting
  dat <- medianImpute(dat$xtrain, dat$xtest, dat$ytrain, dat$ytest)
  
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
  
  # fit random forest with caret
  itemsR2 <- fitRF.caret(dat$xtrain, dat$ytrain, dat$xtest, dat$ytest, 
                           tuneGrid = tuneGridRF)
  facetsR2 <- fitRF.caret(facetsTrain, dat$ytrain, facetsTest, dat$ytest, 
                            tuneGrid = tuneGridRF)
  domainsR2 <- fitRF.caret(domainsTrain, dat$ytrain, domainsTest, dat$ytest, 
                             tuneGrid = tuneGridRF)
  
  return(cbind(domainsR2, facetsR2, itemsR2))
}

################################################################################
# fit model 
################################################################################
# -> missing data imputation for train and test separately
# -> z-Standardizing for train and test separately
# -> z-Standardizing for items, facets and domains

## choose between models:
# glmnet: only without leakage; thus, stick to Hang & Mottus cv! (i.e., fixed alpha and varying lambda)
# caret: without leakage + cross validation of lambda and alpha! 
# factorScores: use factor scores (Bartlett) instead of sum scores as predictors for facets and domains

## choose to include interactions between sex and every other predictor or not:
# ... by setting inter = FALSE or inter = TRUE

## for factor score calculation we need a factor model
# multiple group CFA cannot be estimated for test data due to small sample sizes 
#   -> 240 items and only 132 observations in one group 
factorNames <- list(first = facet.names, 
                    second = domain.names)

facetModel <- genModel(factorNames, p)

myENET <- function(x1, y, pTrain = .75, model = "caret", 
                   factorScores = FALSE, inter = FALSE){
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
  if (facetScores) {
    # fit CFA with second order factors based on covariance matrix 
    fitTrain <- lavaan::cfa(facetModel, 
                            data = dat$xtrain, 
                            se = "none")
    
    fitTest <- lavaan::cfa(facetModel, 
                           data = dat$xtest, 
                           se = "none")
    
    # predict factor scores based in cfa model of the data
    fsTrain <- lavPredict(fitTrain, method = "Bartlett") 
    fsTest <- lavPredict(fitTest, method = "Bartlett") 
    
    # extract only facet or only domain scores
    facetsTrain <- fsTrain[,colnames(fsTrain) %in% facet.names]
    facetsTest <- fsTest[,colnames(fsTest) %in% facet.names]
    domainsTrain <- fsTrain[,colnames(fsTrain) %in% domain.names]
    domainsTest <- fsTest[,colnames(fsTest) %in% domain.names]
    
  } else {
    # unweighted sum scores to calculate facets and domains
    facetsTrain <- getSumScores(dat$xtrain, facet.names)
    facetsTest <- getSumScores(dat$xtest, facet.names)
    
    domainsTrain <- getSumScores(dat$xtrain, domain.names)
    domainsTest <- getSumScores(dat$xtest, domain.names)
  }

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
  
  if (model == "glmnet") {
    # fit elastic net regression with standard tuning of lambda and fixed alpha (close to ridge)
    itemsR2 <- fitENET.glmnet(dat$xtrain, dat$ytrain, dat$xtest, dat$ytest, inter = inter)
    facetsR2 <- fitENET.glmnet(facetsTrain, dat$ytrain, facetsTest, dat$ytest, inter = inter)
    domainsR2 <- fitENET.glmnet(domainsTrain, dat$ytrain, domainsTest, dat$ytest, inter = inter)
  
  } else if (model == "caret") {
    # adjust variable names 
    #   data.frame as input to caret needs identical variable names in train and test
    #   data (which need to be fixed, otherwise covarTrain != covarTest)
    covNames <- colnames(dat$xtrain)[!(colnames(dat$xtrain) %in% item.names)]
    colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
    colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
    
    # fit glmnet with caret
    itemsR2 <- fitENET.caret(xtrain, ytrain, xtest, ytest, 
                             tuneGrid = tuneGrid, inter = inter)
    facetsR2 <- fitENET.caret(facetsTrain, ytrain, facetsTest, ytest, 
                              tuneGrid = tuneGrid, inter = inter)
    domainsR2 <- fitENET.caret(domainsTrain, ytrain, domainsTest, ytest, 
                               tuneGrid = tuneGrid, inter = inter)
  }


  # square correlations of predicted and true criterion values to arrive at R² values
  cbind(domainsR2, facetsR2, itemsR2)
}

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
# do analysis
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


# compare results with different predictor sets
pFactor <- rbindResults(resR2, "pFactor")
intern <- rbindResults(resR2, "intern")
psych <- rbindResults(resR2, "psych")
substance <- rbindResults(resR2, "substance")

predRes <- rbind(pFactor, intern, psych, substance)

idxCondVec <- c("onlySex", "none")

predRes$idxCond <- plyr::mapvalues(predRes$idxCond, 
                                   from = unique(predRes$idxCond), 
                                   to = idxCondVec)

# save(predRes, file = "Dataset/caretResultsPerSample.rda")
str(predRes)
colnames(predRes) <- c("domainsR2", "facetsR2", "itemsR2", "idxCond", "measure")

col2num <- c("domainsR2", "facetsR2", "itemsR2")
predRes[col2num] <- lapply(predRes[col2num], as.numeric)

col2fac <- c("idxCond", "measure")
predRes[col2fac] <- lapply(predRes[col2fac], factor)

predRes$diffDomains <- predRes$itemsR2 - predRes$domainsR2
predRes$diffFacets <- predRes$itemsR2 - predRes$facetsR2

range(predRes$diffDomains)
range(predRes$diffFacets)

plotRes <- aggregate(cbind(domainsR2, facetsR2, itemsR2, 
                           diffDomains, diffFacets) ~ idxCond + measure, 
                     data = predRes, function(x) {cbind(mean(x), 
                                                        sd(x),
                                                        quantile(x, 0.025),
                                                        quantile(x, 0.975))})

plotRes <- do.call(data.frame, plotRes)
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.1", "_M")
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.2", "_SE")
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.3", "_q025")
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.4", "_q975")

plotRes <- tidyr::pivot_longer(plotRes, 
                               cols = !c(idxCond, measure),
                               names_to = c("levels", "DV"), 
                               names_sep = "_",
                               values_to = "values")

plotRes <- tidyr::pivot_wider(plotRes,
                              names_from = DV, 
                              values_from = values)

str(plotRes)

# add original values
original <- cbind(domainsR2 = c(0.50, 0.03, 0.61, 0.03, 0.11, 0.04, 0.32, 0.03),
                  facetsR2 = c(0.51, 0.03, 0.67, 0.02, 0.22, 0.04, 0.41, 0.03), 
                  itemsR2 = c(0.53, 0.02, 0.69, 0.02, 0.27, 0.04, 0.48, 0.03), 
                  idxCond = rep("original", 4*2),
                  measure = rep(c("pFactor", "intern", "psych", "substance"), each = 2),
                  rowNames = rep(c("M", "SE"), times = length(criterionVec)))

original <- tidyr::pivot_longer(data.frame(original), 
                                cols = !c(idxCond, measure, rowNames),
                                names_to = "levels", 
                                values_to = "values")

original <- tidyr::pivot_wider(original,
                               names_from = rowNames, 
                               values_from = values)

plotRes$q025 <- NULL
plotRes$q975 <- NULL
plotRes <- rbind(plotRes[,], original)

str(plotRes)

col2num <- c("M", "SE")
plotRes[col2num] <- lapply(plotRes[col2num], as.numeric)

col2fac <- c("idxCond", "measure", "levels")
plotRes[col2fac] <- lapply(plotRes[col2fac], factor)

plotRes$idxCond <- factor(plotRes$idxCond, 
                          levels = c("original", "none", "onlyAge", "onlySex", "sexANDage"),
                          labels = c("original", "replicate", "+Age", "+Sex", "+Age+Sex"))

colValVec <- c("black", '#339900', '#0000FF', '#9900CC', '#FF0000')

# plot results
(pCompare <- ggplot(plotRes[plotRes$levels != "diffDomains" &
                              plotRes$levels != "diffFacets", ], 
                    aes(x = levels, y = M, 
                        group = idxCond, colour = idxCond, linetype = idxCond)) +
    geom_point(size = 2) +
    geom_line() +
    scale_color_manual(values = colValVec) +
    ylab("average test R²") +
    geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), width = 0.2) +
    facet_wrap(~ measure, scales = "free") + 
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

# # only leakage fixed (but still without alpha tuning)
# save(plotRes, file = "Dataset/plotNOleakResults.rda")
# ggplot2::ggsave(filename = paste0(plotFolder, "/noLeakRepResults.png"),
#                 plot = pCompare,
#                 width = 10.5,
#                 height = 9.41,
#                 units = "in")

# leakage fixed and tuning fixed
save(plotRes, file = "Dataset/plotCaretResults.rda")
ggplot2::ggsave(filename = paste0(plotFolder, "/caretRepResults.png"),
                plot = pCompare,
                width = 10.5,
                height = 9.41,
                units = "in")

# Take Away:
# leakage did not a lot of damage to the results
# adding age as another predictor does not change the R2
# adding sex improves the prediction for intern and substance
#   especially for substance, the improvement by adding sex as a predictor increases 
#   with increasing abstractness of the predictor set (improvement domains > facets > items)
#   -> this could be a sign of measurement invariance as the source of this increase 
#       in prediction 
# however, pFactor or psych are not affected by adding sex as another predictor
#   -> how do these criterion variables differ from each other? 
#   -> is sex correlated with substance or internalizing problems but not with the 
#       higher order pFactor or thought disorders?
# -> under which circumstances does is measurement invariance a plausible mechanism
#     to improve the prediction of single items?

# sex uncorrelated with criterion without any effect
# all correlations of age and any criterion are < 0.05
cor(sex, data$p_factor)
cor(sex, data$psych)
# sex correlated with criterion with effect on prediction
cor(sex, data$substance)
cor(sex, data$intern)

# man kann experimentell mit echten Daten nicht trennen, ob ...
#     a. ... Outcome und Covariate korreliert sind und Covariate in Items widergespiegelt wird
#             -> was man hier beobachtet, da Covariate zu mehr Anstieg in R2 führt in Domänen/Facetten 
#                 als in Items
#     b. ... Outcome und Covariate nicht korreliert sind aber Covariate in Items abgebildert wird bzw. MI verletzt
#             -> wenn Outcome und Covariate unkorreliert sind sollte kein Effekt bestehen außer durch
#                 Messinvarianzverletzung vermittelt
#             -> genau das kann man nur durch Simulation testen
#     c. ... Outcome und Covariate korreliert sind aber Covariate nicht in Items widergespiegelt wird
#             -> "normaler" Prädiktor

# scenario:
# andere Variablen, die mit Outcome korreliert sind und irgendwie in Items abgebildet 
#   werden als Erklärung für bessere prediction
# will man Informationen aus Items einfach nutzen?
#    pro: man muss sich keine Gedanken über weitere zu erhebende Variablen machen
#    contra: u.U. schwieriger zu erkennen, ob etwas durch stabile unveränderbare 
#             Quellen zustande kommt oder veränderbare Eigenschaften? Welcher Anteil
#             von Information im Item hat welchen Einfluss? Keine Trennung von 
#             unterschiedlichen Quellen?
# andere Variablen explizit erheben: 
#     pro: 

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
