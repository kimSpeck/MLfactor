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

setwd("reanalyses/HangMottus2023")
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

#extract factors scores for furt"The tuning parameter grid should have columns"her analysis
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
#   -> thus, reliability might not be the critical feature to chose among predictors!
#       rather, despite lower reliability items are selected more frequently due to more 
#       variance
# Frage danach, ob Domänen mehr als ihre Teile/Items sind. Nein. Domänen sagen nichts
#     über die Einzelitems hinaus aus.
# -> das ist ja der Witz an einem Faktorenmodell. Der Faktor schlägt sich in Items 
#     nieder, aber sagt nichts über Items hinaus aus... 
# to test if more reliable predictors are chosen more frequently by the ENET, we 
#   need to do a simulation that does not confound reliability with explained variance 
#   in the outcome (check behaviour of alpha and lambda in this simulation, too)
#   -> take away: when comparing performance of different models (differing predictor sets)
#       which substantive conclusions are valid and which behaviour is only due to 
#       mathematical consequences in the behavior of the ENET?
#   -> N, simulate R², reliability of the predictors (value, mix: all equal, stepwise differences, ...), correlation of predictors 

fullENET <- function(x1, y, pTrain = 0.75, model = "caret") {
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
    predTMP <- predict(model.glmnet, testData, s = "lambda.min")
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
    model.caret <- train(ytrain ~ .,
                         fullData,
                         method = "glmnet", #elastic net regularized regression
                         metric = "RMSE",
                         #preProc = "scale",
                         trControl = trainControl(method = "cv", 
                                                  number = 10,
                                                  selectionFunction = "best"), 
                         tuneGrid = tuneGrid,
                         verbose = FALSE) #model training using ten-fold cross-validation 
    
    estBeta <- coef(model.caret$finalModel, model.caret$bestTune$lambda)
    selectedVars <- rownames(estBeta)[which(estBeta != 0)]
    
    predTMP <- predict(model.caret, testData, s = "lambda.min")
    predCor <- cor(predTMP, dat$ytest)
  }

  nItems <- sum(selectedVars %in% item.names)
  nFacets <- sum(selectedVars %in% facet.names)
  nDomains <- sum(selectedVars %in% domain.names)
  
  return(list(R2 = predCor**2, estBeta = estBeta, selectedVars = selectedVars,
              nItems = nItems, nFacets = nFacets, nDomains = nDomains))

}

# criterion to predict with ENET
modelVec <- c("glmnet", "caret") # which model is used, thus, cross validation of alpha or not
criterionVec <- list(data$p_factor, data$intern, data$psych, data$substance)
nReps <- 100

nCoresSampling <- 30
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

resBeta <- lapply(seq_along(modelVec), function(iModel) {

  TMP <- lapply(seq_along(criterionVec), function(iCrit) {
    
    tmp <- parLapply(cl, seq_len(nReps), function(iRep) {
      #lapply(seq_len(nReps), function(iRep) {
      fullENET(items, y = criterionVec[[iCrit]], 
               pTrain = 0.75, model = modelVec[iModel])
    })
  })
  names(TMP) <- c("pFactor", "intern", "psych", "substance")
  
  TMP
})

names(resBeta) <- modelVec

stopCluster(cl)

# save output to run something else but analyse later
save(resBeta, file = "Dataset/variableSelectionENET.rda")

# Auswertung später machen
# generell: es werden nicht alle Domänen gewählt; keine Präferenz für Domänen im
#           Vergleich zu Facetten oder Items
# wie oben erwähnt: nicht überraschend! 

################################################################################
# fit random forest for items, facets and domains (separately)
################################################################################
# are items still superior? Why? 
# run random forest regression for items, for facets and for domains separately

# to do: 
# fix the tuning grid

# # default for regression is p/3
# maxMtry <- (ncol(items)/3)
# # maxMtry <- sqrt(ncol(items)) # for classification
# 
# tuneGridRF <- expand.grid(# number of variables randomly sampled as candidates at each split
#                         mtry = c(1, 10, 30, 45, 60, 80), 
#                         # number of trees to grow (how big is ensemble to average across) 
#                         ntree = c(100, 300, 600, 1000, 2000), 
#                         # minimum size of terminal nodes
#                         nodesize = c(5, 25, 50, 100))
# 
# tundGridRF <- data.frame(mtry = c(1, 10, 30, 45, 60, 80))

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
                           tuneGrid = data.frame(mtry = c(10, 40, 80)))
  facetsR2 <- fitRF.caret(facetsTrain, dat$ytrain, facetsTest, dat$ytest, 
                            tuneGrid = data.frame(mtry = c(2, 5, 10)))
  domainsR2 <- fitRF.caret(domainsTrain, dat$ytrain, domainsTest, dat$ytest, 
                             tuneGrid = data.frame(mtry = c(1, 2)))
  
  return(cbind(domainsR2, facetsR2, itemsR2))
}

# add covariates that might result in measurement invariance
sex <- ifelse(adata$sex == 1, 0, 1) # as dummy
dataVec <- list(NULL, "sex") # thus reduce data set further in every iteration

criterionVec <- list(data$p_factor, data$intern, data$psych, data$substance)
nReps <- 100

nCoresSampling <- 30
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

resRF <- lapply(seq_along(dataVec), function(iData) {
  
  # add sex to the data
  items <- cbind(items, sex)

  # remove extra predictors for each iteration (depending on condition)
  items[,dataVec[iData]] <- list(NULL) 
  
  TMP <- lapply(seq_along(criterionVec), function(iCrit) {
    
    tmp <- parLapply(cl, seq_len(nReps), function(iRep) {
      randomForestPrediction(items, y = criterionVec[[iCrit]], pTrain = 0.75)
    })
    do.call(rbind, tmp)
  
  })  
  names(TMP) <- c("pFactor", "intern", "psych", "substance")
  
  # remove all extra variables before the next iteration for different data set starts
  items[,"sex"] <- list(NULL)
  
  TMP
})
names(resRF) <- c("sex", "none")

stopCluster(cl)

save(resRF, file = "Dataset/randomForestResults_withInter.rda")

# to do: 
# run random forests with sex as an additional predictor too; this is only without sex as an additional predictor
# compare results with different predictor sets

# # load data without covariate sex
# load("Dataset/randomForestResults.rda")


predResRF <- data.frame(rbind(cbind(resRF[["pFactor"]], rep("pFactor", times = 100)),
                              cbind(resRF[["intern"]], rep("intern", times = 100)),
                              cbind(resRF[["psych"]], rep("psych", times = 100)),
                              cbind(resRF[["substance"]], rep("substance", times = 100))))

# pFactor <- rbindResults(resRF, "pFactor")
# intern <- rbindResults(resRF, "intern")
# psych <- rbindResults(resRF, "psych")
# substance <- rbindResults(resRF, "substance")
# 
# predResRF <- rbind(pFactor, intern, psych, substance)

# idxCondVec <- c("onlySex", "none")
# idxCondVec <- condNames
# 
# predResRF$idxCond <- plyr::mapvalues(predResRF$idxCond, 
#                                      from = unique(predResRF$idxCond), 
#                                      to = idxCondVec)

# colnames(predResRF) <- c("domains", "facets", "items", "idxCond", "measure")

colnames(predResRF) <- c("domains", "facets", "items", "measure")

col2num <- c("domains", "facets", "items")
predResRF[col2num] <- lapply(predResRF[col2num], as.numeric)

# col2fac <- c("idxCond", "measure")
col2fac <- c("measure")
predResRF[col2fac] <- lapply(predResRF[col2fac], factor)

plotRes <- aggregate(cbind(domains, facets, items) ~ measure, 
                     data = predResRF, function(x) {cbind(mean(x), 
                                                          sd(x))})

plotRes <- do.call(data.frame, plotRes)
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.1", "_M")
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.2", "_SE")

plotRes <- tidyr::pivot_longer(plotRes, 
                               cols = !measure,
                               names_to = c("levels", "DV"), 
                               names_sep = "_",
                               values_to = "values")

plotRes <- tidyr::pivot_wider(plotRes,
                              names_from = DV, 
                              values_from = values)

(ggplot(plotRes, 
       aes(x = levels, y = M, 
           group = measure,
           colour = measure, linetype = measure)) +
  geom_point(size = 2) +
  geom_line() +
  #scale_color_manual(values = colValVec) +
  ylab("average test R²") +
  geom_errorbar(aes(ymin = M - 2*SE, ymax = M + 2*SE), width = 0.2) +
  theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
        panel.background = element_rect(color = "white", fill = "white"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15)))

# auch für Random Forest Regression verschwindet der Vorteil von Items nicht
#   Item Vorteil im Vergleich zu Composites nicht durch komplexe Interaktionen etc. 
#   erklärbar
# run this again with covariate (which already includes interactions in this case)

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
  if (factorScores) {
    # to do: 
    #   - automatically detect Heywood Case (lv variance is negative)!
    #     build try Catch and tell model what to do if Heywood case occurs
    #   - Bartlett methods return only 0s as factor scores for higher order factors. 
    #     why?
    #   - running code with factor scores throws this error:
    #     first error: x has missing values; consider using makeX() to impute them
    #     why?
    
    # fit CFA with second order factors based on covariance matrix 
    fitTrain <- lavaan::cfa(facetModel, 
                            data = dat$xtrain, 
                            se = "none")
    
    idxHeywood.train <- which(diag(fitTrain@Model@GLIST$psi) < 0)
    if (length(idxHeywood.train) > 0) {
      print(paste0("train: Heywood Case in variances of: ", 
                   lavNames(fitTrain, type = "lv")[idxHeywood.train]))
    }
      
    fitTest <- lavaan::cfa(facetModel, 
                           data = dat$xtest, 
                           se = "none")
    
    idxHeywood.test <- which(diag(fitTest@Model@GLIST$psi) < 0)
    if (length(idxHeywood.test) > 0) {
      print(paste0("train: Heywood Case in variances of: ", 
                 lavNames(fitTest, type = "lv")[idxHeywood.test]))
    }
    
    # predict factor scores based in cfa model of the data
    fsTrain <- lavPredict(fitTrain)
    fsTest <- lavPredict(fitTest)
    # # estimation of factor scores based on Bartlett method lead to only zeros in
    # #     the higher order factors (i.e., domains)
    # # why though? https://github.com/yrosseel/lavaan/blob/56b8a112d75793a7f8dff61ff4bb6eca412b9bab/R/lav_predict.R#L1059
    # fsTrain <- lavPredict(fitTrain, method = "Bartlett") 
    # fsTest <- lavPredict(fitTest, method = "Bartlett") 
    
    # extract only facet or only domain scores
    facetsTrain <- data.frame(fsTrain[,colnames(fsTrain) %in% facet.names])
    facetsTest <- data.frame(fsTest[,colnames(fsTest) %in% facet.names])
    domainsTrain <- data.frame(fsTrain[,colnames(fsTrain) %in% domain.names])
    domainsTest <- data.frame(fsTest[,colnames(fsTest) %in% domain.names])
    
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
    itemsR2 <- fitENET.glmnet(dat$xtrain, dat$ytrain, dat$xtest, dat$ytest, 
                              var.names = item.names, inter = inter)
    facetsR2 <- fitENET.glmnet(facetsTrain, dat$ytrain, facetsTest, dat$ytest, 
                               var.names = facet.names, inter = inter)
    domainsR2 <- fitENET.glmnet(domainsTrain, dat$ytrain, domainsTest, dat$ytest, 
                                var.names = domain.names, inter = inter)
  
  } else if (model == "caret") {
    # adjust variable names 
    #   data.frame as input to caret needs identical variable names in train and test
    #   data (which need to be fixed, otherwise covarTrain != covarTest)
    covNames <- colnames(dat$xtrain)[!(colnames(dat$xtrain) %in% item.names)]
    colnames(facetsTrain) <- colnames(facetsTest) <- c(facet.names, covNames) 
    colnames(domainsTrain) <- colnames(domainsTest) <- c(domain.names, covNames) 
    
    # fit glmnet with caret
    itemsR2 <- fitENET.caret(dat$xtrain, dat$ytrain, dat$xtest, dat$ytest, 
                             var.names = item.names, tuneGrid = tuneGrid, inter = inter)
    facetsR2 <- fitENET.caret(facetsTrain, dat$ytrain, facetsTest, dat$ytest, 
                              var.names = facet.names, tuneGrid = tuneGrid, inter = inter)
    domainsR2 <- fitENET.caret(domainsTrain, dat$ytrain, domainsTest, dat$ytest, 
                               var.names = domain.names, tuneGrid = tuneGrid, inter = inter)
  }


  # square correlations of predicted and true criterion values to arrive at R² values
  cbind(domainsR2, facetsR2, itemsR2)
}

################################################################################
# run analysis
################################################################################
# add covariates that might result in measurement invariance
sex <- ifelse(adata$sex == 1, 0, 1) # as dummy
dataVec <- list(NULL, "sex") # thus reduce data set further in every iteration

# ! if there is no covariate we cannot add inter = TRUE since there are no interactions
# therefore remove these conditions (otherwise inter = FALSE condition would be run twice)
conditionGrid <- expand.grid(model = c("glmnet", "caret"), 
                               factorScores = c(FALSE, TRUE),
                               inter = c(FALSE, TRUE),
                               covariate = dataVec) # ! NULL is with covariate; "sex" means that "sex" is removed
conditionGrid <- conditionGrid[!(conditionGrid$covariate == "sex" & 
                                   conditionGrid$inter == TRUE), ]
conditionGrid <- conditionGrid[conditionGrid$factorScores == FALSE, ]

# criterion to predict with ENET
criterionVec <- list(data$p_factor, data$intern, data$psych, data$substance)
nReps <- 100

nCoresSampling <- 30
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

resR2 <- lapply(seq_len(dim(conditionGrid)[1]), function(iCond) {
  print(paste0("now running: condition ", iCond, "; which means model ", conditionGrid[iCond, "model"],
        ", factorScores = ", conditionGrid[iCond, "factorScores"], 
        " and interaction = ", conditionGrid[iCond, "inter"]))
  
  # add sex to the data
  items <- cbind(items, sex)
  
  # remove extra predictors for each iteration (depending on condition)
  items[,conditionGrid[[iCond, "covariate"]]] <- list(NULL)
  
  TMP <- lapply(seq_along(criterionVec), function(iCrit) {
    
    tmp <- parLapply(cl, seq_len(nReps), function(iRep) {
      #lapply(seq_len(nReps), function(iRep) {
      
      myENET(items, y = criterionVec[[iCrit]], pTrain = 0.75,
             model = conditionGrid[iCond, "model"], 
             factorScores = conditionGrid[iCond, "factorScores"], 
             inter = conditionGrid[iCond, "inter"])
    })
    do.call(rbind, tmp)
  })
  names(TMP) <- c("pFactor", "intern", "psych", "substance")
  
  # remove all extra variables before the next iteration for different data set starts
  items[,"sex"] <- list(NULL)
  
  TMP
})

stopCluster(cl)

condNames <- paste0("model.", ifelse(conditionGrid[, "model"] == "glmnet", "noLeak", "alphaCV"), 
                    "_fs.", conditionGrid[, "factorScores"],
                    "_inter.", conditionGrid[, "inter"],
                    "_covariate.", ifelse(conditionGrid[, "covariate"] == "NULL", "sex", "none"))
names(resR2) <- condNames

# # save data 
# # 100 repetitions
# # noLeakage modelled with glmnet vs. model including cross validation of alpha in caret
# # with or without interaction of sex with every item, facet or domain, respectively
# save(resR2, file = "Dataset/resultsPerSample_sex_noLeakVSalphaCV_inter.rda")

################################################################################
# anaylse results
################################################################################
# load data and start from here 
load("Dataset/resultsPerSample_sex_noLeakVSalphaCV_inter.rda")

# compare results with different predictor sets
pFactor <- rbindResults(resR2, "pFactor")
intern <- rbindResults(resR2, "intern")
psych <- rbindResults(resR2, "psych")
substance <- rbindResults(resR2, "substance")

predRes <- rbind(pFactor, intern, psych, substance)

# idxCondVec <- c("onlySex", "none")
idxCondVec <- condNames

predRes$idxCond <- plyr::mapvalues(predRes$idxCond, 
                                   from = unique(predRes$idxCond), 
                                   to = idxCondVec)

# save(predRes, file = "Dataset/resultsPerSample_nReps5.rda")
str(predRes)
colnames(predRes) <- c("domains", "facets", "items", "idxCond", "measure")

col2num <- c("domains", "facets", "items")
predRes[col2num] <- lapply(predRes[col2num], as.numeric)

col2fac <- c("idxCond", "measure")
predRes[col2fac] <- lapply(predRes[col2fac], factor)

predRes$diffDomains <- predRes$items - predRes$domains
predRes$diffFacets <- predRes$items - predRes$facets

range(predRes$diffDomains)
range(predRes$diffFacets)

plotRes <- aggregate(cbind(domains, facets, items, 
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
original <- cbind(domains = c(0.50, 0.03, 0.61, 0.03, 0.11, 0.04, 0.32, 0.03),
                  facets = c(0.51, 0.03, 0.67, 0.02, 0.22, 0.04, 0.41, 0.03), 
                  items = c(0.53, 0.02, 0.69, 0.02, 0.27, 0.04, 0.48, 0.03), 
                  idxCond = rep("model.original_fs.FALSE_inter.FALSE_covariate.none", 4*2),
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

plotRes <- tidyr::separate(plotRes, idxCond, 
                           into = c("model", "factorScores", "interaction", "covariate"),
                           sep = "_",
                           remove = T)

plotRes$model <- stringr::str_remove(plotRes$model, "model.")
plotRes$factorScores <- stringr::str_remove(plotRes$factorScores, "fs.")
plotRes$interaction <- stringr::str_remove(plotRes$interaction, "inter.")
plotRes$covariate <- stringr::str_remove(plotRes$covariate, "covariate.")

plotRes <- tidyr::unite(plotRes, "modelConfig", c(model, covariate), 
                        sep = ".", remove = F)
unique(plotRes$modelConfig)
plotRes$modelConfig <- factor(plotRes$modelConfig, 
                          levels = c("original.none", "alphaCV.sex", "noLeak.sex", "alphaCV.none", "noLeak.none"))

colValVec <- c("black", '#339900', '#0000FF', '#9900CC', '#FF0000')

unique(plotRes$factorScores) # ignore factorScores
unique(plotRes$interaction) # ignore factorScores

# plot results
(pCompare <- ggplot(plotRes[plotRes$levels != "diffDomains" &
                              plotRes$levels != "diffFacets", ], 
                    aes(x = levels, y = M, 
                        group = interaction(modelConfig, interaction),
                        colour = modelConfig, linetype = interaction)) +
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

# # original data as correlation instead of R²
# save(plotRes, file = "Dataset/plotResults_sex_noLeakVSalphaCV_inter.rda")
# ggplot2::ggsave(filename = paste0(plotFolder, "/sex_noLeakVSalphaCV_inter.png"),
#                 plot = pCompare,
#                 width = 10.5,
#                 height = 9.41,
#                 units = "in")

# # with original data as squared correlations 
# ggplot2::ggsave(filename = paste0(plotFolder, "/sex_noLeakVSalphaCV_inter_originalR2.png"),
#                 plot = pCompare,
#                 width = 10.5,
#                 height = 9.41,
#                 units = "in")

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
cor(sex, data$p_factor) # 0.09802503
cor(sex, data$psych) # -0.02821687
# sex correlated with criterion with effect on prediction
cor(sex, data$substance) # -0.2542565
cor(sex, data$intern) # 0.2829462

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

