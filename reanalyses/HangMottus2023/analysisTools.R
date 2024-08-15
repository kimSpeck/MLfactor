################################################################################
# general setup functions
################################################################################
# function to create folder for output, plots, logfiles
createFolder <- function(folderPath) {
  if (!file.exists(folderPath)){
    dir.create(folderPath)
  }
}

################################################################################
# analysis helper functions
################################################################################

splitData <- function(x1, y, pTrain) {
  s = sample(nrow(x1), pTrain * nrow(x1)) 
  # if(is.null(x2)) x2 = x1 # copy x1 to get the test data if selection by index
  
  xtrain <- x1[s,] # item data with covariates
  ytrain <- y[s]
  
  xtest <- x1[-s,] # item data with covariates
  ytest <- y[-s]
  
  return(list(xtrain = xtrain, xtest = xtest,
              ytrain = ytrain, ytest = ytest))
}

medianImpute <- function(xtrain, xtest, ytrain, ytest) {
  # imputation with median = median replace (impute missing values with median values)! 
  # do this separately for training and test data!
  # all missing values are removed
  # missing values already removed from covariates 
  xtrain = data.frame(xtrain) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  xtest = data.frame(xtest) %>% 
    mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))
  
  return(list(xtrain = xtrain, xtest = xtest,
              ytrain = ytrain, ytest = ytest))
}

getSumScores <- function(data, var.names) {
  # calculate sum scores of var.names for data
  # generate empty data.frame to fill in 
  tmp <- as.data.frame(matrix(ncol = length(var.names), 
                              nrow = nrow(data))) %>% `names<-`(var.names)
  
  # calculate unweighted sum scores 
  for(i in seq_along(var.names)) {
    tmp[,var.names[i]] = data.frame(data) %>% 
      select(contains(var.names[i])) %>% 
      rowMeans(na.rm = TRUE)
  }
  
  return(tmp)
}

# # calculate facet and domain values for test and train separately
# facetsTrain <- as.data.frame(matrix(ncol = 30, nrow = nrow(xtrain))) %>% `names<-`(facet.names)
# facetsTest <- as.data.frame(matrix(ncol = 30, nrow = nrow(xtest))) %>% `names<-`(facet.names)
# domainsTrain <- as.data.frame(matrix(ncol = 5, nrow = nrow(xtrain))) %>% `names<-`(domain.names)
# domainsTest <- as.data.frame(matrix(ncol = 5, nrow = nrow(xtest))) %>% `names<-`(domain.names)
# 
# # unweighted sum scores to calculate facets and domains
# for(i in seq_along(facet.names)) {
#   facetsTrain[,facet.names[i]] = data.frame(xtrain) %>% select(contains(facet.names[i])) %>% rowMeans(na.rm = TRUE)
#   facetsTest[,facet.names[i]] = data.frame(xtest) %>% select(contains(facet.names[i])) %>% rowMeans(na.rm = TRUE)
# }
# 
# for (i in domain.names) {
#   domainsTrain[,i] = data.frame(xtrain) %>% select(contains(i)) %>% rowMeans(na.rm = TRUE)
#   domainsTest[,i] = data.frame(xtest) %>% select(contains(i)) %>% rowMeans(na.rm = TRUE)
# }

genInteractionData <- function(data, var.names) {
  # generate model formula for model including interactions with ONE covariate
  #   ! this code cannot generate interactions with multiple covariates
  # then, generate the data matrix including all specified interactions
  
  # find covariate in the data and save variable & covariate name
  covIdx <- colnames(data) %in% var.names
  covariate <- colnames(data)[!covIdx]
  
  # basic model (only predictors)
  model <- paste0("~ (0 +", paste(colnames(data)[covIdx],
                                  sep = "", collapse = "+"), ")") # all variables 
  
  # include interactions by multiplying every predictor in the model with the covariate
  model <- paste0(model, "*", covariate)
  # add covariate as linear predictor (main effect)
  model <- paste0(model, " + " , covariate)
  
  newData <- model.matrix(as.formula(model), data.frame(data))
  return(newData)
} 

fitENET.glmnet <- function(xtrain, ytrain, xtest, ytest, var.names, inter = FALSE){
  ###
  # train elastic net regression on training data, predict test data and calculate R2
  ## input:
  # xtrain      - [data frame] 
  # ytrain      - [vector]
  # xtest       - [data frame] predictor matrix (i.e., data frame) for train data
  # ytest       - [vector] criterion/target test data
  # tuneGrid    - [list]
  ## output:
  # predCor**2  - [scalar] R² as the squared correlation of the predicted and  
  #                         observed values of the test data set
  ###
  # is there a covariate in the data?
  existCovariate <- any(!(colnames(xtrain) %in% var.names))
  
  # if we want to include interactions and if there is a covariate in the data
  if (inter & existCovariate) { 
    # add interactions to the predictor matrix (train and test)
    xtrain <- genInteractionData(xtrain, var.names)
    xtest <- genInteractionData(xtest, var.names)
  }
  
  # fit elastic net regression with standard tuning of lambda and fixed alpha (close to ridge)
  model <- cv.glmnet(xtrain, ytrain, nfolds = 10, alpha = .05)
  
  # predict criterion values from the model for the test data
  predTMP <- predict(model, xtest, s = "lambda.min")
  predCor <- cor(predTMP, ytest)
  
  return(predCor**2)
}

# modelItems <- cv.glmnet(dat$xtrain, dat$ytrain, alpha = .05)
# modelFacets <- cv.glmnet(as.matrix(facetsTrain), dat$ytrain, alpha = .05)
# modelDomains <- cv.glmnet(as.matrix(domainsTrain), dat$ytrain, alpha = .05)
# 
# # predict criterion values from the model for the test data
# predItems <- predict(modelItems, dat$xtest, s = "lambda.min")
# itemsR2.b <- cor(predItems, dat$ytest)
# 
# predFacets <- predict(modelFacets, as.matrix(facetsTest), s = "lambda.min")
# facetsR2.b <- cor(predFacets, dat$ytest)
# 
# predDomains <- predict(modelDomains, as.matrix(domainsTest), s = "lambda.min")
# domainsR2.b <- cor(predDomains, dat$ytest)

fitENET.caret <- function(xtrain, ytrain, xtest, ytest, var.names, tuneGrid, inter = F){
  ###
  # train elastic net regression on training data, predict test data and calculate R2
  ## input:
  # trainData     - [data frame] 
  # xtest         - [data frame] predictor matrix (i.e., data frame) for train data
  # ytest         - [data frame] criterion/target test data
  # tuneGrid      - [list]
  ## output:
  # tmp**2        - [scalar] R² as the squared correlation of the predicted and  
  #                          observed values of the test data set
  ###
  # is there a covariate in the data?
  existCovariate <- any(!(colnames(xtrain) %in% var.names))
  
  # if we want to include interactions and if there is a covariate in the data
  if (inter & existCovariate) { 
    # add interactions to the predictor matrix (train and test)
    xtrain <- genInteractionData(xtrain, var.names)
    xtest <- genInteractionData(xtest, var.names)
  }
  
  # fit model (including cross-validation)
  trainData <- cbind(xtrain, ytrain = ytrain)
  model <- train(ytrain ~ .,
                 trainData,
                 method = "glmnet", #elastic net regularized regression
                 metric = "RMSE",
                 #preProc = "scale",
                 trControl = trainControl(method = "cv", 
                                          number = 10,
                                          selectionFunction = "best"), 
                 tuneGrid = tuneGrid,
                 verbose = FALSE) #model training using ten-fold cross-validation 
  
  predTMP <- predict(model, xtest, s = "lambda.min")
  predCor <- cor(predTMP, ytest)
  return(predCor**2)
}

fitRF.caret <- function(xtrain, ytrain, xtest, ytest, tuneGrid) {
  # fit model (including cross-validation)
  trainData <- cbind(xtrain, ytrain = ytrain)
  
  model.rf <- train(ytrain ~., 
                    trainData, 
                    method = 'rf', # random forest
                    # metric = 'RMSE', 
                    trControl = trainControl(method='cv', # repetition is done manually
                                             number=10, # nfolds
                                             selectionFunction = "best"),
                    # tuneGrid = tuneGrid,
                    tuneGrid = data.frame(mtry = 80,
                                          ntree = 500),
                    verbose = FALSE) 
  
  predTMP <- predict(model.rf, xtest)
  predCor <- cor(predTMP, ytest)
  return(predCor**2)
}

# for analysis with factor scores
genModel <- function(factorNames, items) {
  ###
  # generate higher order factor model in lavaan syntax
  ## input:
  # factorNames   - [data frame] 
  # items         - [data frame] predictor matrix (i.e., data frame) for train data
  ## output:
  # cfaModel      - [scalar] R² as the squared correlation of the predicted and  
  #                          observed values of the test data set
  ###
  varsPerFac <- sapply(factorNames[['first']], function(iFacet) {
    tmp = colnames(items %>% select(contains(iFacet)))
    paste0(sapply(tmp, function(iVar) {
      paste0(iVar, " +")}, simplify = "array"), collapse = " ")
  })
  cfaModel <- c()
  for (iFactor in factorNames[['first']]) {
    (cfaModel <- paste0(cfaModel, "\n", 
                        iFactor, " =~ ", 
                        substr(varsPerFac[[iFactor]], 1, nchar(varsPerFac[[iFactor]])-2)))
  }
  facPerDom <- sapply(factorNames[['second']], function(iDomain) {
    tmp = stringr::str_subset(factorNames[['first']], iDomain)
    paste0(sapply(tmp, function(iFac) {
      paste0(iFac, " +")}, simplify = "array"), collapse = " ")
  })
  for (iDomain in factorNames[['second']]) {
    (cfaModel <- paste0(cfaModel, "\n", 
                        iDomain, " =~ ", 
                        substr(facPerDom[[iDomain]], 1, nchar(facPerDom[[iDomain]])-2)))
  }
  return(cfaModel)
}

################################################################################
# analysis helper functions
################################################################################
rbindResults <- function(data, result) {
  ###
  # rbind data from nested list + add simulated condition info (N_pTrash)
  # choose one result matrix by name out of different results lists
  ## input:
  # data    - [list] nested list with N_pTrash as lists
  # result  - [character] name of specific outcome
  ## output:
  #         - [data.frame] result in one data frame including index for N_Trash Labels
  ###
  tmpRes <- do.call(rbind, lapply(seq_along(data), function(subList) {
    cbind(data[[subList]][[result]], 
          idxCond = rep(subList, nrow(data[[subList]][[result]])),
          measure = result)
  }))
  return(data.frame(tmpRes))
}

rowNames2col <- function(data, variableName = "rowNames"){
  ###
  # add rownames as another variable to the data frame/matrix 
  ## input:
  # data          - [data frame] data to add row names variable to
  # variableName  - [chr] variable name for new row names variable
  ## output:
  # data          - [data frame] data with row names variable
  ###
  # rownames to column in data frame
  data$rowNames <- rownames(data)
  # remove numbering (rownames need to be unique values, thus suffices .1, .2, ...)
  data$rowNames <- stringr::str_replace(data$rowNames, "\\.[:digit:]{1,}$", "")
  # name new column
  names(data)[names(data) == 'rowNames'] <- variableName
  return(data)
}

