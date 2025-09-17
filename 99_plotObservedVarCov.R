# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

library(tidyverse)

dataFolder <- "data"
conditionFolder <- paste0("/simData_MIvio", MIvio, "_miItems", miItems,
                          "_P", P, "_Q", Q, "_cAlpha", cAlpha, "_varT", varT, 
                          "_corGO", corGO, "_R2", R2)

fileVec <- list.files(paste0(dataFolder, conditionFolder), pattern = "^sample_test", 
                      include.dirs = FALSE, recursive = FALSE)


# get result object dimensions
testSample <- readRDS(paste0(resFolder, "/", fileVec[1]))
colsPerformVec <- ncol(testSample[[1]][["performVec"]])
colsCondVec <- length(parseFilename(fileVec[1]))

# careful! this might change! (see 02_fitData.R)!
uniqueModels <- c("items", "score", "comb", "items", "comb", "items", "comb", "comb")

# fitted models x training sample size (N) x number of sample (nTrain) 
rowsPerFile <- length(uniqueModels) * length(params$dgp$trainN) * params$sim$nTrain
nRowPerform <- rowsPerFile * length(fileVec)