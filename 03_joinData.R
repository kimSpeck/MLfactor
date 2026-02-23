# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

# source("utils/simTools.R") 
# source("utils/fitTools.R")
library(tidyverse)

# resFolder <- "results"
resFolder <- "extraResults"
depFolder <- "dependentMeasures"
createFolder(paste0(resFolder, "/", depFolder))

fileVec <- list.files(resFolder, pattern = "^resFile_", 
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

################################################################################
# extract performance measures per sample
################################################################################

performMeasures <- matrix(NA, nrow = nRowPerform, 
                          ncol = (colsPerformVec + colsCondVec))

for (iFile in seq_along(fileVec)) {
  
  iRes <- readRDS(paste0(resFolder, "/", fileVec[iFile]))
  condVec <- matrix(unlist(parseFilename(fileVec[iFile])), ncol = colsCondVec, nrow = 1)
  
  # # check if file name information maps unambiguously to specific condition
  # # -> it does
  # params$dgp$grid[which(params$dgp$grid$MIvio == "specific" &
  #                         params$dgp$grid$miItems == 0 &
  #                         params$dgp$grid$P == 4 &
  #                         params$dgp$grid$Q == 1 &
  #                         params$dgp$grid$cAlpha == 0.6 &
  #                         params$dgp$grid$varT == 0.5 &
  #                         params$dgp$grid$corGO == 0 &
  #                         params$dgp$grid$R2 == 0.3), ]
  
  performMeasures[(((iFile-1)*rowsPerFile)+1):(iFile*rowsPerFile),
                  (colsCondVec + 1):(colsCondVec + colsPerformVec)] <- do.call(
    rbind, lapply(iRes, `[[`, "performVec"))
  performMeasures[(((iFile-1)*rowsPerFile)+1):(iFile*rowsPerFile),
                  seq_len(colsCondVec)] <- matrix(condVec, ncol = ncol(condVec), 
                                             nrow = rowsPerFile, byrow = T)
}

performMeasures <- data.frame(performMeasures)
# str(performMeasures)
names(performMeasures) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2", 
                            "sample", "N", "model", "predSet", 
                            "RMSE_CVtest", "Rsquared_CVtest", "MAE_CVtest", 
                            "RMSE_train", "Rsquared_train", "MAE_train",
                            "RMSE_HOtest", "Rsquared_HOtest", "MAE_HOtest")
col2num <- c("R2", "RMSE_CVtest", "Rsquared_CVtest", "MAE_CVtest", "RMSE_train", "Rsquared_train", 
             "MAE_train", "RMSE_HOtest", "Rsquared_HOtest", "MAE_HOtest")
performMeasures[col2num] <- sapply(performMeasures[col2num], as.numeric)

col2fac <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO", "sample", "N", "model", "predSet")
performMeasures[col2fac] <- sapply(performMeasures[col2fac], factor)

# save file with performance measures for every sample
save(performMeasures, file = paste0(resFolder, "/", depFolder, "/performPerSample.rda"))

################################################################################
# calculate measures of central tendency and dispersion 
################################################################################
# load(paste0(resFolder, "/", depFolder, "/performPerSample.rds"))
performMeasures <- tidyr::pivot_longer(performMeasures,
                                       cols = c("RMSE_CVtest", "Rsquared_CVtest", "MAE_CVtest", 
                                                "RMSE_train", "Rsquared_train", "MAE_train", 
                                                "RMSE_HOtest", "Rsquared_HOtest", "MAE_HOtest"),
                                       names_to = "measure_type",
                                       values_to = "values")

performMeasures$miItemsPlot <- ifelse(performMeasures$miItems %in% c(2, 4), "P/2", performMeasures$miItems)
performMeasures$N <- factor(performMeasures$N, levels = c(100, 300, 1000))

performMeasures <- tidyr::unite(data = performMeasures, col = "predSet_model", c(predSet, model), sep = "_")
performMeasures <- tidyr::pivot_wider(performMeasures, names_from = predSet_model,
                                      values_from = values)

# calculate difference between items and scores (LM) for all models and in every sample
#   items - score, thus, positive values indicates more explained variance based on items as predictors
performMeasures$diff_LM <- performMeasures$items_LM - performMeasures$score_LM
performMeasures$diff_Ridge <- performMeasures$items_Ridge - performMeasures$score_LM
performMeasures$diff_ENET <- performMeasures$items_ENET - performMeasures$score_LM

performMeasures <- tidyr::pivot_longer(performMeasures, 
                                       cols = !c("MIvio", "miItems", "P", "Q", "cAlpha", "varT",
                                                 "corGO", "R2", "sample", "N", "measure_type",
                                                 "miItemsPlot"),
                                       names_to = "model", values_to = "values")

# # save file with performance measures for every sample
# save(performMeasures, file = paste0(resFolder, "/", depFolder, "/performPerSample_differences.rda"))

performStats <- performMeasures %>% 
  group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, measure_type, model) %>% 
  summarize(M = mean(values, na.rm = T),
            SD = sd(values, na.rm = T)) 

# # save file with only data on central tendency and variance of the estimates
# save(performStats, file = paste0(resFolder, "/", depFolder, "/performStats.rda"))

rm(performStats, performMeasures)
gc()
################################################################################
# extract regression coefficient choice
################################################################################
colsCoefVec <- ncol(testSample[[1]][["regCoef"]])

coefMatrix <- matrix(NA, nrow = nRowPerform, ncol = (colsCoefVec + colsCondVec))

for (iFile in seq_along(fileVec)) {
  
  iRes <- readRDS(paste0(resFolder, "/", fileVec[iFile]))
  condVec <- matrix(unlist(parseFilename(fileVec[iFile])), ncol = colsCondVec, nrow = 1)
  
  coefMatrix[(((iFile-1)*rowsPerFile)+1):(iFile*rowsPerFile),
                  (colsCondVec + 1):(colsCondVec + colsCoefVec)] <- do.call(
                    rbind, lapply(iRes, `[[`, "regCoef"))
  coefMatrix[(((iFile-1)*rowsPerFile)+1):(iFile*rowsPerFile),
                  seq_len(colsCondVec)] <- matrix(condVec, ncol = ncol(condVec), 
                                                  nrow = rowsPerFile, byrow = T)
}

coefMatrix <- data.frame(coefMatrix)
# str(coefMatrix)
names(coefMatrix) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2", 
                       "sample", "N", "model", "predSet", "intercept", "X1", "X2", 
                       "X3", "X4", "X5", "X6", "X7", "X8", "score")

col2num <- c("R2", "intercept", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "score")
coefMatrix[col2num] <- sapply(coefMatrix[col2num], as.numeric)

col2fac <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO", "sample", "N", "model", "predSet")
coefMatrix[col2fac] <- sapply(coefMatrix[col2fac], factor)

# save file with performance measures for every sample
save(coefMatrix, file = paste0(resFolder, "/", depFolder, "/coefPerSample.rda"))

# LM: is score coefficient estimated in any model? perfect linear combination of the items!
#   -> no! score is always a perfect linear combination of the items! regression coefficient
#     therefore always NA
which(!is.na(coefMatrix[which(coefMatrix$model == "LM" &
                   coefMatrix$predSet == "comb"), "score"]))
all(is.na(coefMatrix[which(coefMatrix$model == "LM" &
                   coefMatrix$predSet == "comb"), "score"]))

rm(coefMatrix)
gc()
################################################################################
# extract hyperparameter choices
################################################################################
colsHyperVec <- ncol(testSample[[1]][["hyperParams"]])

# careful! this might change! (see 02_fitData.R)!
uniqueModels <- c("items", "comb", "items", "comb", "comb")

# fitted models x training sample size (N) x number of sample (nTrain) 
rowsPerFileHyper <- length(uniqueModels) * length(params$dgp$trainN) * params$sim$nTrain
nRowHyper <- rowsPerFileHyper * length(fileVec)
hyperSelection <- matrix(NA, nrow = nRowHyper, 
                          ncol = (colsHyperVec + colsCondVec))

for (iFile in seq_along(fileVec)) {
  
  iRes <- readRDS(paste0(resFolder, "/", fileVec[iFile]))
  condVec <- matrix(unlist(parseFilename(fileVec[iFile])), ncol = colsCondVec, nrow = 1)
  
  hyperSelection[(((iFile-1)*rowsPerFileHyper)+1):(iFile*rowsPerFileHyper),
             (colsCondVec + 1):(colsCondVec + colsHyperVec)] <- do.call(
               rbind, lapply(iRes, `[[`, "hyperParams"))
  hyperSelection[(((iFile-1)*rowsPerFileHyper)+1):(iFile*rowsPerFileHyper),
             seq_len(colsCondVec)] <- matrix(condVec, ncol = ncol(condVec), 
                                             nrow = rowsPerFileHyper, byrow = T)
}

hyperSelection <- data.frame(hyperSelection)
# str(hyperSelection)
names(hyperSelection) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2", 
                       "sample", "N", "model", "predSet", "alpha", "lambda")

col2num <- c("alpha", "lambda")
hyperSelection[col2num] <- sapply(hyperSelection[col2num], as.numeric)

col2fac <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO", "sample", "N", "model", "predSet", "R2")
hyperSelection[col2fac] <- sapply(hyperSelection[col2fac], factor)

# save file with performance measures for every sample
save(hyperSelection, file = paste0(resFolder, "/", depFolder, "/hyperPerSample.rda"))

