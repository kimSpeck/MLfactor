# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

library(tidyverse)

plotFolder <- "plots"
createFolder(plotFolder)

plotVarCovFolder <- "varCovPlots"
# test it
# plotVarCovFolder <- "varCovPlotsExtra"
createFolder(paste0(plotFolder, "/", plotVarCovFolder))

dataFolder <- "extraData"
conditionFolderVec <- list.files(dataFolder, include.dirs = TRUE, recursive = FALSE)

covList <- vector(mode = "list", length(conditionFolderVec))
corList <- vector(mode = "list", length(conditionFolderVec))
names(covList) <- names(corList) <- conditionFolderVec

for (iCond in seq_along(conditionFolderVec)) {
  fileVec <- list.files(paste0(dataFolder, "/", conditionFolderVec[iCond]), 
                        pattern = "^sample_test", 
                        include.dirs = FALSE, recursive = FALSE)
  load(paste0(dataFolder, "/", conditionFolderVec[iCond], "/", fileVec))
  condVec <- parseFilename(paste0(conditionFolderVec[iCond], ".rds"))
  
  covData <- cov(data[,c(paste0("X", seq_len(as.numeric(condVec[["P"]]))), "Y", "group")])
  corData <- cor(data[,c(paste0("X", seq_len(as.numeric(condVec[["P"]]))), "Y", "group")])
  
  covList[[conditionFolderVec[iCond]]] <- covData
  corList[[conditionFolderVec[iCond]]] <- corData
  
  covDF <- reshape2::melt(covData)
  colnames(covDF) <- c("dim1", "dim2", "varCov")
  
  col <- colorRampPalette(c("#FFFFFF", "#FFCCCC", "#FF6666", 
                            "#FF0000", "#990000", "#4D0000"))
  # Plot similar to corrplot
  pCov <- ggplot(covDF, aes(x = dim1, y = dim2, fill = varCov)) +
    geom_tile() +
    scale_fill_gradientn(colours = col(200),
                         limits = c(min(covDF$varCov), max(covDF$varCov))) +
    geom_text(aes(label = round(varCov, 2)), size = 3) +
    scale_y_discrete(limits = rev(levels(covDF$dim2))) +  
    coord_fixed() +
    theme_minimal() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) 
  
  ggsave(filename = paste0(plotFolder, "/", plotVarCovFolder, 
                                    "/pVarCovObs_MIvio", condVec[["MIvio"]], 
                                    "_miItems", condVec[["miItems"]],  
                                    "_P", condVec[["P"]], "_cAlpha", condVec[["cAlpha"]], 
                                    "_varT", condVec[["varT"]], "_corGO", condVec[["corGO"]], 
                                    "_R2", condVec[["R2"]], ".png"),
                  plot = pCov,
                  width = 11.56,
                  height = 9.77,
                  units = "in")
}

# createFolder("results/dependentMeasures")
# save(corList, file = paste0("results/dependentMeasures/observedCorrelationPerCondition.rda"))
# save(covList, file = paste0("results/dependentMeasures/observedCovariancesPerCondition.rda"))

# test it
createFolder("extraResults/dependentMeasures")
save(corList, file = paste0("extraResults/dependentMeasures/observedCorrelationPerCondition.rda"))
save(covList, file = paste0("extraResults/dependentMeasures/observedCovariancesPerCondition.rda"))
