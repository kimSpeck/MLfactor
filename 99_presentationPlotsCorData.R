# plot presentation plots for dppd slides
# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

library(tidyverse)

resFolder <- "results"
depFolder <- "dependentMeasures"

plotFolder <- "plots"
createFolder(plotFolder)
plotPresentFolder <- "present"
createFolder(paste0(plotFolder, "/", plotPresentFolder))

resFile <- paste0(resFolder, "/", depFolder, "/observedCorrelationPerCondition.rda")
load(resFile)

################################################################################
# plot correlation matrices
################################################################################

nameVec <- paste0("simData_MIvio", params$dgp$grid$MIvio, "_miItems", params$dgp$grid$miItems, 
       "_P", params$dgp$grid$P, "_Q1_cAlpha", params$dgp$grid$cAlpha, "_varT", 
       params$dgp$grid$varT, "_corGO", params$dgp$grid$corGO, "_R2", params$dgp$grid$R2)
subGrid <- params$dgp$grid[which(params$dgp$grid$P == 4 &
                                   params$dgp$grid$cAlpha == 0.8 &
                                   params$dgp$grid$varT == 0.25),]
subNameVec <- unique(paste0("simData_MIvio", subGrid$MIvio, "_miItems", subGrid$miItems, 
                     "_P4", "_Q1_cAlpha0.8", "_varT0.25", "_corGO", subGrid$corGO, 
                     "_R2", subGrid$R2))
subIDX <- which(names(corList) %in% subNameVec)
subCorData <- corList[subIDX]


nCondVec <- dim(parseFilename(paste0(names(subCorData)[[1]], ".rds")))[2]
nRowCor <- dim(reshape2::melt(subCorData[[1]]))[1]
nColCor <- dim(reshape2::melt(subCorData[[1]]))[2]
plotCorData <- matrix(NA, nrow = (nRowCor * length(subIDX)), ncol = (nCondVec + nColCor))
dim(plotCorData)

for (iList in seq_along(names(subCorData))) {
  plotCorData[(1+((iList-1)*nRowCor)):(iList*nRowCor), 
              c(seq_len(nCondVec))] <- matrix(rep(unlist(parseFilename(paste0(names(subCorData)[iList], ".rds"))), nRowCor), 
                                              byrow = T, ncol = nCondVec, nrow = nRowCor)
  plotCorData[(1+((iList-1)*nRowCor)):
                (iList*nRowCor), (nCondVec+1):(nCondVec+nColCor)] <- as.matrix(reshape2::melt(subCorData[[names(subCorData)[iList]]]))
}

colnames(plotCorData) <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2", "dim1", "dim2", "varCov")
plotCorData <- as.data.frame(plotCorData)
str(plotCorData)
plotCorData$varCov <- as.numeric(plotCorData$varCov)
plotCorData$MIvio <- factor(plotCorData$MIvio, 
                            levels = c("specific", "scalar", "metric"))
plotCorData$dim2 <- factor(plotCorData$dim2, 
                           levels = c("group", "Y", "X4", "X3", "X2", "X1"))
                           # levels = c("X1", "X2", "X3", "X4", "Y", "group"))
plotCorData$dim1 <- factor(plotCorData$dim1, 
                           #levels = c("group", "Y", "X4", "X3", "X2", "X1"))
                           levels = c("X1", "X2", "X3", "X4", "Y", "group"))

# col <- colorRampPalette(c("#FFFFFF", "#FFCCCC", "#FF6666", "#FF0000", "#990000", "#4D0000"))

# c("#FFFFFF", "#FFF0E0", "#FFD1A3", "#FFB366", "#FF944D", "#FF751A")
col <- colorRampPalette(c("#FFFFFF", "#E5F0FF", "#B3D1FF", "#80B3FF", "#4D94FF", "#1A75FF", "#0055CC"))

# overview: Plot similar to corrplot
for (iCor in unique(plotCorData$corGO)) {
  for (iR2 in unique(plotCorData$R2)) {
    pCorTmp <- ggplot(plotCorData[plotCorData$corGO == iCor &
                                    plotCorData$R2 == iR2, ], 
                      aes(x = dim1, y = dim2, fill = varCov)) +
      geom_tile() +
      scale_fill_gradientn(colours = col(200),
                           limits = c(min(plotCorData$varCov), max(plotCorData$varCov))) +
      geom_text(aes(label = round(varCov, 2)), size = 3) +
      # scale_y_discrete(limits = rev(levels(plotCorData$dim2))) +  
      coord_fixed() +
      theme_minimal() +
      xlab("") + ylab("") +
      facet_grid(miItems ~ MIvio) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            panel.grid = element_blank()) 
    
    ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                      "/pObservedCor_P4_cAlpha0.8_varT0.25_corGO", 
                                      iCor, "_R2", iR2, ".png"),
                    plot = pCorTmp,
                    width = 11.56,
                    height = 9.77,
                    units = "in")
  }
}

for (iMIvio in unique(plotCorData$MIvio)) {
  for (iCor in unique(plotCorData$corGO)) {
    if (iMIvio == "specific") {
      col <- colorRampPalette(c("#FFFFFF", "#E6F9E6", "#B3E6B3", "#80CC80", "#4DA64D", "#267326"))
    } else {
      col <- colorRampPalette(c("#FFFFFF", "#FFD6D6", "#FF9999", "#FF6666", "#CC3333", "#990000"))
      # col <- colorRampPalette(c("#FFFFFF", "#FFE5E5", "#FFB2B2", "#FF8080", "#FF4D4D", "#E60000"))
    }
    pCorSoloTMP <- ggplot(plotCorData[plotCorData$corGO == iCor &
                         plotCorData$R2 == 0.8 &
                         plotCorData$miItems == 2 &
                         plotCorData$MIvio == iMIvio, ], 
           aes(x = dim1, y = dim2, fill = varCov)) +
      geom_tile() +
      scale_x_discrete(position = "top") +
      scale_fill_gradientn(colours = col(200),
                           limits = c(min(plotCorData$varCov), max(plotCorData$varCov))) +
      geom_text(aes(label = round(varCov, 2)), size = 5) +
      coord_fixed() +
      theme_minimal() +
      xlab("") + ylab("") +
      guides(fill = "none") +
      theme(axis.text.x = element_text(hjust = 1),
            panel.grid = element_blank())
    
    (pCorSoloTMP <- themeFunctionPaper(pCorSoloTMP))
    ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                      "/pObservedCor_MIvio", iMIvio, 
                                      "_P4_cAlpha0.8_varT0.25_corGO", iCor, 
                                      "_R20.8.png"),
                    plot = pCorSoloTMP,
                    width = 11.56,
                    height = 9.77,
                    units = "in")
  }
}

