# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R") # general functions for simulation management
source("utils/simTools.R") # general functions for data simulation

# functions to generate data
source("utils/sampleData.R")

# check if all packages are available; install packages that are not available
# default is to call functions from packages explicitly (with ::) instead of loading
#   packages into search space/attaching everything
packageVec <- c("parallel", "mvtnorm", "magic", "lavaan", "semTools", "tidyr")
sapply(packageVec, function(iPackage) {
  if (!require(iPackage, character.only = TRUE)) {
    install.packages(iPackage)
  }
})

library(magic)
library(lavaan)
library(tidyverse)

timeStampFolder <- format(Sys.time(), "%d%m%y_%H%M%S")

grid <- params$dgp$grid
str(grid) 

plotFolder <- "plots"
createFolder(plotFolder)
plotSigmaFolder <- "sigma"
createFolder(paste0(plotFolder, "/", plotSigmaFolder))

##### test it #####
# test simulation condition that crushes SEM estimation
# condVec <- which(params$dgp$grid$P == 4 &
#                    params$dgp$grid$cAlpha == 0.8 &
#                    params$dgp$grid$miItems == 1 &
#                    params$dgp$grid$corGO == 0 &
#                    params$dgp$grid$R2 == 0.8 &
#                    params$dgp$grid$varT == 0.25 &
#                    params$dgp$grid$MIvio == "metric")
# 
# iCond <- condVec[1]
# P <- params$dgp$grid$P[iCond]
# Q <- params$dgp$grid$Q[iCond]
# cAlpha <- params$dgp$grid$cAlpha[iCond]
# varT <- params$dgp$grid$varT[iCond]
# miItems <- params$dgp$grid$miItems[iCond]
# corGO <- params$dgp$grid$corGO[iCond]
# R2 <- params$dgp$grid$R2[iCond]
# varStable <- params$dgp$grid$varStable[iCond]
# flT <- params$dgp$grid$flT[iCond]
# flS <- params$dgp$grid$flS[iCond]
# varE <- params$dgp$grid$varE[iCond]
# MIvio <- params$dgp$grid$MIvio[iCond]
# b0 <- params$dgp$grid$b0[iCond]
# b1 <- params$dgp$grid$b1[iCond]
# betaT <- params$dgp$grid$betaT[iCond]
# varS <- params$dgp$grid$varS[iCond]
# betaS <- params$dgp$grid$betaS[iCond]
# sampleSeed <- grid$sampleSeed[iCond]

sigmaList <- vector(mode = "list", dim(grid)[1])

# ... always simulate nGroups = 2 in every condition but keep parameters identical! 
for (iCond in seq_len(dim(grid)[1])) {
  
  ##### parameter values #####
  P <- grid$P[iCond]
  Q <- grid$Q[iCond]
  cAlpha <- grid$cAlpha[iCond]
  varT <- grid$varT[iCond]
  miItems <- grid$miItems[iCond]
  corGO <- grid$corGO[iCond]
  R2 <- grid$R2[iCond]
  varStable <- grid$varStable[iCond]
  flT <- grid$flT[iCond]
  flS <- grid$flS[iCond]
  varE <- grid$varE[iCond]
  MIvio <- grid$MIvio[iCond]
  b0 <- grid$b0[iCond]
  b1 <- grid$b1[iCond]
  betaT <- grid$betaT[iCond]
  varS <- grid$varS[iCond]
  betaS <- grid$betaS[iCond]
  sampleSeed <- grid$sampleSeed[iCond]
  
  ##### generate matrices #####
  # generate lambda matrix
  lambdaSEM <- defineLambda(flT = flT, flS = flS, nGroups = params$dgp$nGroups, 
                            nItems = P, miItems = miItems, 
                            # for MIvio %in% c("specific", "scalar") miWeight is automatically flT
                            miWeight = flT + b1, nOutcome = params$dgp$nOutcome, 
                            nFactors = Q)
  
  nLatentSEM <- ncol(lambdaSEM)
  
  # generate covariance matrix for the factors
  #   variance of the outcome is unknwon but we only need residual variance of the
  #   outcome since y is an endogenuous variable in the structural model
  MIflag <- if (MIvio != "specific") T else F
  covSEM <- definePhi(nLatentSEM, varT = varT, varStable = varStable, varS = varS,
                      nItems = P, MIvio = MIflag, miItems = miItems, 
                      varNames = colnames(lambdaSEM))
  
  # generate residual variance matrix
  # I remove the specific item variance and instead model the MI violation with exactly
  #   the same variance as for the specific item variance
  #   thus, I do not need to adjust the residual variance of the MI violating item, right?
  thetaSEM <- defineTheta(varE = varE, nItems = P, 
                          nGroups = params$dgp$nGroups, nOutcome = params$dgp$nOutcome)
  
  ##### structure model #####
  # generate beta matrix
  betaSEM <- defineBeta(nLatentSEM, betaT = betaT, betaS = betaS, miItems = miItems, 
                        nItems = P, varNames = colnames(lambdaSEM), nFactors = 1)
  
  ##### complete SEM #####
  # calculate model implied covariance matrix for ...
  #     simulating from SEM by sampling from multivariate normal distribution 
  sigmaSEM <- getSigma(lambdaSEM, betaSEM, covSEM, thetaSEM)
  
  sigmaList[[iCond]] <- sigmaSEM
  
  # generate mean vector for data sampling
  # simulate group-outcome correlations by mean differences in the outcome
  #   -> calculate difference in outcome for target correlation assuming balanced design with 2 groups! 
  if (corGO != 0) {
    Y1g1Idx <- which(colnames(sigmaSEM) == "Y1_g1")
    Y1g2Idx <- which(colnames(sigmaSEM) == "Y1_g2")
    varWithinY <- 0.5 * sigmaSEM[Y1g1Idx, Y1g1Idx] + 0.5 * sigmaSEM[Y1g2Idx, Y1g2Idx]
    meanYg2 <- getMeanYg2(targetRho = corGO, params$dgp$varG, varWithinY)
  } else {
    # depending on correlation between group and outcome; if cor = 0: mYg2 = mOut
    meanYg2 <- params$dgp$mOut
  }
  
  M <- defineMeanVec(mItems = params$dgp$mItems, mOut = params$dgp$mOut, 
                     nItems = P, nGroups = params$dgp$nGroups, 
                     nOutcome = params$dgp$nOutcome, mXg2 = b0, mYg2 = meanYg2,
                     miItems = miItems)
  
  col <- colorRampPalette(c("#FFFFFF", "#FFCCCC", "#FF6666", "#FF0000", "#990000", 
                            "#4D0000"))
  # fileName <- paste0(plotFolder, "/", plotSigmaFolder, 
  #                   "/pSigma_MIvio", MIvio, "_miItems", miItems,  
  #                   "_P", P, "_cAlpha", cAlpha, "_varT", varT, 
  #                   "_corGO", corGO, "_R2", R2, ".png")
  
  sigmaDF <- reshape2::melt(sigmaSEM)
  colnames(sigmaDF) <- c("dim1", "dim2", "varCov")
  
  # Plot similar to corrplot
  pCorr <- ggplot(sigmaDF, aes(x = dim1, y = dim2, fill = varCov)) +
    geom_tile() +
    scale_fill_gradientn(colours = col(200),
                         limits = c(min(sigmaDF$varCov), max(sigmaDF$varCov))) +
    geom_text(aes(label = round(varCov, 2)), size = 3) +
    scale_y_discrete(limits = rev(levels(sigmaDF$dim2))) + 
    coord_fixed() +
    theme_minimal() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) 
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/", plotSigmaFolder, 
                                    "/pSigma_MIvio", MIvio, "_miItems", miItems,  
                                    "_P", P, "_cAlpha", cAlpha, "_varT", varT, 
                                    "_corGO", corGO, "_R2", R2, ".png"),
                  plot = pCorr,
                  width = 11.56,
                  height = 9.77,
                  units = "in")
  
  # png(file = fileName, width = 11.56, height = 9.77, type = "cairo")
  # corrplot::corrplot(sigmaSEM, method = "square", is.corr = FALSE, col = col(200),
  #                    addCoef.col = "black", # Add coefficient of correlation
  #                    tl.col="black", tl.srt=45, tl.cex = 1) #Text label color and rotation
  # dev.off()
  # 
  # pTMP <- {corrplot::corrplot(sigmaSEM, method = "square", is.corr = FALSE, col = col(200),
  #                    addCoef.col = "black", # Add coefficient of correlation
  #                    tl.col="black", tl.srt=45) #Text label color and rotation
  # recordPlot()}
  # ggsave(filename = fileName, plot = replayPlot(pTMP))
  
  # params$dgp$grid[which(params$dgp$grid$MIvio == "specific" &
  #                         params$dgp$grid$miItems == 0 &
  #                         params$dgp$grid$P == 4 &
  #                         params$dgp$grid$Q == 1 &
  #                         params$dgp$grid$cAlpha == 0.6 &
  #                         params$dgp$grid$varT == 0.5 &
  #                         params$dgp$grid$corGO == 0 &
  #                         params$dgp$grid$R2 == 0.3), ]
}

names(sigmaList) <- paste0("simData_MIvio", grid$MIvio, "_miItems", grid$miItems, 
                           "_P", grid$P, "_Q1_cAlpha", grid$cAlpha, "_varT", grid$varT, 
                           "_corGO", grid$corGO, "_R2", grid$R2)
createFolder("results/dependentMeasures")
save(sigmaList, file = paste0("results/dependentMeasures/modelImpliedCovariancesPerCondition.rda"))
