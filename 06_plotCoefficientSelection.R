# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

# load packages
library(tidyverse)
library(see) # violin half plots 

resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/coefPerSample.rda")
load(resFile)

plotFolder <- "plots"
createFolder(plotFolder)
plotCoefFolder <- "coefficients"
createFolder(paste0(plotFolder, "/", plotCoefFolder))

################################################################################
# plot (all) regression coefficients as violin half plots
################################################################################
# plot only N = 1000
# add LM regression coefficients to the plots 
# plot all combinations of cAlpha, varT, R2 and corGO
coefMatrix$miItemsPlot <- ifelse(coefMatrix$miItems %in% c(2, 4), "P/2", coefMatrix$miItems)
coefMatrix$N <- factor(coefMatrix$N, levels = c(100, 300, 1000))
coefMatrix$MIvio <- factor(coefMatrix$MIvio, levels = c("specific", "metric", "scalar"))
coefMatrix$model <- factor(coefMatrix$model, levels = c("itemENET", "ENET", "Ridge", "LM"))

coef_pDiffSub <- coefMatrix[coefMatrix$cAlpha == 0.8 & # for ENET
                        #coefMatrix$model %in% c("ENET", "itemENET", "Ridge") & # for ENET
                        coefMatrix$varT == 0.25, ] # due to scaling numerically bigger

coef_pDiffSub <- tidyr::pivot_longer(coef_pDiffSub, 
                               cols = !c(MIvio, miItemsPlot, miItems, P, Q, cAlpha, varT, corGO,
                                         R2, sample, N, model, predSet),
                               names_to = "coefficient", values_to = "value")
coef_pDiffSub$coefficient <- factor(coef_pDiffSub$coefficient, 
                                    levels = c("X1", "X2", "X3", "X4", "X5", "X6", 
                                               "X7", "X8", "score", "intercept"))

pVec <- unique(coef_pDiffSub$P)
miItemsVec <- unique(coef_pDiffSub$miItemsPlot)
corGOvec <- unique(coef_pDiffSub$corGO)
# r2vec <- unique(coef_pDiffSub$R2)
predSetVec <-  c("items", "comb")


itemCol = c("#5C6BC0", "#42A5F5", "#29B6F6", "#26C6DA", "#26A69A", "#66BB6A", "#9CCC65", "#D4E157")
miCol <- c("#FFA726", "#FF7043", "#EF5350", "#EC407A")
scoreCol <- c("#000000")


for (iP in pVec) {
  for (iPredSet in predSetVec) {
    for (iMI in seq_along(miItemsVec)) {
      for (iCor in corGOvec) {
        
        print(paste(iP, iPredSet, iMI, iCor))
        
        if (iP == "4") {
          nMIitems = c(0, 1, 2)[iMI]
        } else {
          nMIitems = c(0, 1, 4)[iMI]
        }
        
        subTmp <- coef_pDiffSub[coef_pDiffSub$P == iP &
                                  # coef_pDiffSub$R2 == iR2 &
                                  coef_pDiffSub$predSet == iPredSet &
                                  coef_pDiffSub$miItems == nMIitems &
                                  coef_pDiffSub$corGO == iCor &
                                  coef_pDiffSub$coefficient != "intercept",] 
        
        # nMIitems <- ifelse(iMI == "P/2", as.numeric(iP)/2, as.numeric(iMI))
        
        if (iPredSet == "comb") {
          colVec <- c(itemCol[seq_len(as.numeric(iP) - nMIitems)], 
                      miCol[seq_len(nMIitems)],
                      scoreCol)  
        } else if (iPredSet != "comb") {
          colVec <- c(itemCol[seq_len(as.numeric(iP) - nMIitems)], 
                      miCol[seq_len(nMIitems)])
        }
        
        pTmp <- ggplot(subTmp, 
               aes(x = model, y = value, fill = coefficient, color = coefficient)) +
          coord_flip() + 
          geom_violinhalf(width=2.1, alpha = 0.5) +
          xlab("") +
          ylab("Regression Coefficients") +
          scale_color_manual(values = colVec) + 
          scale_fill_manual(values = colVec) +
          geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
          #geom_vline(aes(xintercept = 0), color = "black", alpha = 0.5) +
          facet_grid(R2 ~ MIvio)
        
        pTmp <- themeFunctionPaper(pTmp)
        
        ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCoefFolder, "/pCoef_P", 
                                          iP, "_MIitems", nMIitems, "_", iPredSet, "_cor", iCor, ".png"),
                        plot = pTmp,
                        width = 13.63,
                        height = 12.07,
                        units = "in")
        
      }
    }
  }
}

##### plot one presentation plot #####
nMIitems = c(0, 1, 2)[2]

subTmp <- coef_pDiffSub[coef_pDiffSub$P == 4 &
                          coef_pDiffSub$R2 == 0.8 &
                          coef_pDiffSub$predSet == "items" &
                          coef_pDiffSub$miItems == 1 &
                          coef_pDiffSub$corGO != 0.2 &
                          coef_pDiffSub$coefficient != "intercept",] 


colVec <- c(itemCol[seq_len(as.numeric(4) - nMIitems)], 
            miCol[seq_len(nMIitems)])

pTmp <- ggplot(subTmp, 
               aes(x = model, y = value, fill = coefficient, color = coefficient)) +
  coord_flip() + 
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients") +
  scale_color_manual(values = colVec) + 
  scale_fill_manual(values = colVec) +
  geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
  #geom_vline(aes(xintercept = 0), color = "black", alpha = 0.5) +
  facet_grid(corGO ~ MIvio)

pTmp <- themeFunctionPaper(pTmp)

ggplot2::ggsave(filename = paste0(plotFolder, "/present/pCoef_P4", "_MIitems", 
                                  nMIitems, "_items.png"),
                plot = pTmp,
                width = 13.63,
                height = 12.07,
                units = "in")

subTmp <- coef_pDiffSub[coef_pDiffSub$P == 4 &
                          coef_pDiffSub$R2 == 0.8 &
                          coef_pDiffSub$predSet == "comb" &
                          coef_pDiffSub$miItems == 1 &
                          coef_pDiffSub$corGO != 0.2 &
                          coef_pDiffSub$coefficient != "intercept",] 

colVec <- c(itemCol[seq_len(as.numeric(4) - nMIitems)], 
            miCol[seq_len(nMIitems)],
            scoreCol)  

pTmp <- ggplot(subTmp, 
               aes(x = model, y = value, fill = coefficient, color = coefficient)) +
  coord_flip() + 
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients") +
  scale_color_manual(values = colVec) + 
  scale_fill_manual(values = colVec) +
  geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
  #geom_vline(aes(xintercept = 0), color = "black", alpha = 0.5) +
  facet_grid(corGO ~ MIvio)

pTmp <- themeFunctionPaper(pTmp)

# # Item 4 is the one with the manipulation
# colVec <- c( "#006400", "#FF7F50", "#B22222") 
# ggplot(testSub[testSub$MIvio == "metric" &
#                   testSub$predSet == "items",], 
#        aes(x = X1, y = X4, group = N, 
#            color = N)) +
#   #geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
#   #geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
#   geom_point(size = 1.5, alpha = 0.5) +
#   #geom_line(alpha = 0.4, linewidth = 0.75) +
#   #scale_shape_manual(values = shapeVec) +
#   scale_color_manual(values = colVec) +
#   #scale_linetype_manual(values = lineTypeVec) +
#   #ylab(paste0("difference in ", iMT, " (items - LM score)")) +
#   facet_wrap(~ cAlpha + corGO)

# Does ENET or itemENET keep items with MI violations or specific effects more often 
#   than without manipulation? 
# maybe remove ENET and itemENET if alpha = 0, because ENET would not have been able 
#   to do variable selection!
# Is there a difference between the manipulations in item selection?
