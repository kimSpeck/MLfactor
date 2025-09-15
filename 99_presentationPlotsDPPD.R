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

resFile <- paste0(resFolder, "/", depFolder, "/performStats.rda")
load(resFile)

################################################################################
# plot test plot
################################################################################
themeFunctionPaper <- function(plotObject, guides = T){
  pTmp <- plotObject + theme(
    panel.grid.major = element_line(linewidth = 0.15, linetype = 'solid', color = "lightgrey"), 
    panel.grid.minor = element_line(linewidth = 0.1, linetype = 'solid', color = "lightgrey"),
    panel.background = element_rect(color = "white", fill = "white"),
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    #legend.position = c(.85, .39), 
    legend.position = "bottom", 
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    legend.key.width = unit(2, "cm"),
    legend.box = "horizontal")
  if (guides == F) {
    pTmp <- pTmp + guides(color = "none", shape = "none", linetype = "none")
  }
  return(pTmp)
}

################################################################################
# plot specific effect
################################################################################
# plot only LM results, first
# R² = 0.8, corGO = 0 (aber essentially egal), varT = 0.25
# ... p = 4 vs. p = 8 (without correlation) and cAlpha = 0.8 vs. cAlpha = 0.6

# comb, items, score
shapeVec <- c(13, 20, 3)
colVec <- c("#003399", "#006400", "#006400") 
lineTypeVec <- c("solid", "dashed", "dotted")

diffSub <- performStats[performStats$model %in% c("diff_LM", "diff_ENET", "diff_Ridge"),]
# str(diffSub)
diffSub$model <- factor(diffSub$model, 
                        levels = c("diff_LM", "diff_Ridge", "diff_ENET"), 
                        labels = c("LM", "Ridge", "ENET"))

diffSpecific <- diffSub[diffSub$MIvio == "specific" &
                          diffSub$model == "LM" &
                          diffSub$R2 == 0.8 &
                          diffSub$corGO == 0 &
                          diffSub$varT == 0.25, ]

cAlphaVec = unique(diffSpecific$cAlpha)
pVec = unique(diffSpecific$P)

for (iP in pVec) {
  for (iCA in cAlphaVec) {
    pSpecificTmp <- ggplot(diffSpecific[diffSpecific$cAlpha == iCA &
                                          diffSpecific$P == iP &
                                          diffSpecific$measure_type == "Rsquared_HOtest",], 
                           aes(x = N, y = M, group = miItemsPlot, color = miItemsPlot, linetype = miItemsPlot)) +
      geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
      geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.5, width = 0.15) + 
      geom_point(size = 5) +
      geom_line(alpha = 0.8, linewidth = 1.66) +
      # scale_shape_manual(values = shapeVec) +
      scale_color_manual(values = colVec) +
      scale_linetype_manual(values = lineTypeVec) +
      scale_y_continuous(limits = c(-0.075, 0.2),
                         breaks = seq(0, 0.2, by = 0.1)) +
      ylab("") + xlab("") +
      guides(color = "none", 
             linetype = "none")
    
    (pSpecificTmp <- themeFunctionPaper(pSpecificTmp))
    
    ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                      "/pSpecific_Rsquared_HOtest", 
                                      "_cAlpha", iCA, "_p", iP, ".png"),
                    plot = pSpecificTmp,
                    width = 11.56,
                    height = 9.77,
                    units = "in")
    
  }
}

pSpecific2x2 <- ggplot(diffSpecific[diffSpecific$measure_type == "Rsquared_HOtest",], 
       aes(x = N, y = M, group = miItemsPlot, color = miItemsPlot, linetype = miItemsPlot)) +
  geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.5, width = 0.15) + 
  geom_point(size = 3) +
  geom_line(alpha = 0.8, linewidth = 1) +
  # scale_shape_manual(values = shapeVec) +
  scale_color_manual(values = colVec) +
  scale_linetype_manual(values = lineTypeVec) +
  scale_y_continuous(limits = c(-0.075, 0.2),
                     breaks = seq(0, 0.2, by = 0.1)) +
  ylab("") + xlab("") +
  facet_grid(P ~ cAlpha) +
  guides(color = "none", 
         linetype = "none")

(pSpecific2x2 <- themeFunctionPaper(pSpecific2x2))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                  "/pSpecific_Rsquared_HOtest_cAlphaxP.png"),
                plot = pSpecific2x2,
                width = 11.56,
                height = 9.77,
                units = "in")

################################################################################
# plot scalar effect
################################################################################
# plot only LM results, first
# R² = 0.8, corGO = 0 (aber essentially egal), varT = 0.25
# ... p = 4 vs. p = 8 (without correlation) and cAlpha = 0.8 vs. cAlpha = 0.6

colVec <- c("#003399", "#c9211e", "#c9211e") 
lineTypeVec <- c("solid", "dashed", "dotted")

diffScalar <- diffSub[diffSub$MIvio == "scalar" &
                          diffSub$model == "LM" &
                          diffSub$R2 == 0.8 &
                          diffSub$varT == 0.25, ]

diffScalar$miItemsPlot <- factor(diffScalar$miItemsPlot, 
                                 levels = c("0", "1", "P/2"))

cAlphaVec = unique(diffScalar$cAlpha)
pVec = unique(diffScalar$P)
corVec = unique(diffScalar$corGO)

for (iP in pVec) {
  for (iCA in cAlphaVec) {
    for (iCor in corVec) {
      pScalarTmp <- ggplot(diffScalar[diffScalar$cAlpha == iCA &
                                        diffScalar$P == iP &
                                        diffScalar$corGO == iCor &
                                        diffScalar$measure_type == "Rsquared_HOtest",], 
                           aes(x = N, y = M, group = miItemsPlot, color = miItemsPlot, linetype = miItemsPlot)) +
        geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
        geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.5, width = 0.15) + 
        geom_point(size = 5) +
        geom_line(alpha = 0.8, linewidth = 1.66) +
        # scale_shape_manual(values = shapeVec) +
        scale_color_manual(values = colVec) +
        scale_linetype_manual(values = lineTypeVec) +
        scale_y_continuous(limits = c(-0.075, 0.2),
                           breaks = seq(0, 0.2, by = 0.1)) +
        ylab("") + xlab("") +
        guides(color = "none", 
               linetype = "none")
      
      (pScalarTmp <- themeFunctionPaper(pScalarTmp))
      
      ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                        "/pScalar_Rsquared_HOtest", 
                                        "_cAlpha", iCA, "_p", iP, "_cor", iCor, ".png"),
                      plot = pScalarTmp,
                      width = 11.56,
                      height = 9.77,
                      units = "in")
      
    }
  }
}


for (iP in pVec) {
  pScalar2x3 <- ggplot(diffScalar[diffScalar$P == iP &
                                    diffScalar$measure_type == "Rsquared_HOtest",], 
                       aes(x = N, y = M, group = miItemsPlot, color = miItemsPlot, linetype = miItemsPlot)) +
    geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
    geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.5, width = 0.15) + 
    geom_point(size = 3) +
    geom_line(alpha = 0.8, linewidth = 1) +
    # scale_shape_manual(values = shapeVec) +
    scale_color_manual(values = colVec) +
    scale_linetype_manual(values = lineTypeVec) +
    scale_y_continuous(limits = c(-0.075, 0.2),
                       breaks = seq(0, 0.2, by = 0.1)) +
    ylab("") + xlab("") +
    facet_grid(corGO ~ cAlpha) +
    guides(color = "none", 
           linetype = "none")
  
  (pScalar2x3 <- themeFunctionPaper(pScalar2x3))
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/", plotPresentFolder, 
                                    "/pScalar_Rsquared_HOtest_cAlphaxCor_P", iP,".png"),
                  plot = pScalar2x3,
                  width = 11.56,
                  height = 9.77,
                  units = "in")
}
