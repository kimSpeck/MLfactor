# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

library(tidyverse)

resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/performStats.rda")
load(resFile)

plotFolder <- "plots"
createFolder(plotFolder)
plotPerformanceFolder <- "performance"
createFolder(paste0(plotFolder, "/", plotPerformanceFolder))

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

# comb, items, score
shapeVec <- c(13, 20, 3)
colVec <- c( "#006400", "#FF7F50", "#B22222") 
lineTypeVec <- c("solid", "dotted", "dashed")

################################################################################
# difference measures on train, cvTest and hoTest
################################################################################
diffSub <- performStats[performStats$model %in% c("diff_LM", "diff_ENET", "diff_Ridge"),]
# str(diffSub)
diffSub$model <- factor(diffSub$model, 
                        levels = c("diff_LM", "diff_Ridge", "diff_ENET"), 
                        labels = c("LM", "Ridge", "ENET"))

# iterate over these vectors ...  
MIvioVec <- unique(diffSub$MIvio)
varTvec <- unique(diffSub$varT)
measureTypeVec <- unique(diffSub$measure_type)

# plot differences in all measure types, for all MIviolations and every varT condition
for (iMIvio in MIvioVec) {
  for (iVarT in varTvec) {
    for (iMT in measureTypeVec) {
      pTmp <- ggplot(diffSub[diffSub$MIvio == iMIvio & 
                               diffSub$varT == iVarT & 
                               diffSub$measure_type == iMT,], 
             aes(x = N, y = M, group = interaction(miItemsPlot, model), 
                 color = miItemsPlot,
                 # shape = predSet,
                 linetype = model)) +
        geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
        geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
        geom_point(size = 3) +
        geom_line(alpha = 0.4, linewidth = 0.75) +
        scale_shape_manual(values = shapeVec) +
        scale_color_manual(values = colVec) +
        scale_linetype_manual(values = lineTypeVec) +
        ylab(paste0("difference in ", iMT, " (items - LM score)")) +
        facet_grid(P + corGO ~ interaction(MIvio, cAlpha, R2))
      
      pTmp <- themeFunctionPaper(pTmp)
      
      plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                         "/pDiff_", iMT, "_", iMIvio, "_varT", iVarT, ".png")
      
      ggplot2::ggsave(filename = plotName,
                      plot = pTmp,
                      width = 13.63,
                      height = 9.60,
                      units = "in")
    }
  }
}

# in general: differences are more pronounced in varT = 0.25 instead of varT = 0.5
# the differences between items and scores do not depend on hoTest vs. Train performance data; 
#   no bias or dramatically decreased standard deviations due to simulating only one test sample
#   overfitting in training is apparent but other than that, the results in terms of 
#     difference measure i.e., effect of the MI violation/specific item effect is the same

# plot differences for
# ... measure_type = R2_HOtest 
# ... varT = 0.25
# ... and all MIviolations in one plot
# ... plot R2 separately


diffSub$MIvio <- factor(diffSub$MIvio, 
                        levels = c("specific", "metric", "scalar"))

# plot every combination of varT and R2 with MI/specific effect manipulation
#   directly next to each other
varTvec <- unique(diffSub$varT)
R2vec <- unique(diffSub$R2)

for (iVarT in varTvec) {
  for (iR2 in R2vec) {
    pTmp <- ggplot(diffSub[diffSub$varT == iVarT & 
                             diffSub$R2 == iR2 &
                             diffSub$measure_type == "Rsquared_HOtest",], 
                   aes(x = N, y = M, group = interaction(miItemsPlot, model), 
                       color = miItemsPlot,
                       # shape = predSet,
                       linetype = model)) +
      geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
      geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
      geom_point(size = 3) +
      geom_line(alpha = 0.4, linewidth = 0.75) +
      scale_shape_manual(values = shapeVec) +
      scale_color_manual(values = colVec) +
      scale_linetype_manual(values = lineTypeVec) +
      ylab(paste0("difference in ", iMT)) +
      facet_grid(P + corGO ~ interaction(MIvio, cAlpha))
    
    pTmp <- themeFunctionPaper(pTmp)
    
    plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                       "/pDiff_Rsquared_HOtest", "_R2", iR2, "_varT", iVarT, ".png")
    
    ggplot2::ggsave(filename = plotName,
                    plot = pTmp,
                    width = 13.63,
                    height = 9.60,
                    units = "in")
    
  }
}

# cAlpha is essentially irrelevant; sieht aus wie gespiegelte Spalten! 
# R2 hat die größte Auswirkung auf den Effekt der metrischen MI violation, weil 
#   die Manipulation nur indirekt über den Faktor und dessen Effekt auf das Outcome wirkt
# metric MI violation as indirect effect on the outcome via the factor itself
# scalar MI violation as indirect effect on the outcome via correlated grouping variable 

# plot central results graphic for R2
# ... varT = 0.25
# ... R2 = 0.8 
# ... cAlpha = 0.8 

r2_hoTest_R20.8_varT0.25_cAlpha0.8 <- ggplot(diffSub[diffSub$varT == 0.25 & 
                                                       diffSub$R2 == 0.8 &
                                                       diffSub$cAlpha == 0.8 &
                                                       diffSub$measure_type == "Rsquared_HOtest",], 
               aes(x = N, y = M, group = interaction(miItemsPlot, model), 
                   color = miItemsPlot,
                   # shape = predSet,
                   linetype = model)) +
  geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
  geom_point(size = 3) +
  geom_line(alpha = 0.4, linewidth = 0.75) +
  scale_shape_manual(values = shapeVec) +
  scale_color_manual(values = colVec) +
  scale_linetype_manual(values = lineTypeVec) +
  # scale_y_continuous(limits = c(-0.1, 0.2), 
  #                    breaks = seq(0, 0.2, by = 0.1)) +
  ylab("difference in R²") +
  facet_grid(P + corGO ~ MIvio)

(r2_hoTest_R20.8_varT0.25_cAlpha0.8 <- themeFunctionPaper(r2_hoTest_R20.8_varT0.25_cAlpha0.8))

plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                   "/pDiff_Rsquared_HOtest_R20.8_varT0.25_cAlpha0.8.png")

ggplot2::ggsave(filename = plotName,
                plot = r2_hoTest_R20.8_varT0.25_cAlpha0.8,
                width = 11.25,
                height = 11.02,
                units = "in")

################################################################################
# 
################################################################################

# r2_hoTestSub <- unite(data = r2_hoTestSub, col = "predSet_model", c(predSet, model), sep = "_")
# r2_hoTestSub_wide <- tidyr::pivot_wider(r2_hoTestSub, names_from = predSet_model,
#                                         values_from = values)
# r2_hoTestSub_wide$R2_diff_LM <- r2_hoTestSub_wide$items_LM - r2_hoTestSub_wide$score_LM
# r2_hoTestSub_wide$R2_diff_Ridge <- r2_hoTestSub_wide$items_Ridge - r2_hoTestSub_wide$score_LM
# r2_hoTestSub_wide$R2_diff_ENET <- r2_hoTestSub_wide$items_ENET - r2_hoTestSub_wide$score_LM
# r2_hoTestSub_wide$R2_relDiff_LM <- r2_hoTestSub_wide$R2_diff_LM / r2_hoTestSub_wide$R2
# 
# r2_hoTest_plotDiff <- r2_hoTestSub_wide[r2_hoTestSub_wide$varT == 0.25,] %>% 
#   group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N) %>% 
#   summarize(M_diffLM = mean(R2_diff_LM, na.rm = T),
#             SD_diffLM = sd(R2_diff_LM, na.rm = T),
#             M_diffRidge = mean(R2_diff_Ridge, na.rm = T),
#             SD_diffRidge = sd(R2_diff_Ridge, na.rm = T),
#             M_diffENET = mean(R2_diff_ENET, na.rm = T),
#             SD_diffENET = sd(R2_diff_ENET, na.rm = T),
#             M_relDiffLM = mean(R2_relDiff_LM, na.rm = T),
#             SD_relDiffLM = sd(R2_relDiff_LM, na.rm = T)) 
# r2_hoTest_plotDiff <- tidyr::pivot_longer(r2_hoTest_plotDiff, 
#                                           cols = colnames(r2_hoTest_plotDiff)[grep(pattern = "iff", colnames(r2_hoTest_plotDiff))],
#                                           names_to = "model", values_to = "values")
# r2_hoTest_plotDiff <- separate(r2_hoTest_plotDiff, model, into = c("what", "model"), sep = "_")
# r2_hoTest_plotDiff <- tidyr::pivot_wider(r2_hoTest_plotDiff, names_from = what, 
#                                          values_from = values)
# 
# r2_hoTest_plotDiff$model <- factor(r2_hoTest_plotDiff$model, 
#                                    levels = c("diffLM", "diffRidge", "diffENET", "relDiffLM"), 
#                                    labels = c("LM", "Ridge", "ENET", "relDiffLM"))

# # R2 on holdout test sample
# # -> one sample on which every model is tested
# r2_hoTestSub <- performMeasures[performMeasures$measure_type == "Rsquared_HOtest",]
# r2_hoTestSub <- performMeasures[performMeasures$measure_type == "Rsquared_train",]
# r2_hoTestSub$miItemsPlot <- ifelse(r2_hoTestSub$miItems %in% c(2, 4), "P/2", r2_hoTestSub$miItems)
# r2_hoTestSub$N <- factor(r2_hoTestSub$N, levels = c(100, 300, 1000))
# 
# r2_hoTestStats_varT0.25 <- r2_hoTestSub[r2_hoTestSub$varT == 0.25,] %>% 
#   group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
#   summarize(M = mean(values, na.rm = T),
#             SD = sd(values, na.rm = T))
# 
# R2_hoTest_var0.25_LM <- ggplot(r2_hoTestStats_varT0.25[r2_hoTestStats_varT0.25$model %in% c("LM"),], 
#                                aes(x = N, y = M, group = interaction(miItemsPlot, model, predSet), 
#                                    color = miItemsPlot,
#                                    shape = predSet,
#                                    linetype = model)) +
#   geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
#   geom_point(size = 3) +
#   geom_line(alpha = 0.4, linewidth = 0.75) +
#   scale_shape_manual(values = shapeVec) +
#   scale_color_manual(values = colVec) +
#   facet_grid(cAlpha + corGO ~ interaction(MIvio, P, R2))
# themeFunctionPaper(R2_hoTest_var0.25_LM)
# 
# 
# 
# # pR2_hoTest_relDiff <- ggplot(r2_hoTest_plotDiff, 
# #                           aes(x = N, y = M_relDiff, group = interaction(miItemsPlot, model), 
# #                               color = miItemsPlot,
# #                               # shape = predSet,
# #                               linetype = model)) +
# #   geom_errorbar(aes(ymin = M_relDiff - 2*SD_relDiff, ymax = M_relDiff + 2*SD_relDiff), alpha = 0.2, width = 0.1) + 
# #   geom_point(size = 3) +
# #   geom_line(alpha = 0.4, linewidth = 0.75) +
# #   scale_shape_manual(values = shapeVec) +
# #   scale_color_manual(values = colVec) +
# #   facet_grid(cAlpha + corGO ~ interaction(MIvio, P, R2))
# # themeFunctionPaper(pR2_hoTest_relDiff)
# 
# ggplot(r2_hoTestStats_varT0.25[which((
#   r2_hoTestStats_varT0.25$model == "LM" & r2_hoTestStats_varT0.25$predSet == "score") |
#     (r2_hoTestStats_varT0.25$model %in% c("ENET", "itemENET"))),], 
#   aes(x = N, y = M, group = interaction(miItemsPlot, model, predSet), 
#       color = miItemsPlot,
#       shape = predSet,
#       linetype = model)) +
#   geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
#   geom_point(size = 2) +
#   geom_line(alpha = 0.6) +
#   scale_shape_manual(values = shapeVec) +
#   facet_grid(corGO + cAlpha ~ interaction(MIvio, P, R2))
# 
# # R2 for training sample (i.e., full training Data)
# r2_trainSub <- performMeasures[performMeasures$measure_type == "Rsquared_train",]
# r2_trainSub$miItemsPlot <- ifelse(r2_trainSub$miItems %in% c(2, 4), "P/2", r2_trainSub$miItems)
# r2_trainSub$N <- factor(r2_trainSub$N, levels = c(100, 300, 1000))
# 
# r2_trainStats_varT0.25 <- r2_trainSub[r2_trainSub$varT == 0.25,] %>% 
#   group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
#   summarize(M = mean(values, na.rm = T),
#             SD = sd(values, na.rm = T))
# 
# ggplot(r2_trainStats_varT0.25[r2_trainStats_varT0.25$model %in% c("LM"),], 
#        aes(x = N, y = M, group = interaction(miItemsPlot, model, predSet), 
#            color = miItemsPlot,
#            shape = predSet,
#            linetype = model)) +
#   geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
#   geom_point(size = 2) +
#   geom_line(alpha = 0.6) +
#   scale_shape_manual(values = shapeVec) +
#   facet_grid(corGO + cAlpha ~ interaction(MIvio, P, R2))
