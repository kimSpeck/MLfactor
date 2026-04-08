# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

library(tidyverse)

resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/hyperPerSample.rda")
load(resFile)

plotFolder <- "plots"
createFolder(plotFolder)
plotHyperFolder <- "hyperPlots"
createFolder(paste0(plotFolder, "/", plotHyperFolder))

hyperSelection$miItemsPlot <- ifelse(hyperSelection$miItems %in% c(2, 4), "P/2", hyperSelection$miItems)
hyperSelection$N <- factor(hyperSelection$N, levels = c(100, 300, 1000))

################################################################################
# How are tuned alpha values distributed? 
################################################################################
# this illustrates tuned alpha values in every replication via a violin plot
# (plot-design is equivalent to Delta R^2 plots of the results section)
plotAlpha <- hyperSelection[hyperSelection$model == "ENET" & 
                              hyperSelection$predSet == "items", ]

plotAlpha <- plotAlpha[which(!grepl("^NF0\\.", plotAlpha$MIvio)), ]
varTvec <- unique(plotAlpha$varT)
nVec <- unique(plotAlpha$N)
cAlphaVec <- unique(plotAlpha$cAlpha)
plotAlpha$MIvio <- factor(plotAlpha$MIvio, 
                          levels = c("metric", "scalar", "specific"), 
                          labels = c("metric", "scalar", "nuance"))

plotAlpha$corGO_lab <- paste0("r[GY]==", plotAlpha$corGO)
plotAlpha$R2_lab    <- paste0("R[sim]^2==", plotAlpha$R2)

shapeVec <- c(3, 4)
colVecViridis <- c("#337538", "#0072B2", "#CC79A7")
dodge_w <- 1

for (iVarT in varTvec) {
  for (iN in nVec) {  
    for (iAlpha in cAlphaVec) {
      pTmp <- ggplot(plotAlpha[plotAlpha$varT == iVarT &
                                 plotAlpha$N == iN &
                                 plotAlpha$cAlpha == iAlpha,], 
                     aes(x = MIvio, y = alpha, 
                         group = interaction(P, MIvio, miItemsPlot), 
                         color = miItemsPlot, 
                         fill = miItemsPlot)) +
        geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
        geom_violin(
          aes(linetype = P), 
          position = position_dodge(width = dodge_w),
          width = 1.25,
          trim = TRUE,
          alpha = 0.80,
          linewidth = 0.65
        ) +
        stat_summary(aes(shape = P),
                     fun = mean,
                     geom = "point",
                     color = "black",
                     position = position_dodge(width = dodge_w),
                     size = 1.7,
                     stroke = 0.5
        ) +
        scale_color_manual(values = colVecViridis, name = "MI items") +
        scale_fill_manual(values = colVecViridis, name = "MI items") +
        scale_shape_manual(values = shapeVec, name = "P") +
        guides(linetype = "none") +
        labs(y = expression("tuned"~alpha), x = "source of item-specific variance") +
        facet_grid(corGO_lab ~ R2_lab, labeller = label_parsed)
      
      pTmp <- themeFunctionPaper(pTmp) + 
        coord_cartesian(ylim = c(0, 1)) +
        scale_y_continuous(
          breaks = seq(0, 1, by = 0.5))
      
      plotName <- paste0(plotFolder, "/", plotHyperFolder, 
                         "/alpha0_N", iN, "_varT", iVarT, 
                         "_cAlpha", iAlpha, "_violinPlots.png")
      
      ggplot2::ggsave(filename = plotName,
                      plot = pTmp,
                      width = 13.63,
                      height = 9.60,
                      units = "in")
    }
  }
}  
################################################################################
# How often is ENET essentially a ridge regression, due to alpha == 0?
################################################################################
# this illustrates the number of replications that chose exactly alpha = 0!
# essentially: How many replications estimated Ridge regression due to tuning?

hyperSelection$alpha <- factor(hyperSelection$alpha)

countAlpha <- hyperSelection %>% 
  group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
  count(alpha, .drop = FALSE) 

# # check
# all(countAlpha[countAlpha$model == "Ridge", "n"] == 100)

# remove ridge (no information in this data)
countAlpha <- countAlpha[countAlpha$model != "Ridge", ]

# # lots of conditions with all samples alpha = 0
# View(countAlpha[which(countAlpha$n == 100), ])

# # no chance to see the pattern or any kind of numbers here...
# View(countAlpha[which(countAlpha$predSet == "items" & 
#                         countAlpha$alpha == 0), ])

plotRidge <- countAlpha[which(!grepl("^NF0\\.", countAlpha$MIvio) &
                                countAlpha$predSet == "items" & 
                                countAlpha$alpha == 0), ]
varTvec <- unique(plotRidge$varT)
nVec <- unique(plotRidge$N)
cAlphaVec <- unique(plotRidge$cAlpha)
plotRidge$nPC <- plotRidge$n / params$sim$nTrain
plotRidge$MIvio <- factor(plotRidge$MIvio, 
                        levels = c("metric", "scalar", "specific"), 
                        labels = c("metric", "scalar", "nuance"))
shapeVec <- c(3, 4)
colVec <- c("#337538", "#0072B2", "#CC79A7") 
lineTypeVec <- c("dashed", "solid")

range(plotRidge$nPC)

# # test it: items_1000_0.25_0.8
# iVarT <- 0.25
# iN <- 1000
# iAlpha <- 0.8
for (iVarT in varTvec) {
  for (iN in nVec) {  
    for (iAlpha in cAlphaVec) {
      pTmp <- ggplot(plotRidge[plotRidge$varT == iVarT &
                                 plotRidge$N == iN &
                                 plotRidge$cAlpha == iAlpha,], 
                     aes(x = MIvio, y = nPC, 
                         group = interaction(miItemsPlot, P), 
                         color = miItemsPlot,
                         shape = P, 
                         fill = P,
                         linetype = P)) +
        geom_point(aes(fill = P), size = 4, alpha = 0.7) +
        geom_line(alpha = 0.7, linewidth = 1) +
        scale_shape_manual(values = shapeVec) +
        scale_color_manual(values = colVec) +
        scale_linetype_manual(values = lineTypeVec) +
        ylab(paste0("R² difference (items - score)")) +
        facet_grid(corGO ~ R2)
      
      pTmp <- themeFunctionPaper(pTmp)
      
      plotName <- paste0(plotFolder, "/", plotHyperFolder, 
                         "/alpha0_N", iN, "_varT", iVarT, 
                         "_cAlpha", iAlpha, "_linePlots.png")
      
      ggplot2::ggsave(filename = plotName,
                      plot = pTmp,
                      width = 13.63,
                      height = 9.60,
                      units = "in")
    }
  }
}  

# only look at lambda if there are any specific research questions involved, etc. 
# ################################################################################
# # How strong is penalty? lambda values? 
# ################################################################################
# hyperSelection$lambda <- factor(hyperSelection$lambda)
# 
# countLambda <- hyperSelection %>% 
#   group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
#   count(lambda, .drop = FALSE) 
# 
# countLambda
# 
# ggplot(countLambda[countLambda$model == "ENET" & 
#                      countLambda$varT == 0.25, ], 
#        aes(x = factor(lambda), y = n, fill = interaction(miItemsPlot, MIvio))) +
#   geom_col(position = position_dodge(width = 0.7)) +
#   scale_y_continuous(limits = c(-10, 100), breaks = seq(0, 100, 50)) +
#   facet_grid(cAlpha + corGO ~ N + P + R2) +
#   ylab("Count per level")
