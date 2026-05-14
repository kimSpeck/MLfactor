# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

library(tidyverse)
library(legendry)

resFolder <- "results"
depFolder <- "dependentMeasures"

plotFolder <- "plots"
createFolder(plotFolder)

plotPerformanceFolder <- "performance"
createFolder(paste0(plotFolder, "/", plotPerformanceFolder))

resFile <- paste0(resFolder, "/", depFolder, "/performPerSample_differences.rda")
load(resFile)

unique(performMeasures$model)

################################################################################
# difference plots (paper selection)
################################################################################
diffSub <- performMeasures[performMeasures$measure_type == "Rsquared_HOtest" &
                             performMeasures$model %in% c("diff_ENET", "diff_group_ENET"),]

rm(performMeasures)
gc()

diffSub$model <- sub("^.{5}", "", diffSub$model)
diffSub <- tidyr::separate(diffSub, model, into = c("predSet", "predMeth"), 
                           sep = "_", fill = "left")
diffSub$predSet <- ifelse(is.na(diffSub$predSet), "items", diffSub$predSet)

diffSub$R2 <- factor(diffSub$R2)

# nuance factor conditions
if (any(grepl("NF", diffSub$MIvio))){
  diffSub$corNF <- ifelse(grepl("NF", diffSub$MIvio),
                          substr(diffSub$MIvio, 3, nchar(as.character(diffSub$MIvio[1]))),
                          0)
  diffSub$MIvio <- ifelse(grepl("NF", diffSub$MIvio), "NF", diffSub$MIvio)
}

diffSub$MIvio <- factor(diffSub$MIvio, 
                        levels = c("metric", "scalarReal0.4", "scalarReal0.6", 
                                   "scalar", "specific"), 
                        labels = c("metric", "scalar (s)", "scalar (m)", "scalar", "nuance"))

diffSub$corGO_lab <- paste0("r[GY]==", diffSub$corGO)
diffSub$R2_lab    <- paste0("R[sim]^2==", diffSub$R2)
dodge_w <- 1

diffSub <- diffSub[diffSub$MIvio != "NF",]
diffSub$miItemsPlot <- factor(diffSub$miItemsPlot, 
                              levels = c("0", "1", "P/2"), 
                              labels = c("0", "1", "P/2"))

# effectively splitting labels with the regex sep = "[^[:alnum:]]+" (i.e., any non-alphanumeric run)
# any slash, percent sign, parentheses, spaces, Unicode punctuation, etc. will be treated as a split point
# we can't change sep in guide_axis_nested()
diffSub$x_nested <- interaction(diffSub$miItemsPlot, diffSub$MIvio, sep = "|", 
                                lex.order = FALSE)

levels(diffSub$x_nested)
diffSub$effectSize <- ifelse(diffSub$MIvio %in% c("metric", "nuance"), "", 
                             ifelse(diffSub$MIvio == "scalar", "xl", 
                                    ifelse(diffSub$MIvio == "scalar (s)", "s", "m")))

diffSub$MIvioPlot <- factor(diffSub$MIvio, 
                            levels = c("metric", "scalar (s)", "scalar (m)",
                                       "scalar", "nuance"), 
                            labels = c("metric", "scalar", "scalar", "scalar", "nuance"))
diffSub$x_nested2 <- interaction(diffSub$miItemsPlot, diffSub$MIvioPlot, diffSub$effectSize, 
                                 sep = "|", lex.order = FALSE)
levels(diffSub$x_nested2)
unique(diffSub$x_nested2)

varTvec <- unique(diffSub$varT)
predSetVec <- unique(diffSub$predSet)
nVec <- unique(diffSub$N)
cAlphaVec <- unique(diffSub$cAlpha)
Pvec <- unique(diffSub$P)

# # test it: items_1000_0.25_0.8
# iVarT <- 0.25
# iPredSet <- "items" # "group"
# iN <- 300
# iAlpha <- 0.8
# iP = 4

shapeVec <- c(3, 4)
colVecViridis <- c("#addd8e", "#1f9e89", "#1d4f91")
# colVecViridis <- c("#337538", "#0072B2", "#CC79A7")
# colVec <- c( "#006400", "#FF7F50", "#B22222") 

# range(diffSub[diffSub$MIvio != "NF" & 
#                 diffSub$predSet == "items" , ]$values)
# # -0.1668824  0.2610356
# 
# range(diffSub[diffSub$MIvio != "NF" & 
#                 diffSub$N == 300 & 
#                 diffSub$predSet == "items" , ]$values)
# # -0.02515503  0.26101686
# 
# range(diffSub[diffSub$MIvio != "NF" & 
#                 diffSub$predSet != "items" , ]$values)
# # -0.3472198  0.2715301
# 
# range(diffSub[diffSub$MIvio != "NF" & 
#                 diffSub$N == 300 & 
#                 diffSub$predSet != "items" , ]$values)
# # -0.2801884  0.2610176 

# identical scaling for all result plots (to a certain degree)
ylims_items    <- c(-0.17, 0.28)
ylims_itemsLargeN    <- c(-0.06, 0.28)
ylims_group <- c(-0.35, 0.28)

for (iVarT in varTvec) {
  for (iPredSet in predSetVec) {
    for (iN in nVec) {  
      for (iAlpha in cAlphaVec) {
        for (iP in Pvec) {
          pTmp <- ggplot(diffSub[diffSub$varT == iVarT &
                                   diffSub$N == iN &
                                   diffSub$P == iP &
                                   diffSub$predSet == iPredSet &
                                   diffSub$cAlpha == iAlpha,], 
                         aes(x = x_nested, y = values,
                             # group = interaction(P, MIvio, miItemsPlot), 
                             group = interaction(MIvio, miItemsPlot), 
                             color = miItemsPlot, fill = miItemsPlot)) +
            geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
            geom_violin(
              #aes(linetype = P), 
              position = position_dodge(width = dodge_w),
              width = 1.25,
              trim = TRUE,
              alpha = 0.90,
              linewidth = 0.65
            ) +
            stat_summary(#aes(shape = P),
                         fun = mean,
                         shape = 4,
                         geom = "point",
                         color = "black",
                         position = position_dodge(width = dodge_w),
                         size = 1.7,
                         stroke = 0.5
            ) +
            scale_color_manual(values = colVecViridis, name = "MI items") +
            scale_fill_manual(values = colVecViridis, name = "MI items") +
            scale_shape_manual(values = shapeVec, name = "P") +
            scale_x_discrete(
              breaks = levels(diffSub$x_nested),
              labels = sub("\\|.*$", "", levels(diffSub$x_nested)),
              guide = guide_axis_nested(
                key = key_range_manual(
                  start = c(1, 4, 6, 8, 11, 14, 17, 19, 21, 24, 27, 30, 32, 34, 37),
                  end   = c(3, 5, 7, 10, 13, 16, 18, 20, 23, 26, 29, 31, 33, 36, 39),
                  name  = c("metric", "\ns", "scalar\nm", "\nxl", "nuance",
                            "metric", "\ns", "scalar\nm", "\nxl", "nuance",
                            "metric", "\ns", "scalar\nm", "\nxl", "nuance")
                )
              )
            ) +
            # scale_x_discrete(
            #   breaks = levels(diffSub$x_nested2),
            #   labels = sub("\\|.*$", "", levels(diffSub$x_nested2)),
            #   guide = guide_axis_nested(
            #     key = list(
            #       # Level 1: metric / scalar / nuance
            #       legendry::key_range_manual(
            #       start = c(1, 4, 11, 14, 17, 24, 27, 30, 37),
            #       end   = c(3, 10, 13, 16, 23, 26, 29, 36, 39),
            #       name  = c("metric", "scalar", "nuance", 
            #                 "metric", "scalar", "nuance", 
            #                 "metric", "scalar", "nuance")
            #       ), 
            #       # Level 2: scalar subtypes (s / m / xl) inside the scalar block
            #       legendry::key_range_manual(
            #         start = c(4, 6, 8, 17, 19, 21, 30, 32, 34),
            #         end   = c(5, 7, 10, 18, 20, 23, 31, 33, 36),
            #         name  = c("s", "m", "xl", 
            #                   "s", "m", "xl", 
            #                   "s", "m", "xl")
            #       )
            #     )
            #   )
            # ) +
            guides(linetype = "none", color = "none", fill = "none") +
            labs(y = expression(Delta*R^2~"(items - score)"), x = "source of item-specific variance") +
            facet_grid(corGO_lab ~ R2_lab, labeller = label_parsed)
          
          #pTmp
          
          if (iPredSet == "items") {
            if (iN > 100) {
              ylims <- ylims_itemsLargeN
            } else {
              ylims <- ylims_items
            }
          } else {
            ylims <- ylims_group
          }
          
          pTmp <- pTmp + coord_cartesian(ylim = ylims) +
            scale_y_continuous(
              breaks = seq(floor(ylims[1]*5)/5, ceiling(ylims[2]*5)/5, by = 0.2)
            )
          
          (pTmp <- themeFunctionPaper(pTmp))
          
          plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                             "/pDiffR2_Ho_", iPredSet, "_N", iN, "_varT", iVarT, 
                             "_alpha", iAlpha, "_P", iP, "_violinPlots.png")
          
          ggplot2::ggsave(filename = plotName,
                          plot = pTmp,
                          width = 13.63,
                          height = 9.60,
                          units = "in")
        }
      }
    }
  }  
} 

################################################################################
# difference plots (exploratory)
################################################################################
# items - combination (items + scores) for ENET
# does not show differences in most conditions with larger uncertainty in differences
# P = 8 shows large uncertainty with averages across replications showing slightly 
#   positive values
# for N = 1000: P = 8 with one manipulated item still indicates predictive gains 
#   across all sources of item-specific variance <- why?
load(resFile)

diffSub <- performMeasures[performMeasures$measure_type == "Rsquared_HOtest" &
                             performMeasures$model == "diff_comb_ENET",]

rm(performMeasures)
gc()

diffSub$model <- "ENET"
diffSub$predSet <- "comb"

diffSub$R2 <- factor(diffSub$R2)
diffSub$corNF <- ifelse(grepl("NF", diffSub$MIvio), 
                        substr(diffSub$MIvio, 3, nchar(as.character(diffSub$MIvio[1]))), 
                        0)
diffSub$MIvio <- ifelse(grepl("NF", diffSub$MIvio), "NF", diffSub$MIvio)
diffSub$MIvio <- factor(diffSub$MIvio, 
                        levels = c("metric", "scalar", "specific", "NF"), 
                        labels = c("metric", "scalar", "nuance", "NF"))

diffSub$corGO_lab <- paste0("r[GY]==", diffSub$corGO)
diffSub$R2_lab    <- paste0("R[sim]^2==", diffSub$R2)
dodge_w <- 1

diffSub <- diffSub[diffSub$MIvio != "NF",]
diffSub$miItemsPlot <- factor(diffSub$miItemsPlot, 
                              levels = c("0", "1", "P/2"), 
                              labels = c("0", "1", "P/2"))

# effectively splitting labels with the regex sep = "[^[:alnum:]]+" (i.e., any non-alphanumeric run)
# any slash, percent sign, parentheses, spaces, Unicode punctuation, etc. will be treated as a split point
# we can't change sep in guide_axis_nested()
diffSub$x_nested <- interaction(diffSub$miItemsPlot, diffSub$MIvio, sep = "|", 
                                lex.order = FALSE)
levels(diffSub$x_nested)

varTvec <- unique(diffSub$varT)
predSetVec <- unique(diffSub$predSet)
nVec <- unique(diffSub$N)
cAlphaVec <- unique(diffSub$cAlpha)

# # test it: items_1000_0.25_0.8
# iVarT <- 0.25
# iPredSet <- "comb" # "group"
# iN <- 1000
# iAlpha <- 0.8

range(diffSub[diffSub$MIvio != "NF", ]$values)
# -0.1428541  0.2945100

range(diffSub[diffSub$MIvio != "NF" & 
                diffSub$N == 300, ]$values)
# -0.0302232  0.1165763

ylims <- c(-0.15, 0.3)

for (iVarT in varTvec) {
  for (iN in nVec) {  
    for (iAlpha in cAlphaVec) {
      pTmp <- ggplot(diffSub[diffSub$varT == iVarT &
                               diffSub$MIvio != "NF" &
                               diffSub$N == iN &
                               diffSub$cAlpha == iAlpha,], 
                     # aes(x = MIvio, y = values, 
                     aes(x = x_nested, y = values, 
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
        scale_x_discrete(
          breaks = levels(diffSub$x_nested),
          labels = sub("\\|.*$", "", levels(diffSub$x_nested)),
          guide = guide_axis_nested(
            key = key_range_manual(
              start = c(1, 4, 7, 10, 13, 16, 19, 22, 25),
              end   = c(3, 6, 9, 12, 15, 18, 21, 24, 27),
              name  = c("metric", "scalar", "nuance", 
                        "metric", "scalar", "nuance", 
                        "metric", "scalar", "nuance")
            )
          )
        ) +
        guides(linetype = "none", color = "none", fill = "none") +
        labs(y = expression(Delta*R^2~"(items - comb)"), x = "source of item-specific variance") +
        facet_grid(corGO_lab ~ R2_lab, labeller = label_parsed)
      
      pTmp <- pTmp + coord_cartesian(ylim = ylims) +
        scale_y_continuous(
          breaks = seq(floor(ylims[1]*5)/5, ceiling(ylims[2]*5)/5, by = 0.2)
        )
      
      (pTmp <- themeFunctionPaper(pTmp))
      
      plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                         "/pDiffR2_Ho_comb", "_N", iN, "_varT", iVarT, 
                         "_alpha", iAlpha, "_violinPlots.png")
      
      ggplot2::ggsave(filename = plotName,
                      plot = pTmp,
                      width = 13.63,
                      height = 9.60,
                      units = "in")
    }
  }
}  
