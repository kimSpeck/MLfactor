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
# descriptive statistics
################################################################################
dataSub <- performMeasures[performMeasures$measure_type == "Rsquared_HOtest" &
                  performMeasures$model == "diff_ENET",]

range(dataSub$values)

stats <- dataSub %>% 
  group_by(MIvio, P, cAlpha, varT, corGO, R2, N, miItemsPlot) %>% 
  summarize(M = mean(values, na.rm = T), 
            SD = sd(values, na.rm = T), 
            min = min(values, na.rm = T), 
            max = max(values, na.rm = T)) %>% 
  arrange(desc(M))

stats$varS <- ifelse(stats$varT == 0.25, 0.5, 0.25)
stats <- stats[!(stats$MIvio %in% c("scalarReal0.4", "scalarReal0.6")),]
dodge_w <- 1

pMeanDifferences <- ggplot(stats, 
      aes(x = factor(corGO), y = M,
          group = interaction(factor(varS), factor(R2)), shape = factor(varS), 
          color = factor(R2))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  labs(y = expression(Delta*R^2~"(items - score)"), x = expression(rho[GY])) + 
  facet_grid(miItemsPlot ~ MIvio, )

(pMeanDifferences <- themeFunctionPaper(pMeanDifferences))

# calculate descriptive statistics
funList <- list(mean = mean, median = median, sd = sd, min = min, max = max, 
                p05 = function(data) quantile(data, probs = .05),
                q25 = function(data) quantile(data, probs = .25),
                q75 = function(data) quantile(data, probs = .75), 
                p95 = function(data) quantile(data, probs = .95))

# overall
descriptive <- lapply(names(funList), function(iFun) {
  agg <- aggregate(M ~ 1, 
                   data = stats, 
                   FUN = funList[[iFun]])
  names(agg)[names(agg) == "M"] <- iFun
  agg
})

statsVec <- unlist(descriptive)
statsVec

# varT = 0.5
descriptiveVarT <- lapply(names(funList), function(iFun) {
  agg <- aggregate(M ~ 1, 
                   data = stats[stats$varT == 0.5,], 
                   FUN = funList[[iFun]])
  names(agg)[names(agg) == "M"] <- iFun
  agg
})

(statsVecVarT <- unlist(descriptiveVarT))


ggplot(stats, aes(x = M)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    color = "white",
    fill = "grey70"
  ) +
  geom_density(linewidth = 1) +
  labs(
    x = expression(Delta*R^2~"condition-level mean"),
    y = "Density"
  ) +
  theme_classic()


dim(stats[stats$M > 0.01, ])[1]/dim(stats)[1]
dim(stats[stats$M > 0.02, ])[1]/dim(stats)[1]
dim(stats[stats$M > 0.05, ])[1]/dim(stats)[1]
dim(stats[stats$M > 0.1, ])[1]/dim(stats)[1]
dim(stats[stats$M > 0.15, ])[1]/dim(stats)[1]

dim(stats[stats$MIvio != "specific" & stats$M > 0.01, ])[1]/dim(stats)[1]
dim(stats[stats$MIvio != "specific" & stats$M > 0.02, ])[1]/dim(stats)[1]
dim(stats[stats$MIvio != "specific" & stats$M > 0.05, ])[1]/dim(stats)[1]
dim(stats[stats$MIvio != "specific" & stats$M > 0.1, ])[1]/dim(stats)[1]
dim(stats[stats$MIvio != "specific" & stats$M > 0.15, ])[1]/dim(stats)[1]

range(stats[!(stats$MIvio %in% c("specific", "scalar")), "M"])
plot(stats[!(stats$MIvio %in% c("specific", "scalar")),]$M)
stats[!(stats$MIvio %in% c("specific", "scalar")) & stats$M > 0.05, ]
dim(stats[!(stats$MIvio %in% c("specific", "scalar")) & stats$M > 0.05, ])[1]/dim(stats)[1]

ggplot(stats[stats$MIvio != "specific",], aes(x = M)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    color = "white",
    fill = "grey70"
  ) +
  geom_density(linewidth = 1) +
  labs(
    x = expression(Delta*R^2~"condition-level mean"),
    y = "Density"
  ) +
  theme_classic()

plot(stats$M)
plot(stats[stats$MIvio != "specific",]$M)

################################################################################
# difference plots (paper selection)
################################################################################
# performance measures (holdout test set) to plot
measureVec <- c("Rsquared_HOtest", "R2pred_HOtest", "RMSE_HOtest")
stopifnot(all(measureVec %in% unique(performMeasures$measure_type)))

diffSub <- performMeasures[performMeasures$measure_type %in% measureVec &
                             performMeasures$model %in% c("diff_ENET", "diff_group_ENET"),]

rm(performMeasures)
gc()

diffSub$model <- sub("^.{5}", "", diffSub$model)
diffSub <- tidyr::separate(diffSub, model, into = c("predSet", "predMeth"), 
                           sep = "_", fill = "left")
diffSub$predSet <- ifelse(is.na(diffSub$predSet), "items", diffSub$predSet)

diffSub$R2 <- factor(diffSub$R2)

# # nuance factor conditions
# if (any(grepl("NF", diffSub$MIvio))){
#   diffSub$corNF <- ifelse(grepl("NF", diffSub$MIvio),
#                           substr(diffSub$MIvio, 3, nchar(as.character(diffSub$MIvio[1]))),
#                           0)
#   diffSub$MIvio <- ifelse(grepl("NF", diffSub$MIvio), "NF", diffSub$MIvio)
# }
# 
# diffSub <- diffSub[!(diffSub$MIvio %in% c("scalarReal0.4", "scalarReal0.6")), ]
diffSub$MIvio <- factor(diffSub$MIvio, 
                        levels = c("metric", "scalar", "specific"), 
                        labels = c("Metric", "Scalar", "Nuance"))



diffSub$corGO_lab <- paste0("r[GY]==", diffSub$corGO)
diffSub$R2_lab    <- paste0("R[F]^2==", diffSub$R2)
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

# analytic population Delta R² (99_checkProof.R): one non-invariant item, OLS, N -> infinity, without G
load(paste0(resFolder, "/", depFolder, "/deltaR2_analytic.rda"))
deltaR2Analytic$MIvio <- factor(deltaR2Analytic$MIvio, levels = c("metric", "scalar", "specific"),
                                labels = c("Metric", "Scalar", "Nuance"))
deltaR2Analytic$x_nested <- factor(paste(deltaR2Analytic$miItemsPlot, deltaR2Analytic$MIvio, sep = "|"),
                                   levels = levels(diffSub$x_nested))
deltaR2Analytic$corGO_lab <- paste0("r[GY]==", deltaR2Analytic$corGO)
deltaR2Analytic$R2_lab    <- paste0("R[F]^2==", deltaR2Analytic$R2)

# generate vectors to iterate over
varTvec <- unique(diffSub$varT)
predSetVec <- unique(diffSub$predSet)
nVec <- unique(diffSub$N)
cAlphaVec <- unique(diffSub$cAlpha)
Pvec <- unique(diffSub$P)

# # test it: items_1000_0.25_0.8
# iVarT <- 0.5
# iPredSet <- "items" # "group"
# iN <- 300
# iAlpha <- 0.8
# iP = 4

shapeVec <- c(3, 4)
colVecViridis <- c("#addd8e", "#1f9e89", "#1d4f91")
# colVecViridis <- c("#337538", "#0072B2", "#CC79A7")
# colVec <- c( "#006400", "#FF7F50", "#B22222") 

aggregate(diffSub, values ~ measure_type + varT + N + predSet, range)

##### robustness checks #####
# condition means across replications, one row per condition and predictor set
condMeans <- diffSub %>%
  filter(measure_type %in% c("Rsquared_HOtest", "R2pred_HOtest", "RMSE_HOtest"), !is.na(MIvio)) %>%
  group_by(measure_type, predSet, MIvio, miItemsPlot, P, cAlpha, varT, corGO, R2, N) %>%
  summarize(M = mean(values, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = measure_type, values_from = M)

# conditions sorted by absolute difference between Delta R² and Delta R²_pred (condition means)
condDiff <- condMeans %>%
  mutate(diff = Rsquared_HOtest - R2pred_HOtest,
         absDiff = abs(diff)) %>%
  arrange(desc(absDiff))

print(head(condDiff, 20), width = Inf)

# save condition means and differences between Delta R² and Delta R²_pred
save(condDiff, file = paste0(resFolder, "/", depFolder, "/condDiff_R2_R2pred.rda"))

# calculate correlations between different dependent measures to report as robustness checks
(corR2R2pred <- bind_rows(mutate(condMeans, predSet = "all"), condMeans) %>%
    group_by(predSet) %>%
    summarize(nCond = n(),
              r_R2 = cor(Rsquared_HOtest, R2pred_HOtest),
              rho_R2 = cor(Rsquared_HOtest, R2pred_HOtest, method = "spearman"),
              r_RMSE = cor(Rsquared_HOtest, -RMSE_HOtest),
              rho_RMSE = cor(Rsquared_HOtest, -RMSE_HOtest, method = "spearman"),
              maxAbsDiff_R2 = max(abs(Rsquared_HOtest - R2pred_HOtest)),
              .groups = "drop"))

# y-axis labels and plot file names per measure
yLabs <- list(Rsquared_HOtest = expression(Delta*R^2~"(items - score)"),
              R2pred_HOtest = expression(Delta*R[pred]^2~"(items - score)"),
              RMSE_HOtest = expression(Delta*RMSE~"(items - score)"))
plotLabs <- c(Rsquared_HOtest = "R2", R2pred_HOtest = "R2pred", RMSE_HOtest = "RMSE")

# identical scaling for all result plots (to a certain degree)
ylims_items    <- c(-0.17, 0.28)
ylims_varT05    <- c(-0.06, 0.11)
ylims_itemsLargeN    <- c(-0.06, 0.28)
ylims_group <- c(-0.35, 0.28)

for (iMeasure in measureVec) {
  for (iVarT in varTvec) {
    for (iPredSet in predSetVec) {
      for (iN in nVec) {  
        for (iAlpha in cAlphaVec) {
          for (iP in Pvec) {
            pTmp <- ggplot(diffSub[diffSub$measure_type == iMeasure &
                                     diffSub$varT == iVarT &
                                     diffSub$N == iN &
                                     diffSub$P == iP &
                                     diffSub$predSet == iPredSet &
                                     diffSub$cAlpha == iAlpha,], 
                           aes(x = x_nested, y = values,
                               # group = interaction(P, MIvio, miItemsPlot), 
                               group = interaction(MIvio, miItemsPlot), 
                               color = miItemsPlot, fill = miItemsPlot)) +
              geom_vline(xintercept = seq(3.5, 7.5, by = 3), linetype = "dotted",
                         colour = "grey70", linewidth = 0.8) +
              geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
              geom_violin(
                #aes(linetype = P), 
                position = position_dodge(width = dodge_w),
                width = 1.25,
                trim = TRUE,
                alpha = 0.90,
                linewidth = 0.65) +
              stat_summary(#aes(shape = P),
                fun = mean,
                shape = 4,
                geom = "point",
                color = "black",
                position = position_dodge(width = dodge_w),
                size = 1.7,
                stroke = 0.5) +
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
                    name  = c("Metric", "Scalar", "Nuance",
                              "Metric", "Scalar", "Nuance",
                              "Metric", "Scalar", "Nuance")
                  )
                ), drop = FALSE
              ) +
              guides(linetype = "none", color = "none", fill = "none") +
              labs(y = yLabs[[iMeasure]], x = "source of item-specific variance") +
              facet_grid(corGO_lab ~ R2_lab, labeller = label_parsed)
            
            # add analytic population Delta R² (only defined for models without G)
            if (iMeasure == "Rsquared_HOtest" & iPredSet == "items") {
              pTmp <- pTmp +
                # geom_errorbar(data = deltaR2Analytic[deltaR2Analytic$varT == iVarT &
                #                                        deltaR2Analytic$cAlpha == iAlpha &
                #                                        deltaR2Analytic$P == iP, ],
                #               aes(x = x_nested, ymin = deltaR2, ymax = deltaR2), inherit.aes = FALSE,
                #               width = 0.7, linewidth = 0.9, color = "#d7301f")
                geom_point(data = deltaR2Analytic[deltaR2Analytic$varT == iVarT &
                                                    deltaR2Analytic$cAlpha == iAlpha &
                                                    deltaR2Analytic$P == iP, ],
                           aes(x = x_nested, y = deltaR2), inherit.aes = FALSE,
                           shape = 23, size = 3.5, stroke = 0.9, fill = alpha("white", 0.5), color = "#d7301f")
              nLayers <- length(pTmp$layers)
              pTmp$layers <- pTmp$layers[c(1:3, nLayers, 4:(nLayers - 1))] # analytic line before mean cross
            }
            
            pTmp
            if (iMeasure == "Rsquared_HOtest") {
              if (iPredSet == "items") {
                if (iN > 100) {
                  if (iVarT == 0.25) {
                    ylims <- ylims_itemsLargeN
                  } else {
                    ylims <- ylims_varT05
                  }
                } else {
                  ylims <- ylims_items
                }
              } else {
                ylims <- ylims_group
              }
              yBreaks <- seq(floor(ylims[1]*5)/5, ceiling(ylims[2]*5)/5, by = 0.1)
            } else {
              # identical scaling across alpha and P (same measure, predictor set, N and varT)
              ylims <- range(diffSub$values[diffSub$measure_type == iMeasure &
                                              diffSub$predSet == iPredSet &
                                              diffSub$N == iN &
                                              diffSub$varT == iVarT &
                                              !is.na(diffSub$MIvio)], na.rm = TRUE)
              yBreaks <- pretty(ylims)
            }
            
            
            pTmp <- pTmp + coord_cartesian(ylim = ylims) +
              scale_y_continuous(
                breaks = yBreaks
              )
            
            (pTmp <- themeFunctionPaper(pTmp))
            
            plotName <- paste0(plotFolder, "/", plotPerformanceFolder, 
                               "/pDiff", plotLabs[[iMeasure]],"_Ho_", iPredSet, "_N", iN, "_varT", iVarT, 
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
}
