setwd("reanalyses/McClureAmmermanJacobucci2024")
source('setParameters.R')

require(ggplot2)

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

################################################################################
# analyse results
################################################################################
# load results data
load("simResults.rda")

out <- data.frame(do.call(rbind, out))
colnames(out)
# names(out) <-  c("Rep", "N","J","betaT","betaS", "measVar",
#                  "MI", "nGroups", "miItems", "miWeight", "mSpecific",
#                  "sigmaEy", "formVSref", "simR2",
#                  "item_rmse","comp_rmse", "itemCov_rmse","compCov_rmse",  "true_rmse",
#                  "item_R2", "comp_R2", "itemCov_R2", "compCov_R2", "true_R2")

range(out$simR2)

out <- tidyr::pivot_longer(out, 
                           !c(Rep, N, J, betaS, betaT, measVar, MI, nGroups, miItems, 
                              miWeight, mSpecific, sigmaEy, formVSref, simR2),
                           names_to = "measures", values_to = "values")

str(out)

# produce a variable that describes conditions in words instead of numbers 
# ! always think about betaS, too! this is the difference between the two different
#   questions:  "false positives" - can we provoke superiority of the items by violating MI assumptions?
#               "false negatives" - does MI violation cover superiority of the items?
out$condition <- NA
unique(out$betaS)
out$condition[which(out$betaS == 0 &
                      out$nGroups == 1)] <- "noItemEffect"
out$condition[which(out$betaS != 0 & 
                      is.na(out$miWeight))] <- "largeItemEffect"

out$condition[which(out$miWeight == 0)] <- "configuralMIvio"
out$condition[which(out$miWeight == -0.8)] <- "metricMIvio"
out$condition[which(out$miWeight == 1)] <- "scalarMIvio"
out$condition <- factor(out$condition, 
                           levels = c("noItemEffect", "largeItemEffect", 
                                      "configuralMIvio", "metricMIvio", "scalarMIvio"))

# # check that I really do not need miWeight as a variable that produces issues in aggregate
# unique(out$miWeight)
# unique(out$condition)
# miWeightIDX <- which(is.na(out$miWeight))
# nGroupsIDX <- which(out$nGroups == 1)
# all(miWeightIDX == nGroupsIDX) # information from miWeight in nGroups, too

# separate model (= which predictors are in?) and measure (dependent measure)
out <- tidyr::separate(out, measures, c("model", "measures"), sep = "_", remove = T)
out$measures <- factor(out$measures, levels = c("R2", "rmse"))

# dependent measures in wide format to calculate M and SD separately for DVs
out <- tidyr::pivot_wider(out, names_from = measures, values_from = values) 

# calculate mean values and standard deviation of the values
# remove miWeight as NAs in miWeight are dropped and every condition is represented in other variables, too 
plotRes <- aggregate(cbind(R2, rmse) ~ N + J + betaT + betaS + measVar + MI + nGroups +
                       miItems + mSpecific + sigmaEy + formVSref + model + condition, 
                     data = out, function(x) {cbind(mean(x), sd(x))})

plotRes <- do.call(data.frame, plotRes)
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.1", "_M")
colnames(plotRes) <- stringr::str_replace_all(colnames(plotRes), "\\.2", "_SE")

# change variables into factors for plotting categorial variables
cols2fac <- c("N", "J", "betaT", "betaS", "measVar", "MI", "nGroups",
              "miItems", "mSpecific", "sigmaEy", "formVSref", "model", "condition")
plotRes[cols2fac] <- lapply(plotRes[cols2fac], factor)

# separate dependent measure (R2 & RMSE) from descriptive statistic (M & SD)
plotRes <- tidyr::pivot_longer(plotRes, 
                               cols = c("R2_M", "R2_SE", "rmse_M", "rmse_SE"),
                               names_to = c("levels", "DV"), 
                               names_sep = "_",
                               values_to = "values")

plotRes <- tidyr::pivot_wider(plotRes,
                              names_from = DV, 
                              values_from = values)

################################################################################
# plot results
################################################################################
# these are results graphics for different parameter manipulation grids to show
#   different effects on the superiority of items over composites

################################################################################
# # 1. manipulate one item; different measurement invariance violations
################################################################################
# plot R2
unique(plotRes$betaS)
# betaS.cond = 1 == betaS = 0
# betaS.cond = 4 == betaS = 2
# mean difference with mSpecific = 5 increases variance of the specific item effect 
#   and hence the variance of the variable
subset_higherVar <- plotRes[which(plotRes$measVar == 1 & plotRes$betaT == 0.5 &
                ((plotRes$MI == 1 & (plotRes$betaS == 0 | plotRes$betaS == 2)) |
                   (plotRes$MI == 0 & plotRes$mSpecific == 5 & plotRes$betaS == 0 & plotRes$miItems == 1))), ]

(plotMImani_higherVar <- ggplot(subset_higherVar[subset_higherVar$levels == "R2" &
                               subset_higherVar$model != "true", ],
                      aes(y = M, x = N,
                          group = interaction(J, model),
                          shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_higherVar.eps"),
#                 plot = plotMImani_higherVar,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_higherVar.png"),
#                 plot = plotMImani_higherVar,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# only condition with constant variance in specific item effect
subset_equalVar <- plotRes[which(plotRes$measVar == 1 & plotRes$betaT == 0.5 &
                          ((plotRes$MI == 1 & (plotRes$betaS == 0 | plotRes$betaS == 2)) |
                             (plotRes$MI == 0 & plotRes$mSpecific == 1 & plotRes$betaS == 0 & plotRes$miItems == 1))), ]

(plotMImani_equalVar <- ggplot(subset_equalVar[subset_equalVar$levels == "R2" &
                                                 subset_equalVar$model != "true", ],
                                aes(y = M, x = N,
                                    group = interaction(J, model),
                                    shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ condition, scales = "free_y", labeller = label_both, switch = "y"))


# ggplot2::ggsave(filename = paste0(plotFolder, "/01_plotMImanipulations_equalVar.png"),
#                 plot = plotMImani_equalVar,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# wie passt Abwesenheit von Unterschieden in prediction auf Basis von items und
#   composite mit Pratiwi Simulation zusammen?

# skalare Messinvarianzverletzung durch Mittelwertsunterschiede in specific item
#   effects zwischen Gruppen; aber specific item effects ohne Beziehung zu Outcome
# Wie wird durch Verletzung der skalaren MI item prediction > composite prediction?
# Wie ist das übertragbar/vereinbar mit der Vorstellung eines inhaltlichen Nuancen-Faktors?
#   wenn Mittelwertsunterschiede zwischen Items ohne Effekt auf Outcome prediction
#     vergrößern, brauchen wir dann items für Nuancen?
#   ist Nuance nichts anderes als die Variable, die die Mittelwertsunterschiede erlären kann?
#   Nuance als specific item effect als latente Variable?

################################################################################
# 2. add conditions with different deviations in mean value of the specific item effect
################################################################################
subset <- plotRes[which(plotRes$measVar == 1 & plotRes$betaT == 0.5 &  
                          ((plotRes$MI == 1 & (plotRes$betaS == 0 | plotRes$betaS == 2)) |
                             (plotRes$MI == 0 & plotRes$betaS == 0 & plotRes$miItems == 1))),] 

(plotMIspecific <- ggplot(subset[subset$levels == "R2" & subset$model != "true", ], 
                          aes(y = M, x = N, 
                              group = interaction(J, model, mSpecific),
                              shape = J, linetype = model, colour = mSpecific)) +
    geom_point() +
    geom_line() +
    #scale_linetype_manual(values = c("dashed", "dotted", "solid")) + 
    scale_color_manual(values = c("#0033FF", "#0033FF", "#006600", "#FF6600", "#CC0000")) +
    facet_grid(formVSref + sigmaEy ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_overview.eps"),
#                 plot = plotMIspecific,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotMImeanSpecificEffect_overview.png"),
#                 plot = plotMIspecific,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


# the additional manipulation of different means in scalar MI violation only affects the
#   scalar MI  violation condition -> only plot this subset of the data
# mSpecific = 1 means that the variance of the specific item effect is fixed to 0.25
onlyScalarMI <- subset[subset$levels == "R2" &
                         subset$model != "true" &
                         subset$condition == "scalarMIvio",]

(plotMIscalarSpecific <- ggplot(onlyScalarMI,
                                aes(y = M, x = N,
                                    group = interaction(J, model, mSpecific),
                                    shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    # scale_color_manual(values = c("#CC0000", "#FF3300", "#006600")) +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ mSpecific, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean.eps"),
#                 plot = plotMIscalarSpecific,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/02_plotScalarMImean.png"),
#                 plot = plotMIscalarSpecific,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# with any kind of measurement invariance violation:
#   for small sample sizes variance in item based prediction is higher (for small sample sizes)
#   compared to composite based prediction; 
#   since item based prediction is worse for item based prediction supperior predictive 
#   performance of items cannot be effect of selective reporting
# the stronger mean deviations for specific item effects, the higher the bias in the composite
#   more bias for large number of items in the composite although number of items deviating from 
#   the mean of 0 is constant (only 1 item!) (this is counter intuitive!!!)
# the stronger mean deviations for specific item effects, the higher advantages in 
#   the predictive performance (R2)
#     - for specific effect 0 vs. 2 there is seems to be no advantage in item based 
#       prediction
#     - for specific effect 0 vs. {5, 10} there are advantages in predictive performance 
#       for items compared with the composite

################################################################################
# 3. add conditions with varying number of items that are subject to measurement invariance violations
################################################################################
subset <- plotRes[which(plotRes$measVar == 1 & plotRes$betaT == 0.5 &
                          # MI is violated, no specific item effect, constant variance in specific item effect
                          (plotRes$MI == 0 & plotRes$betaS == 0 & plotRes$mSpecific == 1)), ]

unique(subset$miItems)
# # the actual number of items in miItems.cond = {2, 3} differs depending on the number of
# #   items in the composite {J/3, floor(J/2)}
# #   -> backtracking actual item numbers (depending on J) to item number conditions (independent of J)
# subset$miItems.cond <- plyr::mapvalues(subset$miItems,
#                                        from = c(0, 1, 3, 7, 4, 10),
#                                        to = c(0, 1, 2, 2, 3, 3))

str(subset)
(plotMInrItems_overview <- ggplot(subset[subset$levels == "R2" &
                                  subset$model != "true", ],
                         aes(y = M, x = N,
                             group = interaction(J, model, miItems),
                             shape = J, linetype = model, colour = miItems)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_overview.eps"),
#                 plot = plotMInrItems_overview,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_overview.png"),
#                 plot = plotMInrItems_overview,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


onlyConfiguralMI <- subset[subset$levels == "R2" & 
                         subset$model != "true" &
                         subset$condition == "configuralMIvio",]

(plotMI_configuralNrItems <- ggplot(onlyConfiguralMI,
                               aes(y = M, x = N,
                                   group = interaction(J, model, miItems),
                                   shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ miItems, scales = "free_y", labeller = label_both, switch = "y"))

onlyMetricMI <- subset[subset$levels == "R2" & 
                             subset$model != "true" &
                             subset$condition == "metricMIvio",]

(plotMI_metricNrItems <- ggplot(onlyMetricMI,
                                aes(y = M, x = N,
                                    group = interaction(J, model, miItems),
                                    shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ miItems, scales = "free_y", labeller = label_both, switch = "y"))


onlyScalarMI <- subset[subset$levels == "R2" & 
                         subset$model != "true" &
                         subset$condition == "scalarMIvio",]

(plotMI_scalarNrItems <- ggplot(onlyScalarMI,
                               aes(y = M, x = N,
                                   group = interaction(J, model, miItems),
                                   shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy ~ miItems, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_configural.png"),
#                 plot = plotMI_configuralNrItems,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_metric.png"),
#                 plot = plotMI_metricNrItems,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# ggplot2::ggsave(filename = paste0(plotFolder, "/03_plotMInrItems_scalar.png"),
#                 plot = plotMI_scalarNrItems,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")


# the number of manipulated items does not seem to make any difference for configural 
#   and metric measurement invariance violations, i.e., the predictive performance of
#   the item based models do not surpass the composite based condition (unit weighted!)
# only for scalar measurement invariance violations the number of items with deviating
#   means of the specific effect on the items does make a difference in terms of 
#   predictive performance (R2, RMSE and bias)

################################################################################
# 4. add conditions with specific effect on outcome
################################################################################
subset <- plotRes[plotRes$measVar == 1 & plotRes$betaT == 0.5 &
                    ((plotRes$MI == 1 & (plotRes$betaS == 0 | plotRes$betaS == 2)) |
                       (plotRes$MI == 0 & plotRes$mSpecific == 1 &
                          (plotRes$betaS == 0 | plotRes$betaS == 2))), ]

# unique(subset$miItems) # number of MI violating items is constant across all conditions
# unique(subset$mSpecific) # mean value of scalar MI violating items only has one value
# unique(subset$betaS)
# unique(subset$model)

subset <- subset[subset$model != "true",]
subset$betaS <- factor(subset$betaS)
subset <- tidyr::unite(subset, "model_betaS", c(model, betaS), sep = " x ", remove = F)
unique(subset$model_betaS)
subset$model_betaS <- factor(subset$model_betaS, 
                             levels = c("comp x 0", "compCov x 0", "comp x 2", "compCov x 2",
                                        "item x 0", "itemCov x 0", "item x 2", "itemCov x 2"))
# ggplot(subset[subset$levels == "R2", ],
#        aes(y = M, x = N,
#            group = interaction(J, model, model_betaS),
#            shape = J, linetype = model, colour = model_betaS)) +
#   geom_point() +
#   geom_line() +
#   #scale_color_manual(values = c("#006600", "#33CC33", "#CC0000", "#FF3300")) +
#   facet_grid(formVSref + sigmaEy ~ condition, scales = "free_y", labeller = label_both, switch = "y")

(plotMIxBetaS <- ggplot(subset[subset$levels == "R2" & 
                                 (subset$miItems == 0 | subset$miItems == 1),],
                        aes(y = M, x = N,
                            group = interaction(J, model),
                            shape = J, linetype = model, colour = model)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("#CC0000", "#FF3300", "#006600", "#33CCCC")) +
    facet_grid(formVSref + sigmaEy + betaS ~ condition, scales = "free_y", labeller = label_both, switch = "y"))

# # save plots as files
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS_overview.eps"),
#                 plot = plotMIxBetaS,
#                 device = cairo_ps,
#                 dpi = 300,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")
# 
# ggplot2::ggsave(filename = paste0(plotFolder, "/04_plotMIxBetaS_overview.png"),
#                 plot = plotMIxBetaS,
#                 width = 13.08,
#                 height = 12.18,
#                 units = "in")

# zusätzlich einen spezifischen Effekt dazu zu simulieren katapultiert bei scalar
#   MI violations die predictive performance unabhängig von Anzahl der Items weit über 
#   die predictive performance der composite
# bei configural und metric measurement invariance ist bei 21 items die predictive 
#   performance trotzdem nicht viel besser

# combined condition (items + factor) macht keinen Unterschied (unit weighting) 