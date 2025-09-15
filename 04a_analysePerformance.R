# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

# calculate ANOVA for difference in R² to check for effects of manipulated variables in simulation 
library(afex) # for aov_ez()
library(effectsize) # for efect size calculation; generalized eta²
library(emmeans)
library(tidyverse)

resFolder <- "results"
depFolder <- "dependentMeasures"

plotFolder <- "plots"
createFolder(plotFolder)
plotAnovaFolder <- "ANOVA"
createFolder(paste0(plotFolder, "/", plotAnovaFolder))

resFile <- paste0(resFolder, "/", depFolder, "/performPerSample_differences.rda")
load(resFile)

################################################################################
# ANOVA on difference measures
################################################################################
# get subset to calculate ANOVA on
unique(performMeasures$measure_type)
unique(performMeasures$model)

diffSub <- performMeasures[performMeasures$model %in% c("diff_LM", "diff_Ridge", "diff_ENET") &
                  performMeasures$measure_type == "Rsquared_HOtest",]

rm(performMeasures)
gc()

# change variables to factors
col2fac <- c("MIvio", "P", "Q", "cAlpha", "varT", "corGO", "R2", "sample", "N", 
             "measure_type", "miItemsPlot", "model")
diffSub[col2fac] <- lapply(diffSub[col2fac], factor)

################################################################################
# mixed ANOVA for difference in R² between items and score
################################################################################
# mixed ANOVA with ...
# ... id = sample but ID 
#     (independent samples between simulated conditions expect from N)
# ... dv = {diff_LM, diff_Ridge, diff_ENET} in separate mixed ANOVAs
# ... between: 3 x 3 x 2 x 2 x 2 x 3 x 2
#   MIvio (3) {specific, metric, scalar}
#   miItemsPlot (3) {0, 1, P/2}
#   P (2) {4, 8}
#   cAlpha (2) {0.6, 0.8}
#   varT (2) {0.25, 0.5}
#   corGO (3) {0, 0.2, 0.5}
#   R2 (2) {0.3, 0.8}
# ... within:
#   N (3) {100, 300, 1000}

# Type 3 sums of squaress (e.g., Maxwell and Delaney, 2004)
# contr.sum for effects-coding for the categorical variables
# outcome: generalized eta2 (Olejnik and Algina 2003)

diffSub <- tidyr::pivot_wider(diffSub, names_from = model, values_from = values) 

# ID variables according to within & between structure assumed in the data
diffSub <- diffSub %>%
  group_by(N) %>%
  mutate(ID = seq_len(n())) %>%
  ungroup()

outVec <- c("diff_LM", "diff_Ridge", "diff_ENET")
for (iOut in outVec) {
  
  # ANOVA:
  anovaTmp <- aov_ez(id = "ID",
                         dv = iOut,
                         data = diffSub,
                         between = c("MIvio", "P", "cAlpha", "varT", "corGO", "R2",
                                     "miItemsPlot"),
                         within = "N")
  
  eta2Tmp <- eta_squared(
    anovaTmp, # fitted model
    partial = FALSE, # not partial!
    generalized = TRUE, # generalized eta squared
    ci = 0.95,
    verbose = TRUE)
  
  # sort generalized eta-squared results
  # which higher order interactions do we need to illustrate to report simulation results?
  (eta2Tmp.ordered <- eta2Tmp[order(eta2Tmp$Eta2_generalized, decreasing = T),])
  
  eta2Thresh <- 0.05
  plotEta2Tmp <- eta2Tmp.ordered[eta2Tmp.ordered$Eta2_generalized > eta2Thresh,]
  
  sortIdx <- order(plotEta2Tmp$Eta2_generalized, decreasing = T)
  plotEta2Tmp$Parameter <- factor(plotEta2Tmp$Parameter,
                                      levels = plotEta2Tmp$Parameter[sortIdx])
  
  # plot generalized eta^2 for mixed ANOVA with test R^2 as dependent variable and
  # ... model as within factor
  # ... N, pTrash, rel, R2 and lin_inter as between factors
  # comparison in test R2 between models with and without option to extract interaction effects is not fair
  # -> next: plot generalized eta^2 for only ENETw and GBM which both can take interaction effects into account
  (pEta2Tmp <- ggplot(plotEta2Tmp, aes(x = Parameter, y = Eta2_generalized)) +
      geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
      geom_text(aes(label=round(Eta2_generalized, 2)),
                #angle = 90, hjust = 1.5, vjust=0.5,
                angle = 0, vjust=0.5, hjust = 1.5,
                color="white", size=4)+
      ylab("generalisiertes eta^2") +
      coord_flip()) 
  pEta2Tmp <- themeFunctionPaper(pEta2Tmp)
  
  ggplot2::ggsave(filename = paste0(plotFolder, "/", plotAnovaFolder, "/mixedANOVA_", 
                                    iOut, "_thresh", eta2Thresh, ".png"),
                  plot = pEta2Tmp,
                  width = 17.52,
                  height = 10.76,
                  units = "in")
  
  eta2modelTable <- plotEta2Tmp[plotEta2Tmp$Eta2_generalized >= eta2Thresh,]
  
  print(xtable::xtable(eta2modelTable, type = "latex"),
        file = paste0(plotFolder, "/", plotAnovaFolder, "/mixedANOVA_", 
                      iOut, "_thresh", eta2Thresh, ".tex"))
  
}

# MIvio interactions with lots of simulation parameters! 
# -> different manipulations seem to affect different MIvio levels differently
# interpretation of the other manipulated parameters however becomes difficult
# therefore, run ANOVA for every level of the MIvio manipulation separately

################################################################################
# mixed ANOVA for every MIvio separately
################################################################################
# mixed ANOVA for every combination of model and MIvio manipulation
#   model (3) {LM, Ridge, ENET}
#   MIvio (3) {specific, metric, scalar}
# mixed ANOVA with ...
# ... id = sample but ID 
#     (independent samples between simulated conditions expect from N)
# ... dv = {diff_LM, diff_Ridge, diff_ENET} in separate mixed ANOVAs
# ... between: 3 x 2 x 2 x 2 x 3 x 2
#   miItemsPlot (3) {0, 1, P/2}
#   P (2) {4, 8}
#   cAlpha (2) {0.6, 0.8}
#   varT (2) {0.25, 0.5}
#   corGO (3) {0, 0.2, 0.5}
#   R2 (2) {0.3, 0.8}
# ... within:
#   N (3) {100, 300, 1000}

# Type 3 sums of squaress (e.g., Maxwell and Delaney, 2004)
# contr.sum for effects-coding for the categorical variables
# outcome: generalized eta2 (Olejnik and Algina 2003)

# ID variables according to within & between structure assumed in the data
diffSub2 <- diffSub %>%
  group_by(N, MIvio) %>%
  mutate(ID = seq_len(n())) %>%
  ungroup()

for (iMI in unique(diffSub2$MIvio)) {
  for (iOut in outVec) {
    
    # ANOVA:
    anovaTmp <- aov_ez(id = "ID",
                       dv = iOut,
                       data = diffSub2[diffSub2$MIvio == iMI,],
                       between = c("P", "cAlpha", "varT", "corGO", "R2",
                                   "miItemsPlot"),
                       within = "N")
    
    eta2Tmp <- eta_squared(
      anovaTmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    # sort generalized eta-squared results
    # which higher order interactions do we need to illustrate to report simulation results?
    (eta2Tmp.ordered <- eta2Tmp[order(eta2Tmp$Eta2_generalized, decreasing = T),])
    
    eta2Thresh <- 0.05
    plotEta2Tmp <- eta2Tmp.ordered[eta2Tmp.ordered$Eta2_generalized > eta2Thresh,]
    
    sortIdx <- order(plotEta2Tmp$Eta2_generalized, decreasing = T)
    plotEta2Tmp$Parameter <- factor(plotEta2Tmp$Parameter,
                                    levels = plotEta2Tmp$Parameter[sortIdx])
    
    # plot generalized eta^2 for mixed ANOVA with test R^2 as dependent variable and
    # ... model as within factor
    # ... N, pTrash, rel, R2 and lin_inter as between factors
    # comparison in test R2 between models with and without option to extract interaction effects is not fair
    # -> next: plot generalized eta^2 for only ENETw and GBM which both can take interaction effects into account
    (pEta2Tmp <- ggplot(plotEta2Tmp, aes(x = Parameter, y = Eta2_generalized)) +
        geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
        geom_text(aes(label=round(Eta2_generalized, 2)),
                  #angle = 90, hjust = 1.5, vjust=0.5,
                  angle = 0, vjust=0.5, hjust = 1.5,
                  color="white", size=4)+
        ylab("generalisiertes eta^2") +
        coord_flip()) 
    pEta2Tmp <- themeFunctionPaper(pEta2Tmp)
    
    ggplot2::ggsave(filename = paste0(plotFolder, "/", plotAnovaFolder, "/mixedANOVA_", 
                                      iOut, "_MIvio", iMI, "_thresh", eta2Thresh, ".png"),
                    plot = pEta2Tmp,
                    width = 17.52,
                    height = 10.76,
                    units = "in")
    
    eta2modelTable <- plotEta2Tmp[plotEta2Tmp$Eta2_generalized >= eta2Thresh,]
    
    print(xtable::xtable(eta2modelTable, type = "latex"),
          file = paste0(plotFolder, "/", plotAnovaFolder, "/mixedANOVA_", 
                        iOut, "_MIvio", iMI, "_thresh", eta2Thresh, ".tex"))
    
  }
}
