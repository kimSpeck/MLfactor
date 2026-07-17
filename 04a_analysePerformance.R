# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

# calculate ANOVA for difference in R² to check for effects of manipulated variables in simulation 
library(afex) # for aov_ez()
library(effectsize) # for efect size calculation; generalized eta²
library(emmeans)
library(tidyverse)
library(reshape2)
library(ggplot2)

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

diffVars <- unique(performMeasures$model)[grepl("diff", unique(performMeasures$model))]
diffSub <- performMeasures[performMeasures$model %in% diffVars &
                  performMeasures$measure_type == "Rsquared_HOtest",]

rm(performMeasures)
gc()

# change variables to factors
col2fac <- c("MIvio", "P", "Q", "cAlpha", "varT", "corGO", "R2", "sample", "N", 
             "measure_type", "miItemsPlot", "model")
diffSub[col2fac] <- lapply(diffSub[col2fac], factor)

# remove small non-invariance effect sizes (not full factorial)
diffSub <- diffSub[!(diffSub$MIvio %in% c("scalarReal0.4", "scalarReal0.6")), ]

# # get ANOVA result without the large effect size condition varT = 0.25
# unique(diffSub$varT)
# diffSub <- diffSub[diffSub$varT == 0.5, ]

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

subNoNF <- diffSub[!grepl("^NF0\\.", diffSub$MIvio), ]

# outVec <- colnames(subNoNF)[grepl("diff", colnames(subNoNF))]
outVec <- c("diff_LM", "diff_Ridge", "diff_ENET")
for (iOut in outVec) {
  
  # ANOVA:
  anovaTmp <- aov_ez(id = "ID",
                         dv = iOut,
                         # data = diffSub[diffSub$MIvio != "nuanceFactor", ],
                         #data = diffSub[!grepl("^NF0\\.", diffSub$MIvio), ],
                         data = subNoNF,
                         between = c("MIvio", "P", "cAlpha", "varT", "corGO", "R2",
                                     "miItemsPlot"),
                     # between = c("MIvio", "P", "cAlpha", "corGO", "R2",
                     #             "miItemsPlot"),
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

idxNF <- which(grepl("^NF", diffSub2$MIvio))
diffSubNF <- diffSub2[idxNF,]
diffSubNF$corNF <- as.factor(substr(diffSubNF$MIvio, 3, nchar(as.character(diffSubNF$MIvio[1]))))

diffSubNF <- diffSubNF %>%
  group_by(N) %>%
  mutate(ID = seq_len(n())) %>%
  ungroup()


MIvec <- c("metric", "scalar", "specific")    
# outVec <- colnames(diffSub2)[grepl("diff", colnames(diffSub2))]
outVec <- c("diff_LM", "diff_Ridge", "diff_ENET")

for (iMI in MIvec) {
  for (iOut in outVec) {
    
    if (iMI != "NF") {
      # ANOVA with miItemsPlot because multiple levels in factor
      anovaTmp <- aov_ez(id = "ID",
                         dv = iOut,
                         data = diffSub2[diffSub2$MIvio == iMI,],
                         between = c("P", "cAlpha", "varT", "corGO", "R2",
                                     "miItemsPlot"),
                         # between = c("P", "cAlpha", "corGO", "R2",
                         #             "miItemsPlot"),
                         within = "N")
      
    } else {
      
      # ANOVA without miItemsPlot because only one level in factor
      anovaTmp <- aov_ez(id = "ID",
                         dv = iOut,
                         data = diffSubNF,
                         # between = c("P", "cAlpha", "varT", "corGO", "R2", "corNF"),
                         between = c("P", "cAlpha", "corGO", "R2", "corNF"),
                         within = "N")
      
    }
    
    # investigate marginal means later
    anovaObjectName <- paste0("anova_", iMI, "_", iOut)
    assign(anovaObjectName, anovaTmp)
    
    eta2Tmp <- eta_squared(
      anovaTmp, # fitted model
      partial = FALSE, # not partial!
      generalized = TRUE, # generalized eta squared
      ci = 0.95,
      verbose = TRUE)
    
    # prepare merging of eta results
    eta2Tmp[,"MLmeth"] <- iOut
    eta2Tmp[,"MIvio"] <- iMI
    
    # sort generalized eta-squared results
    # which higher order interactions do we need to illustrate to report simulation results?
    (eta2Tmp.ordered <- eta2Tmp[order(eta2Tmp$Eta2_generalized, decreasing = T),])
    
    # produce more results with these objects within this script later (paper tables, ...)
    etaObjectName <- paste0("eta2_", iMI, "_", iOut)
    assign(etaObjectName, eta2Tmp.ordered)
    
    
    # produce overview plots for every ANOVA (combination of MLmeth x MIvio) 
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

##### eta2 table #####
# build tables of results for the paper
#   with results for different ANOVAs from the same prediction method but different
#   MIvios next to each other
# eta2res <- ls(pattern = "eta2_") # check objects that match the pattern
eta2List <- lapply(grep("eta2_", ls(), value = TRUE), get)
eta2Table <- do.call(rbind, eta2List)

# compare generalized eta^2 across different methods and predictor sets 
# calculate rank order coeffiecient (Kendall's Tau)
# 1st: data to wide format
# 2nd: get Kendall's Tau and weighted Kendalls' Tau  
eta2TableCor <- eta2Table[c("Parameter", "Eta2_generalized", "MLmeth", "MIvio")]
eta2TableCor <- unite(eta2TableCor, "MLmeth_MIvio", c(MLmeth, MIvio), remove = TRUE)
eta2TableCor <- tidyr::pivot_wider(eta2TableCor, names_from = MLmeth_MIvio, 
                                   values_from = Eta2_generalized) 

# remove columns with only NAs (no information!)
naCols <- apply(eta2TableCor, 2, function(x) {all(is.na(x))})
eta2TableCor <- eta2TableCor[, !naCols]

# all correlations without nuanceFactor
idxParamCol <- which(colnames(eta2TableCor) == "Parameter")
idxNF <- which(grepl("NF", colnames(eta2TableCor)))
idx <- union(idxParamCol, idxNF)

# test it
kendallTau <- cor(eta2TableCor[, -idx], method = "kendall", use = "pairwise.complete.obs")
corrplot::corrplot.mixed(kendallTau)

# # checK: without weights, the wdm package provides identical estimates as cor!
# wdmKendallTau_NF <- wdm::wdm(eta2TableCor[, idxNF], method = "kendall")
# corrplot::corrplot.mixed(wdmKendallTau_NF) 

# differences items - scores for method x MIvio
grid <- expand.grid(c("LM", "Ridge", "ENET"), c("specific", "metric", "scalar"))
diffCols <- paste0(rep("diff_", dim(grid)[1]), grid[,"Var1"], "_", grid[,"Var2"])
idxDiffs <- which(colnames(eta2TableCor) %in% diffCols)

diffGroupCols <- paste0(rep("diff_group_", dim(grid)[1]), grid[,"Var1"], "_", grid[,"Var2"])
idxGroupDiffs <- which(colnames(eta2TableCor) %in% diffGroupCols)

linearWeights <- seq(from = 1, to = 0, length.out = (dim(eta2TableCor)[1]+1))[1:dim(eta2TableCor)[1]]
plot(seq(dim(eta2TableCor)[1]), linearWeights)
wdmKendallTau_Diffs <- wdm::wdm(eta2TableCor[, idxDiffs], method = "kendall", 
                             weights = linearWeights)
colnames(wdmKendallTau_Diffs) <- rownames(wdmKendallTau_Diffs) <- substring(colnames(wdmKendallTau_Diffs), 6)
corrplot::corrplot(wdmKendallTau_Diffs, order = 'hclust', addrect = 3, 
                   addCoef.col = 'black', tl.srt = 45,
                   tl.cex = 0.9) # focus on your positive range) 

ord <- c("LM_metric", "ENET_metric","Ridge_metric",
         "LM_scalar","ENET_scalar","Ridge_scalar", 
         "LM_specific", "ENET_specific","Ridge_specific")
wdmKendallTau_Diffs <- wdmKendallTau_Diffs[ord, ord]
plotTau <- melt(wdmKendallTau_Diffs, varnames = c("row", "col"), value.name = "tau")
plotTau$row <- factor(plotTau$row, levels = rev(ord))
plotTau$col <- factor(plotTau$col, levels = ord)

# remove suffix for axis labels: "ENET_metric" -> "ENET"
strip_suffix <- function(x) sub("_(metric|scalar|specific)$", "", x)

# helper for header positions (tile centers)
header_y <- length(ord) + 0.9  # above top row of tiles

(pWeightedTau <- ggplot(plotTau, aes(col, row, fill = tau)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", tau)), size = 4) +
  # positive-only gradient (adjust limits to your range)
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B",
                      limits = c(0.45, 1), oob = scales::squish,
                      name = "Weighted \nKendall's \u03C4") +
  coord_equal(clip="off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    axis.text.y = element_text(),
    plot.margin = margin(20, 20, 40, 20)
  ) +
  scale_x_discrete(labels = strip_suffix, # position = "top", hjust = 0.5,
                   expand = expansion(mult = c(0, 0), add = c(0, 0))) +
  scale_y_discrete(labels = strip_suffix, expand = expansion(mult = c(0, 0), add = c(0, 1.2))) +
  # draw thick lines between the 3 blocks (after 3 and 6)
  geom_vline(xintercept = c(3.5, 6.5), linewidth = 1.0, color = "grey20") +
  geom_hline(yintercept = c(3.5, 6.5), linewidth = 1.0, color = "grey20") +
  # add block headers (top)
  annotate("text", x = 2, y = header_y, label = "Metric", fontface = "bold", size = 4) +
  annotate("text", x = 5, y = header_y, label = "Scalar", fontface = "bold", size = 4) +
  annotate("text", x = 8, y = header_y, label = "Nuance", fontface = "bold", size = 4))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotAnovaFolder, "/weightedTau_difference.png"),
                plot = pWeightedTau,
                width = 9.16,
                height = 9.16,
                units = "in")

wdmKendallTau_groupDiffs <- wdm::wdm(eta2TableCor[, idxGroupDiffs], method = "kendall", 
                                weights = linearWeights)
colnames(wdmKendallTau_groupDiffs) <- rownames(wdmKendallTau_groupDiffs) <- substring(colnames(wdmKendallTau_groupDiffs), 6)
corrplot::corrplot(wdmKendallTau_Diffs, order = 'hclust', addrect = 3, 
                   addCoef.col = 'black')

# 
eta2Table <- tidyr::pivot_wider(eta2Table[,c("Parameter", "Eta2_generalized", "MLmeth", "MIvio")], 
                                names_from = c(MLmeth, MIvio),  
                                values_from = Eta2_generalized,
                                names_sep = "_")

# do not remove variables below a certain threshold
eta2Thresh <- 0.05

# generate sorted tables for different subsets of ANOVA results
for (iMLmeth in unique(outVec)){
  idx <- stringr::str_detect(colnames(eta2Table), iMLmeth)
  idx[1] <- TRUE  
  
  TMPeta2 <- eta2Table[,idx]
  
  # sort
  TMPeta2$thresh <- apply(TMPeta2[,!(colnames(TMPeta2) %in% "Parameter")], 1, 
                            function(r) any(r >= eta2Thresh))
  TMPeta2$sumEta2 <- apply(TMPeta2[,!(colnames(TMPeta2) %in% c("Parameter", "thresh"))], 1, sum)
  TMPeta2$M <- apply(TMPeta2[,!(colnames(TMPeta2) %in% c("Parameter", "thresh", "sumEta2"))], 1, mean)
  TMPeta2 <- dplyr::arrange(TMPeta2, desc(M))
  
  colnames(TMPeta2)
  texTable <- TMPeta2[TMPeta2$thresh == TRUE,!(colnames(TMPeta2) %in% c("thresh", "sumEta2", "M"))]
  orderCols <- c("Parameter", paste0(iMLmeth, "_specific"), 
                 paste0(iMLmeth, "_scalar"), paste0(iMLmeth, "_metric"))
  texTable <- texTable[, orderCols]
  
  # formating: round to two
  idxNum <- stringr::str_detect(colnames(texTable), iMLmeth)
  texTable[idxNum] <- lapply(texTable[idxNum], function(x) sprintf("%.2f", x))
  
  # remove leading zero for numbers in (-1,1)
  drop0 <- function(s) {
    gsub("^(-?)0\\.", "\\1.", s) 
  }
  texTable[idxNum] <- lapply(texTable[idxNum], drop0)
  
  print(xtable::xtable(texTable, type = "latex"),
        file = paste0(plotFolder, "/ANOVA/betweenANOVA_", iMLmeth,"_", eta2Thresh, ".tex"), 
        include.rownames = FALSE, sanitize.text.function = identity)
}

# table with LM and ENET next to each other
idx <- !(stringr::str_detect(colnames(eta2Table), "Ridge"))
TMPeta2 <- eta2Table[,idx]
TMPeta2 <-TMPeta2[colnames(TMPeta2)[grepl("Parameter", colnames(TMPeta2)) | 
                              grepl("diff_ENET", colnames(TMPeta2)) | 
                              grepl("diff_LM", colnames(TMPeta2))]]

# sort
TMPeta2$thresh <- apply(TMPeta2[,!(colnames(TMPeta2) %in% "Parameter")], 1, 
                        function(r) any(r >= eta2Thresh))
TMPeta2$sumEta2 <- apply(TMPeta2[,!(colnames(TMPeta2) %in% c("Parameter", "thresh"))], 1, sum, na.rm = T)
TMPeta2$M <- apply(TMPeta2[,!(colnames(TMPeta2) %in% c("Parameter", "thresh", "sumEta2"))], 1, mean, na.rm = T)
TMPeta2 <- dplyr::arrange(TMPeta2, desc(M))

colnames(TMPeta2)
texTable <- TMPeta2[TMPeta2$thresh == TRUE,!(colnames(TMPeta2) %in% c("thresh", "sumEta2", "M"))]
# texTable[colnames(texTable)[grepl("Parameter", colnames(texTable)) | 
#                      grepl("diff_ENET", colnames(texTable)) | 
#                      grepl("diff_LM", colnames(texTable))]]
orderCols <- c("Parameter", 
               paste0("diff_ENET", "_metric"), paste0("diff_ENET", "_scalar"),
               paste0("diff_ENET", "_specific"), 
               paste0("diff_LM", "_metric"), paste0("diff_LM", "_scalar"), 
               paste0("diff_LM", "_specific"))
texTable <- texTable[, orderCols]

# formating: round to two
idxNum <- !(stringr::str_detect(colnames(texTable), "Parameter"))
# texTable[idxNum] <- lapply(texTable[idxNum], function(x) sprintf("%.2f", x))

# remove leading zero for numbers in (-1,1)
drop0 <- function(s) {
  gsub("^(-?)0\\.", "\\1.", s) 
}

# typeset effect sizes > boldThresh in bold + drop leading zero for those effect sizes
boldThresh <- 0.2
texTable[idxNum] <- apply(texTable[idxNum], c(1,2), function(x) {
  ifelse(x >= boldThresh, 
         paste0("\\textbf{", drop0(sprintf("%.2f", x)), "}"), 
         sprintf("%.2f", x))
})

# drop zero for effect sizes <  bold thresh
texTable[idxNum] <- lapply(texTable[idxNum], drop0)

print(xtable::xtable(texTable, type = "latex"),
      file = paste0(plotFolder, "/ANOVA/betweenANOVA_", "ENET_LM","_", eta2Thresh, ".tex"), 
      include.rownames = FALSE, sanitize.text.function = identity)
