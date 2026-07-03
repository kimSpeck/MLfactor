# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

grid <- params[["dgp"]][["grid"]]

library(ggplot2)
library(psych)

# folder to save plots in
plotFolder <- "plots"
createFolder(plotFolder)

plotCheckFolder <- "checkGraphics"
createFolder(paste0(plotFolder, "/", plotCheckFolder))

# load data
resFolder <- "results"
depFolder <- "dependentMeasures"
createFolder(paste0(resFolder, "/", depFolder))

load(paste0(resFolder, "/", depFolder, "/nuanceEffectSizes.rda"))
load(paste0(resFolder, "/", depFolder, "/miEffectSizes_dMACS.rda"))

# dMACS data frame: go to long format
subDMACS <- tidyr::pivot_longer(miTestStats, cols = starts_with("dMACS"), 
                                names_to = "item", names_prefix = "dMACS_", 
                                values_to = "dMACS", values_drop_na = TRUE)

# mark rows with manipulated items
manipulatedItems <- c(which(subDMACS$P == 4 & subDMACS$miItems == 1 & subDMACS$item == "X4"), 
                      which(subDMACS$P == 8 & subDMACS$miItems == 1 & subDMACS$item == "X8"),
                      which(subDMACS$P == 4 & subDMACS$miItems == "50%" & subDMACS$item %in% c("X3", "X4")),
                      which(subDMACS$P == 8 & subDMACS$miItems == "50%" & subDMACS$item %in% c("X5", "X6", "X7", "X8")))

# informative effect sizes: 
#   - dMACS for metric and scalar for manipulated items only
#   - dMACS uninformative for nuances!
#   - instead, Cohen's d informative for nuances
# therefore, ...
#   ... remove unmanipulated items for metric and scalar 
#   ... replace nuance information with Cohen's d

subDMACS <- subDMACS[manipulatedItems, ]
unique(subDMACS$MIvio)
subDMACS <- subDMACS[subDMACS$MIvio != "nuance", ]
subDMACS$miItems <- ifelse(subDMACS$miItems == 1, subDMACS$miItems, 
                           ifelse(subDMACS$miItems == "50%" & subDMACS$P == 4, 
                                  2, 4))

rmColsIDX <- grep("^delta", colnames(subDMACS)) 
subDMACS <- subDMACS[, -c(rmColsIDX)]

subGrid_N$item <- NA
unique(subGrid_N$MIvio)
subGrid_N$MIvio <- "nuance"  
colnames(subDMACS)[colnames(subDMACS) == "dMACS"] <- "d"
colnames(subGrid_N)[colnames(subGrid_N) == "dS_marginal"] <- "d"

commonColumns <- intersect(names(subGrid_N), names(subDMACS))
effectSizeData <- merge(subDMACS, subGrid_N[, commonColumns], by = commonColumns, all = TRUE)
effectSizeData <- effectSizeData[effectSizeData$miItems > 0, ]
unique(effectSizeData$MIvio)

psych::describeBy(d ~ MIvio + R2 + corGO + miItems, data = effectSizeData)

################################################################################
# plot effect sizes (overview)
################################################################################
shapeVec <- c(3, 4)
dodge_w <- 0.5

effectSizeData$varS <- factor(ifelse(effectSizeData$varT == 0.25, 0.5, 0.25))
effectSizeData$miItemsLabels <- paste0("# Items: ", effectSizeData$miItems)

# metric also differes depending on Cronbachs alpha which is not visualized by 
#     any graphic feature 
pEffectSizes <- ggplot(effectSizeData, 
       aes(x = factor(corGO), y = d,
           group = interaction(factor(varS), factor(R2)), shape = factor(varS), 
           color = factor(R2))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  labs(y = expression(d), x = expression(rho[GY])) + 
  facet_grid(miItemsLabels ~ MIvio, )

(pEffectSizes <- themeFunctionPaper(pEffectSizes))

# ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/allEffectSizes_inD.png"),
#                 plot = pEffectSizes,
#                 width = 13.63,
#                 height = 9.60,
#                 units = "in")

################################################################################
# plot effect sizes (paper)
################################################################################
effectSizeData$MIvio <- factor(effectSizeData$MIvio, 
                               levels = c("metric", "scalar (S)", "scalar (M)", "scalar (XL)", "nuance"), 
                               labels = c("Metric", "scalar (S)", "scalar (M)", "Scalar", "Nuance"))
pPaperEffectSizes <- ggplot(effectSizeData[!(effectSizeData$MIvio %in% c("scalar (S)", "scalar (M)")), ], 
                       aes(x = factor(corGO), y = d,
                           group = interaction(factor(varS), factor(R2)), shape = factor(varS), 
                           color = factor(R2))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  labs(y = expression(d), x = expression(rho[GY])) + 
  facet_grid(miItemsLabels ~ MIvio, )

(pPaperEffectSizes <- themeFunctionPaper(pPaperEffectSizes))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/paperEffectSizes_inD.png"),
                plot = pPaperEffectSizes,
                width = 13.63,
                height = 9.60,
                units = "in")

aggregate(effectSizeData[!(effectSizeData$MIvio %in% c("scalar (S)", "scalar (M)")), ], 
          d ~ MIvio + varT, mean)
aggregate(effectSizeData[!(effectSizeData$MIvio %in% c("scalar (S)", "scalar (M)")), ], 
          d ~ MIvio + varT, sd)
aggregate(effectSizeData[!(effectSizeData$MIvio %in% c("scalar (S)", "scalar (M)")), ], 
          d ~ MIvio + varT, median)
aggregate(effectSizeData[!(effectSizeData$MIvio %in% c("scalar (S)", "scalar (M)")), ], 
          d ~ MIvio + varT, range)

