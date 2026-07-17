# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

# load R packages
library(ggplot2)

# result object structure (plots, tables, ...)
plotFolder <- "plots"
createFolder(plotFolder)

plotCheckFolder <- "checkGraphics"
createFolder(paste0(plotFolder, "/", plotCheckFolder))

# path to folder with check data
checkDataFolder <- "checkData"
conditionFolders <- dir(checkDataFolder) # every condition folder 

# rename
fileVec <- list.files(checkDataFolder, pattern = "^checkData_", 
                      include.dirs = FALSE, recursive = FALSE)
# renameFiles <- sub("MIvioNF_", "MIvioNF", fileVec)
# file.rename(paste0(checkDataFolder, "/", fileVec), paste0(checkDataFolder, "/", renameFiles))

# function to parse file name to extract condition parameters
parseFilename <- function(x) {
  # pattern captures signed integers or decimals after each tag
  pat <- paste0(
    "^.*_([^_]+)",                      
    "_miItems([+-]?\\d+)",
    "_P([+-]?\\d+)",
    "_Q([+-]?\\d+)",
    "_cAlpha([+-]?\\d*\\.?\\d+)",
    "_varT([+-]?\\d*\\.?\\d+)",
    "_corGO([+-]?\\d*\\.?\\d+)",
    "_R2([+-]?\\d*\\.?\\d+)",
    "_b0([+-]?\\d*\\.?\\d+)\\.rda$"
    #"_R2([+-]?\\d*\\.?\\d+)\\.rda$"
  )
  m <- stringr::str_match(x, pat)
  if (is.na(m[1,1])) stop("Pattern not matched.")
  out <- as.data.frame(m[, -1, drop = FALSE], stringsAsFactors = FALSE)
  # names(out) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2")
  names(out) <- c("MIvio","miItems","P","Q","cAlpha","varT","corGO","R2", "b0")
  # convert numeric fields (leave MIvio as character)
  num_cols <- setdiff(names(out), "MIvio")
  out[num_cols] <- lapply(out[num_cols], as.numeric)
  out["MIvio"] <- substring(out["MIvio"], 6, nchar(out["MIvio"]))
  out
}
# parseFilename(iFile)

# path to folder with check data
checkDataFolder <- "checkData"
conditionFiles <- dir(checkDataFolder) # every condition folder 
# conditionFiles <- conditionFiles[!grepl("^checkData_MIvioNF", conditionFiles)]

conditionFiles <- conditionFiles[grepl("^checkDataMI_", conditionFiles)]

# prepare result objects as a list with entry for every condition
# one value per sample
sampleStats <- vector(mode = "list", length = length(conditionFiles)) 
# P {4, 8} values per sample; there are NAs!
variableStats <- vector(mode = "list", length = length(conditionFiles))
miTestStats <- vector(mode = "list", length = length(conditionFiles))

# iterate through all conditions to extract data into one object
for (iFolder in seq_len(length(conditionFiles))) {
  # get file name and load data
  # iFile <- list.files(paste0(checkDataFolder, "/", conditionFiles[iFolder]))
  iFile <- paste0(checkDataFolder, "/", conditionFiles[iFolder])
  load(paste0(checkDataFolder, "/", conditionFiles[iFolder]))
  
  # vector of item names
  itemNames <- sapply(seq_len(as.numeric(data[[1]]$condition[["P"]])), function(iP) {
    paste0("X", iP)
  })
  
  # extract simualted conditions
  conditionVars <- parseFilename(iFile) # from file name
  condVec <- which(params$dgp$grid$P == data[[1]]$condition[["P"]] &
                     params$dgp$grid$cAlpha == data[[1]]$condition[["cAlpha"]] &
                     params$dgp$grid$varT == data[[1]]$condition[["varT"]] &
                     params$dgp$grid$R2 == data[[1]]$condition[["R2"]],)
  
  # construct results matrix for scalar values describing the sample
  #   condition vars
  #   Correlation with group, Cronbach's Alpha, simulated vs. actual parameters 
  tmpMat <- matrix(NA, ncol = 29)
  tmpMat[seq_len(length(conditionVars))] <- unlist(conditionVars)
  tmpMat[length(conditionVars) + seq_len(length(data[[1]]$trueVShat))] <- unlist(data[[1]]$trueVShat)
  tmpMat[length(conditionVars) + length(data[[1]]$trueVShat) + 
           seq_len(length(data[[1]]$R2_fit))] <- unlist(data[[1]]$R2_fit)
  tmpMat[length(conditionVars) + length(data[[1]]$trueVShat) + 
           length(data[[1]]$R2_fit) + 1] <- unlist(data[[1]]$var_hat[["Y"]])
  tmpMat[length(conditionVars) + length(data[[1]]$trueVShat) + 
           length(data[[1]]$R2_fit) + 2] <- unlist(data[[1]]$var_hat[["group"]])
  
  sampleStats[[iFolder]] <- tmpMat
  
  # construct results matrix for vectors of values describing the sample
  #   factor loadings, standardized factor loadings, residual variances, variances of manifest variables, R² of items (by factor)
  tmpMat <- matrix(NA, ncol = 19, nrow = 5)
  tmpMat[1, ] <- c(unlist(conditionVars), "lambda", params$dgp$grid[condVec[1], "flT"],
                   if (conditionVars["P"] == 4) {
                     c(data[[1]][["lambda_hat"]], rep(NA, 4))
                   } else {data[[1]][["lambda_hat"]]})
  tmpMat[2, ] <- c(unlist(conditionVars), "lambda.std", 0,
                   if (conditionVars["P"] == 4) {
                     c(t(data[[1]][["flT.std"]]), rep(NA, 4))
                   } else {t(data[[1]][["flT.std"]])})
  tmpMat[3, ] <- c(unlist(conditionVars), "theta", 
                   params$dgp$grid[condVec[1], "varStable"] + params$dgp$grid[condVec[1], "varE"],
                   if (conditionVars["P"] == 4) {
                     c(data[[1]][["theta_hat"]], rep(NA, 4))
                   } else {data[[1]][["theta_hat"]]})
  tmpMat[4, ] <- c(unlist(conditionVars), "varItems", params$dgp$varItems,
                   if (conditionVars["P"] == 4) {
                     c(data[[1]][["var_hat"]][itemNames], rep(NA, 4))
                   } else { data[[1]][["var_hat"]][itemNames]})
  tmpMat[5, ] <- c(unlist(conditionVars), "R2_itemsT", 
                   params$dgp$grid[condVec[1], "flT"]**2 * params$dgp$grid[condVec[1], "varT"],
                   if (conditionVars["P"] == 4) {
                     c(data[[1]][["R2_itemsT"]], rep(NA, 4))
                   } else {data[[1]][["R2_itemsT"]]})
  variableStats[[iFolder]] <- tmpMat
  
  # construct results matrix for measurement invariance evaluation
  tmpMat <- matrix(NA, ncol = 25) # 8 condition values, 16 measures
  tmpMat[seq_len(length(conditionVars))] <- unlist(conditionVars)
  tmpMat[length(conditionVars) + (1:8)] <- unlist(data[[1]]$miTest)[1:8]
  tmpMat[(length(conditionVars) + 8) + (1:8)] <- if (conditionVars["P"] == 4) {
    c(unlist(data[[1]]$miTest)[9:12], rep(NA, 4))
  } else {unlist(data[[1]]$miTest)[9:16]}
  
  miTestStats[[iFolder]] <- tmpMat
}

################################################################################
# check measurement invariance manipulations
################################################################################
miTestStats <- as.data.frame(do.call(rbind, miTestStats))
colnames(miTestStats) <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2", "b0",
                           "deltaChiSquare_metric", "deltaChiSquare_scalar", 
                           "deltaCFI_metric", "deltaRMSEA_metric", "deltaSRMR_metric", 
                           "deltaCFI_scalar", "deltaRMSEA_scalar", "deltaSRMR_scalar", 
                           paste0("dMACS_X", 1:max(params$dgp$P)))


numCols <- colnames(miTestStats)[2:length(colnames(miTestStats))]
miTestStats[numCols] <- sapply(miTestStats[numCols], as.numeric)

range(miTestStats$deltaChiSquare_metric, na.rm = T)
range(abs(miTestStats$deltaCFI_metric), na.rm = T)

range(miTestStats$deltaChiSquare_scalar, na.rm = T)
range(abs(miTestStats$deltaCFI_scalar), na.rm = T)

dodge_w <- 1
# plot delta CFI for scalar non-invariance
ggplot(miTestStats[miTestStats$MIvio == "scalar", ], 
       aes(x = factor(b0), y = abs(deltaCFI_scalar),
           group = interaction(factor(varT), factor(cAlpha)), shape = factor(cAlpha), 
           color = factor(varT))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[F]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = expression(DeltaCFI), x = expression(R[sim]^2)) + 
  facet_grid( ~ P)

# plot dMACS for metric & scalar non-invariance
subDMACS <- tidyr::pivot_longer(miTestStats, cols = starts_with("dMACS"), 
                                names_to = "item", names_prefix = "dMACS_", 
                                values_to = "dMACS", values_drop_na = TRUE)

manipulatedItems <- c(which(subDMACS$P == 4 & subDMACS$miItems == 1 & subDMACS$item == "X4"), 
                      which(subDMACS$P == 8 & subDMACS$miItems == 1 & subDMACS$item == "X8"),
                      which(subDMACS$P == 4 & subDMACS$miItems == 2 & subDMACS$item %in% c("X3", "X4")),
                      which(subDMACS$P == 8 & subDMACS$miItems == 4 & subDMACS$item %in% c("X5", "X6", "X7", "X8")))

(pDMACS <- ggplot(subDMACS[manipulatedItems,], 
                 aes(x = factor(b0), y = abs(dMACS),
                     group = interaction(factor(varT), factor(cAlpha)), shape = factor(cAlpha), 
                     color = factor(varT))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[F]^2)) +
  scale_shape_discrete(name = expression(alpha)) +
  labs(y = expression(d[MACS]), x = "number of manipulated Items") + 
  facet_grid(~ MIvio)
)
pDMACS <- themeFunctionPaper(pDMACS)

# ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/checkMI_dMACS.png"),
#                 plot = pDMACS,
#                 width = 13.63,
#                 height = 9.60,
#                 units = "in")

ggplot(subDMACS[-manipulatedItems,], 
       aes(x = factor(miItems), y = abs(dMACS),
           group = interaction(factor(varT), factor(cAlpha)), shape = factor(cAlpha), 
           color = factor(varT))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[F]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = expression(DeltaCFI), x = expression(R[sim]^2)) + 
  facet_grid( ~ MIvio)
################################################################################
# check simulation parameters 
################################################################################
# check parameters with exactly one value per sample
sampleStats <- as.data.frame(do.call(rbind, sampleStats))
colnames(sampleStats) <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2", "b0",
                           "corGO_sim", "corGO_hat", "sigmaEy", "sigmaEy_hat", "betaT", "betaT_hat",
                           "R2_sim", "R2_hat", "cAlpha_sim", "cAlpha_myFUN", "cAlpha_semTools",
                           "betaT_hat_g1", "betaT_hat_g2", "R2_hat_g1", "R2_hat_g2",
                           "R2_items", "R2_score", "R2_avg", "varY_hat", "varG_hat")

# change variable types
numCols <- colnames(sampleStats)[2:length(colnames(sampleStats))]
sampleStats[numCols] <- sapply(sampleStats[numCols],as.numeric)
facCols <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2")
sampleStats[facCols] <- sapply(sampleStats[facCols],as.factor)
sampleStats$miItemsPlot <- ifelse(sampleStats$miItems %in% c(2, 4), "P/2", sampleStats$miItems)

# check in which conditions the SEM did not converge and therefore produced NA
# 16 out of 432 is fine; ignore NAs and still draw valie conclusions regarding sanity checks 
naIdx <- which(apply(sampleStats, 1, anyNA))
length(naIdx) 
# View(sampleStats[naIdx,])

# check simulated ...
# ... variance of the grouping variable
#     -> solid at 0.25 for every condition
plot(sampleStats$varG_hat)

# check simulated ...
# ... variance of the outcome variable
plot(sampleStats$varY_hat)

################################################################################
# plot standardized factor loadings
################################################################################
variableStats <- as.data.frame(do.call(rbind, variableStats))
colnames(variableStats) <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2", "b0",
                             "measure", "true", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")

# change variable types
numCols <- c("true", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")
variableStats[numCols] <- sapply(variableStats[numCols], as.numeric)
facCols <- c("MIvio", "miItems", "P", "Q", "cAlpha", "varT", "corGO", "R2", "b0")
variableStats[facCols] <- sapply(variableStats[facCols], as.factor)

itemNames <- c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")
variableStats$mean <- rowMeans(variableStats[, itemNames], na.rm = T)
variableStats$sd <- matrixStats::rowSds(as.matrix(variableStats[, itemNames], ncols = 8), na.rm = T)

unique(variableStats$measure)

variableStats$corNF <- ifelse(grepl("NF", variableStats$MIvio), 
                              substr(variableStats$MIvio, 3, nchar(as.character(variableStats$MIvio[1]))), 
                              0)
variableStats$MIvio <- ifelse(grepl("NF", variableStats$MIvio), "NF", variableStats$MIvio)
variableStats$miItemsPlot <- ifelse(variableStats$miItems > 1, "P/2", variableStats$miItems)

# check parameter estimates averaged across P items
#   ... only for miItems = 0!
variableStats$diff <-  variableStats$true - variableStats$mean

# item variances (manifest variable); should be fixed to 1
ggplot(variableStats[variableStats$miItems == "0" &
                       variableStats$measure == "varItems", ],
       aes(x = cAlpha, y = diff, group = interaction(P, R2),
           color = interaction(P, R2), fill = interaction(P, R2))) +
  geom_jitter()



# lambda
ggplot(variableStats[variableStats$miItems == "0" &
                       variableStats$measure == "lambda", ],
       aes(x = R2, y = diff, group = P,
           color = P, fill = P)) +
  geom_jitter() +
  facet_grid(cAlpha ~ varT)

# theta
ggplot(variableStats[variableStats$miItems == "0" &
                       variableStats$measure == "theta", ],
       aes(x = cAlpha, y = diff, group = P,
           color = P, fill = P)) +
  geom_jitter()


# R2 true score on items
ggplot(variableStats[variableStats$miItems == "0" &
                       variableStats$measure == "R2_itemsT", ],
       aes(x = cAlpha, y = diff, group = interaction(P, R2),
           color = interaction(P, R2), fill = interaction(P, R2))) +
  geom_jitter()

# examine standardized factor loadings
variableStats$miItemsLevels <- ifelse(variableStats$miItems > 1, 2, variableStats$miItems)
variableStats$MIvio <- factor(variableStats$MIvio,
                              levels = c("metric", "scalar", "specific", "NF"),
                              labels = c("metric", "scalar", "nuance", "NF"))
ggplot(variableStats[variableStats$measure == "lambda.std", ],
       aes(x = miItemsLevels, y = mean, group = P,
           color = P, shape = P)) +
  geom_jitter() +
  facet_grid(cAlpha ~  MIvio)

# long format of items
stdLambdaSub <- variableStats[variableStats$measure == "lambda.std", ]
Xcols <- grep("^X\\d+$", names(stdLambdaSub), value = TRUE)
stdLambdaSub_long <- tidyr::pivot_longer(stdLambdaSub, Xcols,
                                         names_to = "item", values_to = "flT.std")
stdLambdaSub_long$manipulated <- ifelse(as.numeric(stdLambdaSub_long$miItems) > 0 &
                                          as.numeric(substr(stdLambdaSub_long$item, 2, 2)) >
                                          as.numeric(stdLambdaSub_long$P) - as.numeric(stdLambdaSub_long$miItems), 1, 0)

unique(stdLambdaSub_long$manipulated)
stdLambdaSub_long$manipulated <- factor(stdLambdaSub_long$manipulated)
stdLambdaSub_long$miItemsLevels <- factor(stdLambdaSub_long$miItemsLevels,
                                          levels = c(0, 1, 2),
                                          labels = c("0", "1", "50%"))

(stdFL_perItem <- ggplot(stdLambdaSub_long,
                         aes(x = miItemsLevels, y = flT.std, group = P,
                             color = manipulated, shape = P, fill = manipulated)) +
    geom_jitter() +
    facet_grid(cAlpha ~  MIvio) +
    scale_color_manual(values = c("#D55E00", "#0072B2")) +
    scale_y_continuous(breaks = scales::breaks_width(0.1)) +
    ylab("Standardized Factor Loadings") +
    xlab("Number of Manipulated Items"))

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

(stdFL_perItem <- themeFunctionPaper(stdFL_perItem))
ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/stdFactorLoadings_perItem_appendix.png"),
                plot = stdFL_perItem,
                width = 11.66,
                height = 9.12,
                units = "in")


# mean for MI manipulated vs. unmanipulated items
stdLambdaSub <- variableStats[variableStats$measure == "lambda.std", ]
Xcols <- grep("^X\\d+$", names(stdLambdaSub), value = TRUE)
get_means <- function(P, miItems, data) {
  
  P <- as.integer(P); miItems <- as.integer(miItems)
  
  if (miItems == 0) {
    manip_idx   <- c(integer(0))
    unmanip_idx <- c(paste0("X", 1:P))
  } else if (miItems == 1) {
    manip_idx   <- c(paste0("X", P))
    unmanip_idx <- c(paste0("X", 1:(P-1)))
  } else { # m > 1
    manip_idx   <- c(paste0("X", (P - miItems + 1):P))
    unmanip_idx <- c(paste0("X", 1:(P - miItems)))
  }
  
  manip_mean   <- if (length(manip_idx)   == 0) NA_real_ else mean(data[manip_idx],   na.rm = TRUE)
  unmanip_mean <- if (length(unmanip_idx) == 0) NA_real_ else mean(data[unmanip_idx], na.rm = TRUE)
  
  c(manip = manip_mean, unmanip = unmanip_mean)
}

res <- t(mapply(
  FUN = get_means,
  P = stdLambdaSub$P,
  miItems = stdLambdaSub$miItems,
  data = stdLambdaSub)
)

stdLambdaSub$manipulated <- res[,"manip"]
stdLambdaSub$unManipulated <- res[,"unmanip"]

# varT and R2 without any influence
ggplot(variableStats[variableStats$measure == "lambda.std", ],
       aes(x = miItemsLevels, y = mean, group = P,
           color = P, shape = P, fill = manipulated)) +
  geom_jitter() +
  facet_grid(cAlpha ~  MIvio)


# ################################################################################
# # preliminary R2 comparisons (items vs. scores)
# ################################################################################
# # initial plots to compare R2 for items vs. scores (in linear regression)
# plot(sampleStats[,"R2_items"], sampleStats[,"R2_score"])
# # difference between items and score R2 (positive values: items > scores)
# sampleStats$R2_diff <- sampleStats$R2_items - sampleStats$R2_score
# sampleStats$R2_relDiff <- (sampleStats$R2_items - sampleStats$R2_score)/sampleStats$R2_sim 
# 
# sampleStats$MIvio <- factor(sampleStats$MIvio, 
#                             levels = c("specific", "metric", "scalar"))
# 
# ggplot(sampleStats[sampleStats$cAlpha == 0.6 &
#                      sampleStats$varT == 0.25, ], 
#        aes(x = MIvio, y = R2_diff, group = interaction(miItemsPlot, P), 
#            color = miItemsPlot, fill = miItemsPlot,
#            linetype = P)) +
#   geom_point() + 
#   geom_line() + 
#   facet_grid(R2 ~ corGO)
# 
# ggplot(sampleStats[sampleStats$cAlpha == 0.6 &
#                      sampleStats$varT == 0.5, ], 
#        aes(x = MIvio, y = R2_relDiff, group = interaction(miItemsPlot, P), 
#            color = miItemsPlot, fill = miItemsPlot,
#            linetype = P)) +
#   geom_point() + 
#   geom_line() + 
#   facet_grid(R2 ~ corGO)
# 
