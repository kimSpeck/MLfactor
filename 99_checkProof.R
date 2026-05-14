# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

library(ggplot2)

plotFolder <- "plots"
createFolder(plotFolder)

plotCheckFolder <- "checkGraphics"
createFolder(paste0(plotFolder, "/", plotCheckFolder))

grid <- params[["dgp"]][["grid"]]
################################################################################
# covariance outcome -- residualized scalar non-invariant item
################################################################################
subGrid <- grid[grid$MIvio %in% c("scalar", "scalarReal0.4", "scalarReal0.6") &
                  grid$miItems == 1,]

# calculate var(y)
subGrid$varY <- (subGrid$betaT**2 * subGrid$varT + params$dgp$sigmaEy) / (1 - subGrid$corGO**2)

# check if analytic var(y) ~ var(ŷ): yes! (comparison with var(ŷ) plots)
(yPlot_MI <- ggplot(subGrid, 
                    aes(x = factor(corGO), y = varY, group = factor(R2), color = factor(R2))) +
    geom_jitter(size = 3) +
    scale_colour_grey(start = 0.8, end = 0.2, name = expression(R^2)))

(yPlot_MI <- themeFunctionPaper(yPlot_MI) + 
    labs(x = expression(r[GY]), y = expression(Var(Y))))

# calculate beta_G
subGrid$betaG <- subGrid$corGO * (sqrt(subGrid$varY)/sqrt(params$dgp$varG))

# calculate cov(Ip, Y)
subGrid$covIY <- (subGrid$betaT * subGrid$flT * subGrid$varT) + 
  (subGrid$betaG * subGrid$b0 * params$dgp$varG)

# calculate cov(Ip, S)
subGrid$covIS <- (subGrid$flT**2 * subGrid$varT) + 
 ((subGrid$b0**2 / subGrid$P) * params$dgp$varG) +
 ((1/subGrid$P) * subGrid$varE)

# calculate cov(Y, S)
subGrid$covYS <- (subGrid$betaT * subGrid$flT * subGrid$varT) + 
  (subGrid$betaG * (subGrid$b0/subGrid$P) * params$dgp$varG)

# calculate varS
subGrid$varScoreHetero <- (subGrid$flT**2 * subGrid$varT) + 
  ((subGrid$b0**2 / subGrid$P**2) * params$dgp$varG) + 
  ((1/subGrid$P**2) * (((subGrid$P-1) * (subGrid$varE + subGrid$varStable)) + subGrid$varE))

# # assuming homoscedastic error; this is wrong!
# subGrid$varScoreHomo <- (subGrid$flT**2 * subGrid$varT) + 
#   ((subGrid$b0**2 / subGrid$P**2) * params$dgp$varG) + 
#   (subGrid$varE/subGrid$P)

# calculate cov(Ip, Y|S) - partial covariance
subGrid$covIY_S_Hetero <- subGrid$covIY - ((subGrid$covIS * subGrid$covYS) / subGrid$varScoreHetero)
# # assuming homoscedastic error; this is wrong!
# subGrid$covIY_S_Homo <- subGrid$covIY - ((subGrid$covIS * subGrid$covYS) / subGrid$varScoreHomo)

# calculate the covariance of Y with the non-invariant item residualized by 
#   what the Score already explains in the item 
# calculate gamma 

subGrid$gamma <- subGrid$covIS / subGrid$varScoreHetero
# subGrid$gamma <- (subGrid$flT**2 * subGrid$varT + 
#                     (subGrid$b0**2 / subGrid$P) *  params$dgp$varG +
#                     subGrid$varE/subGrid$P) / 
#   (subGrid$flT**2 * subGrid$varT + 
#      (subGrid$b0**2 / subGrid$P**2) *  params$dgp$varG +
#      (1/subGrid$P**2) * (((subGrid$P-1) * (subGrid$varE + subGrid$varStable)) + subGrid$varE))
  
subGrid$covYresI <- subGrid$covIY - subGrid$gamma * subGrid$covYS

# plot results
shapeVec <- c(3, 4)
dodge_w <- 0.5
subGrid$corGO_lab <- paste0("r[GY]==", subGrid$corGO)
subGrid$R2_lab    <- paste0("R[sim]^2==", subGrid$R2)
subGrid$P_lab <- paste0("P == ", subGrid$P)
subGrid$cAlpha_lab <- paste0("alpha == ", subGrid$cAlpha)

covPlotScalar <- ggplot(subGrid, 
       aes(x = factor(R2), y = covIY_S_Hetero,
           group = interaction(factor(varStable), factor(P)), shape = factor(P), 
           color = factor(varStable))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[specific]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = "cov(Ip, Y | S)", x = expression(R[sim]^2)) +
  #facet_grid(cAlpha_lab ~ corGO_lab, labeller = label_parsed)
  facet_grid(MIvio ~ corGO_lab, labeller = label_parsed)

(covPlotScalar <- themeFunctionPaper(covPlotScalar))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/covIY_scalar.png"),
                plot = covPlotScalar,
                width = 13.63,
                height = 9.60,
                units = "in")


(covPlotScalar2 <- ggplot(subGrid, 
                        aes(x = factor(R2), y = covYresI,
                            group = interaction(factor(varStable), factor(P)), shape = factor(P), 
                            color = factor(varStable))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[specific]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = expression("cov(Y, Ip - "~beta[S]~"S)"), x = expression(R[sim]^2)) +
  facet_grid(cAlpha_lab ~ corGO_lab, labeller = label_parsed))

(covPlotScalar2 <- themeFunctionPaper(covPlotScalar2))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/covIY_scalar2.png"),
                plot = covPlotScalar2,
                width = 13.63,
                height = 9.60,
                units = "in")

################################################################################
# covariance outcome -- residualized metric non-invariant item
################################################################################
subGrid_M <- grid[grid$MIvio == "metric" &
                  grid$miItems == 1,]

# calculate var(y) [same as for scalar non-invariance]
subGrid_M$varY <- (subGrid_M$betaT**2 * subGrid_M$varT + params$dgp$sigmaEy) / (1 - subGrid_M$corGO**2)

# calculate beta_G [same as for scalar non-invariance]
subGrid_M$betaG <- subGrid_M$corGO * (sqrt(subGrid_M$varY)/sqrt(params$dgp$varG))

# calculate cov(Ip, Y)
subGrid_M$covIY <- (subGrid_M$betaT * subGrid_M$varT) * 
  (subGrid_M$flT + subGrid_M$b1 * params$dgp$meanG)

# calculate cov(Ip, S)
subGrid_M$covIS <- (subGrid_M$flT**2 * subGrid_M$varT) + 
  ((subGrid_M$b1**2 / subGrid_M$P) * (params$dgp$meanG**2 * subGrid_M$varT + 
                                        params$dgp$varG * subGrid_M$varT)) +
  (subGrid_M$flT * subGrid_M$b1 * (params$dgp$meanG * subGrid_M$varT) * (1 + (1/subGrid_M$P))) +
  ((1/subGrid_M$P) * subGrid_M$varE)
  

# calculate cov(Y, S)
subGrid_M$covYS <- (subGrid_M$betaT * subGrid_M$varT) * 
  (subGrid_M$flT + ((subGrid_M$b1 / subGrid_M$P) * params$dgp$meanG))

# calculate varS [items are not homoscedastic!]
subGrid_M$varScoreHetero <- (subGrid_M$flT**2 * subGrid_M$varT) +
  ((subGrid_M$b1**2 / subGrid_M$P**2) * (params$dgp$meanG**2 * subGrid_M$varT + 
                                          params$dgp$varG * subGrid_M$varT)) +
  ((1/subGrid_M$P**2) * (((subGrid_M$P-1) * (subGrid_M$varE + subGrid_M$varStable)) + subGrid_M$varE))
  

# calculate cov(Ip, Y|S) - partial covariance
subGrid_M$covIY_S_Hetero <- subGrid_M$covIY - ((subGrid_M$covIS * subGrid_M$covYS) / subGrid_M$varScoreHetero)

# calculate the covariance of Y with the non-invariant item residualized by 
#   what the Score already explains in the item 
# calculate gamma 
subGrid_M$gamma <- subGrid_M$covIS / subGrid_M$varScoreHetero
range(subGrid_M$gamma)

subGrid_M$covYresI <- subGrid_M$covIY - subGrid_M$gamma * subGrid_M$covYS

# plot results
shapeVec <- c(3, 4)
dodge_w <- 0.5
subGrid_M$corGO_lab <- paste0("r[GY]==", subGrid_M$corGO)
subGrid_M$R2_lab    <- paste0("R[sim]^2==", subGrid_M$R2)
subGrid_M$P_lab <- paste0("P == ", subGrid_M$P)
subGrid_M$cAlpha_lab <- paste0("alpha == ", subGrid_M$cAlpha)

covPlotMetric <- ggplot(subGrid_M, 
                        aes(x = factor(R2), y = covIY_S_Hetero,
                            group = interaction(factor(varStable), factor(P)), shape = factor(P), 
                            color = factor(varStable))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[specific]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = "cov(Ip, Y | S)", x = expression(R[sim]^2)) +
  facet_grid(cAlpha_lab ~ corGO_lab, labeller = label_parsed)

(covPlotMetric <- themeFunctionPaper(covPlotMetric))

idxVarE0 <- which(subGrid_M$varE == 0)
idxPlot <- which(subGrid_M$P == 4 & 
        subGrid_M$varStable == 0.5 & 
        subGrid_M$cAlpha == 0.8)

identical(idxVarE0, idxPlot)

subGrid_M$gamma[idxVarE0]
range(subGrid_M$gamma)
range(subGrid_M$covYresI)
which(subGrid_M$covYresI < 0) %in% idxVarE0

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/covIY_metric.png"),
                plot = covPlotMetric,
                width = 13.63,
                height = 9.60,
                units = "in")

(covPlotMetric2 <- ggplot(subGrid_M, 
                          aes(x = factor(R2), y = covYresI,
                              group = interaction(factor(varStable), factor(P)), shape = factor(P), 
                              color = factor(varStable))) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
    geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
    scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[specific]^2)) +
    scale_shape_discrete(name = "P") +
    labs(y = expression("cov(Y, Ip - "~beta[S]~"S)"), x = expression(R[sim]^2)) +
    facet_grid(cAlpha_lab ~ corGO_lab, labeller = label_parsed))

(covPlotMetric2 <- themeFunctionPaper(covPlotMetric2))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/covIY_metric2.png"),
                plot = covPlotMetric2,
                width = 13.63,
                height = 9.60,
                units = "in")

(covPlotMetric <- ggplot(subGrid_M, 
                        aes(x = factor(corGO), y = covIY_S_Hetero,
                            group = interaction(factor(varStable), factor(P)), shape = factor(P), 
                            color = factor(varStable))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[specific]^2)) +
  scale_shape_discrete(name = "P") +
  labs(y = "cov(Ip, Y | S)", x = expression(r[GY])) +
  facet_grid(cAlpha_lab ~ R2_lab, labeller = label_parsed))

(covPlotMetric <- themeFunctionPaper(covPlotMetric))

ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/covIY_metric_clear.png"),
                plot = covPlotMetric,
                width = 13.63,
                height = 9.60,
                units = "in")
