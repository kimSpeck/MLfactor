# load parameters & custom functions 
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

library(ggplot2)

plotFolder <- "plots"
createFolder(plotFolder)

plotCheckFolder <- "checkGraphics"
createFolder(paste0(plotFolder, "/", plotCheckFolder))

grid <- params[["dgp"]][["grid"]]

shapeVec <- c(3, 4)
dodge_w <- 0.5
################################################################################
# calculate effect sizes for nuances
################################################################################
unique(grid$MIvio)
subGrid_N <- grid[grid$MIvio == "specific",]

# calculate variance of the outcome for all numbers of nuances
subGrid_N$varY <- (subGrid_N$betaT**2 * subGrid_N$varT + 
                     subGrid_N$miItems * subGrid_N$varS + params$dgp$sigmaEy) /
  (1-subGrid_N$corGO**2)

# calculate regression coefficient for the grouping variable 
subGrid_N$betaG <- 2*subGrid_N$corGO*sqrt(subGrid_N$varY)

# check multiple formulas & compare betaG results! 
identical(subGrid_N$betaG, 
          subGrid_N$corGO * (sqrt(subGrid_N$varY)/sqrt(params$dgp$varG)))

# correlation of the nuance (S) and the outcome
subGrid_N$corSY <-  (subGrid_N$betaS * subGrid_N$varS) / 
  (sqrt(subGrid_N$varS) * sqrt(subGrid_N$varY))

# (partial) correlation when controlling for the latent variable and the group
subGrid_N$corSY_partial <- sqrt(subGrid_N$varS) / sqrt(subGrid_N$miItems * subGrid_N$varS + params$dgp$sigmaEy)

# effect size marginal effect of a single nuance on the outcome
subGrid_N$dS_marginal <- 2 * subGrid_N$corSY / sqrt(1-subGrid_N$corSY**2)

# effect size marginal effect of a single nuance on the outcome
subGrid_N$dS_conditional <- 2 * subGrid_N$corSY_partial / sqrt(1-subGrid_N$corSY_partial**2)

subGrid_N$betaS.std <- subGrid_N$betaS * (sqrt(subGrid_N$varS) / sqrt(subGrid_N$varY))

ggplot(subGrid_N, 
       aes(x = factor(miItems), y = dS_marginal ,
           group = interaction(factor(varT), factor(R2)), shape = factor(R2), 
           color = factor(varT))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[F]^2)) +
  scale_shape_discrete(name = expression(alpha)) +
  labs(y = expression(d[MACS]), x = "number of manipulated Items") + 
  facet_grid(~ corGO)

ggplot(subGrid_N, 
       aes(x = factor(miItems), y = dS_conditional ,
           group = interaction(factor(varT), factor(R2)), shape = factor(R2), 
           color = factor(varT))) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4, alpha = 0.4) +
  geom_point(position = position_dodge(width = dodge_w), size = 3, alpha = 0.7) +
  scale_colour_grey(start = 0.5, end = 0.25, name = expression(sigma[F]^2)) +
  scale_shape_discrete(name = expression(alpha)) +
  labs(y = expression(d[MACS]), x = "number of manipulated Items") + 
  facet_grid(~ corGO)

subGrid_N$miItemsLabels <- paste0("# nuances: ", subGrid_N$miItems)
subGrid_N$varS <- factor(subGrid_N$varS)
subGrid_N$R2 <- factor(subGrid_N$R2)
subGrid_N$corGO <- factor(subGrid_N$corGO)

pD_marginal <- ggplot(subGrid_N[subGrid_N$miItems > 0, ],
       aes(x = corGO, y = dS_marginal, group = interaction(varS, R2),
           color = R2 , shape = varS)) +
  geom_jitter(size = 3, height = 0, width = 0.05) +
  scale_colour_grey(start = 0.8, end = 0.2, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  scale_x_discrete(name = expression(rho[GY])) +
  scale_y_continuous(name = paste0("Cohen", "'", "s d"),   
                     limits = c(0, 1.5),                   
                     breaks = seq(0, 1.5, by = 0.25),      
                     labels = function(x) sprintf("%.2f", x)) + 
  facet_grid( ~ miItemsLabels)

(pD_marginal <- themeFunctionPaper(pD_marginal))

# ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/nuanceEffectSizes_cohensD.png"),
#                 plot = pD_marginal,
#                 width = 13.63,
#                 height = 9.60,
#                 units = "in")

pD_conditional <- ggplot(subGrid_N[subGrid_N$miItems > 0, ],
       aes(x = corGO, y = dS_conditional, group = interaction(varS, R2),
           color = R2 , shape = varS)) +
  geom_jitter(size = 3, height = 0, width = 0.05) +
  scale_colour_grey(start = 0.8, end = 0.2, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  scale_x_discrete(name = expression(rho[GY])) +
  scale_y_continuous(name = paste0("Cohen", "'", "s d (conditional)"),   
                     limits = c(0, 1.5),                     
                     breaks = seq(0, 1.5, by = 0.25),        
                     labels = function(x) sprintf("%.2f", x)) + 
  facet_grid( ~ miItemsLabels)

(pD_conditional <- themeFunctionPaper(pD_conditional))

# ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/nuanceEffectSizes_cohensD_conditional.png"),
#                 plot = pD_conditional,
#                 width = 13.63,
#                 height = 9.60,
#                 units = "in")

# recalculate standardized regression coefficients 
#   based on analytical outcome instead of estimated outcome variance (99_checkSamples.R)
pBetaS_std <- ggplot(subGrid_N[subGrid_N$miItems > 0, ],
                      aes(x = corGO, y = betaS.std, group = interaction(varS, R2),
                          color = R2 , shape = varS)) +
  geom_jitter(size = 3, height = 0, width = 0.05) +
  scale_colour_grey(start = 0.8, end = 0.2, name = expression(R^2)) +
  scale_shape_discrete(name = expression(sigma[S]^{2})) +
  scale_x_discrete(name = expression(rho[GY])) +
  scale_y_continuous(name = expression(beta[S[std.]]), 
                     limits = c(0, 0.7),                     
                     breaks = seq(0, 0.8, by = 0.2),        
                     labels = function(x) sprintf("%.1f", x)) +  
  facet_grid( ~ miItemsLabels)

(pBetaS_std <- themeFunctionPaper(pBetaS_std))

# ggplot2::ggsave(filename = paste0(plotFolder, "/", plotCheckFolder, "/nuanceEffectSizes_betaSTD_analytical.png"),
#                 plot = pBetaS_std,
#                 width = 13.63,
#                 height = 9.60,
#                 units = "in")

# # save nuance effect sizes to file to create a plot with all effect size distributions
# #   across all simulation conditions
# resFolder <- "results"
# depFolder <- "dependentMeasures"
# createFolder(paste0(resFolder, "/", depFolder))
# save(subGrid_N, file = paste0(resFolder, "/", depFolder, "/nuanceEffectSizes.rda"))
