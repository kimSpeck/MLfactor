# load parameters & custom functions
source("utils/setParam.R") # parameter values
source("utils/utilityFUNS.R")

# laod packages
library(tidyverse)

# folder structure
resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/performPerSample_differences.rda")
load(resFile)

# which conditions are included in the result file
unique(performMeasures$MIvio)
performMeasures <- performMeasures[!(performMeasures$MIvio %in% c("scalarReal0.4", "scalarReal0.6")),]
################################################################################
# settings
################################################################################
# performance measures and models for which the MCSE is calculated
# we reported difference in R² test performance, the MCSE for this measure is relevant
#   for the number of replications that we ran
mcseMeasures <- c("Rsquared_HOtest", "RMSE_HOtest")
unique(performMeasures$model)[grepl("diff", unique(performMeasures$model))]
mcseModels <- c("diff_LM", "diff_Ridge", "diff_ENET", 
                "diff_group_LM", "diff_group_Ridge", "diff_group_ENET")

################################################################################
# MCSE per simulated condition
################################################################################
# every combination of simulation parameters, training sample size N,
#   performance measure and model is one condition with nTrain replications
mcseStats <- performMeasures %>%
  filter(measure_type %in% mcseMeasures,
         model %in% mcseModels) %>%
  group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, measure_type, model) %>%
  summarize(nRep = sum(!is.na(values)), # number of replications that contribute to mean
            # condition mean (e.g., Delta R²)
            M = mean(values, na.rm = T), 
            # empirical SE across replications according to Morris et al. (2019)
            SD = sd(values, na.rm = T), 
            .groups = "drop") 

# number of replications that did yield a resuls (any NAs?)
any(mcseStats[mcseStats$nRep != 100, ]) # no NAs
# min difference is not relevant to determine the number of replications 
#   because we do not interpret lots of these differences (we expect them to be 0)
min(abs(mcseStats[mcseStats$measure_type == "Rsquared_HOtest", "M"]))

# precision criterion for condition means on R2 on holdout set: 
#   2 * MCSE <= .01 (smallest interpreted difference)
# mcseTarget = .005 corresponds to 2·MCSE ≤ .01
mcseTarget_R2 <- 0.005
# target more difficult to determine for RMSE because RMSE depends on the units of Y
#   because the variance of y changes, RMSE has different range of values for different
#   conditions
mcseTarget_RMSE <- 0.005
mcseStats$mcseTarget <- ifelse(mcseStats$measure_type == "Rsquared_HOtest", mcseTarget_R2, 
                               ifelse(mcseStats$measure_type == "RMSE_HOtest", mcseTarget_RMSE, NA)) 

# calculate monte carlo standard error in every conditon 
mcseStats <- mcseStats %>% 
  # How precisely does the simulation estimate E[Delta R²] for this condition?
  mutate(MCSE = SD / sqrt(nRep), # MCSE of the condition mean
         nRepTarget = ceiling((SD / mcseTarget)**2), # replications needed to reach mcseTarget
         aboveTarget = MCSE > mcseTarget)

# explicitly look at data that we analysed and reported
hoStats <- mcseStats %>%
  filter(measure_type == "Rsquared_HOtest", model %in% c("diff_ENET", "diff_group_ENET"))

# median, max, and number of cells above target per model
hoStats %>%
  group_by(model) %>%
  summarize(Mdn = median(MCSE), Max = max(MCSE),
            anyAbove = any(aboveTarget), nAbove = sum(aboveTarget))

# cells above target
hoStats %>% filter(aboveTarget)

# MCSE per cell, one panel per model
hoStats %>%
  group_by(model) %>% mutate(cell = row_number()) %>%
  ggplot(aes(cell, MCSE)) + geom_point() + facet_wrap(~ model)

# save MCSE for every condition
save(mcseStats, file = paste0(resFolder, "/", depFolder, "/mcsePerCondition.rda"))
