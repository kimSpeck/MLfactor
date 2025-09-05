# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

library(tidyverse)

resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/hyperPerSample.rda")
load(resFile)

plotFolder <- "plots"
createFolder(plotFolder)
plotPerformanceFolder <- "performance"
createFolder(paste0(plotFolder, "/", plotPerformanceFolder))

hyperSelection$miItemsPlot <- ifelse(hyperSelection$miItems %in% c(2, 4), "P/2", hyperSelection$miItems)
hyperSelection$N <- factor(hyperSelection$N, levels = c(100, 300, 1000))

# How often is ENET or itemENET essentially a ridge regression, due to alpha == 0?
hyperSelection[which(hyperSelection$model == "itemENET" & 
        hyperSelection$alpha == 0),]

countAlpha <- hyperSelection %>% 
  group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
  count(alpha, .drop = FALSE) 

# check
all(countAlpha[countAlpha$model == "Ridge", "n"] == 100)
# remove ridge (no information in this data)
countAlpha <- countAlpha[countAlpha$model != "Ridge", ]

# lots of conditions with all samples alpha = 0
View(countAlpha[which(countAlpha$n == 100), ])

ggplot(countAlpha[countAlpha$model == "itemENET" & 
                    countAlpha$MIvio == "specific" & 
                    countAlpha$varT == 0.25, ], 
       aes(x = factor(alpha), y = n, fill = miItemsPlot)) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_y_continuous(limits = c(-10, 100), breaks = seq(0, 100, 50)) +
  facet_grid(cAlpha + corGO ~ N + P + R2) +
  ylab("Count per level")

# How strong is penalty? lambda values? 
countLambda <- hyperSelection %>% 
  group_by(MIvio, miItemsPlot, P, Q, cAlpha, varT, corGO, R2, N, model, predSet) %>% 
  count(lambda, .drop = FALSE) 

ggplot(countLambda[countLambda$model == "itemENET" & 
                     countLambda$varT == 0.25, ], 
       aes(x = factor(lambda), y = n, fill = interaction(miItemsPlot, MIvio))) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_y_continuous(limits = c(-10, 100), breaks = seq(0, 100, 50)) +
  facet_grid(cAlpha + corGO ~ N + P + R2) +
  ylab("Count per level")
