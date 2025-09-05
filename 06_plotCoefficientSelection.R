# load parameters & custom functions 
source("utils/setParam.R") # parameter values
params$sim$nTrain <- 100

source("utils/utilityFUNS.R")

# load packages
library(tidyverse)
library(see) # violin half plots 

resFolder <- "results"
depFolder <- "dependentMeasures"

resFile <- paste0(resFolder, "/", depFolder, "/coefPerSample.rda")
load(resFile)

plotFolder <- "plots"
createFolder(plotFolder)
plotPerformanceFolder <- "performance"
createFolder(paste0(plotFolder, "/", plotPerformanceFolder))

coefMatrix$miItemsPlot <- ifelse(coefMatrix$miItems %in% c(2, 4), "P/2", coefMatrix$miItems)
coefMatrix$N <- factor(coefMatrix$N, levels = c(100, 300, 1000))
coefMatrix$MIvio <- factor(coefMatrix$MIvio, levels = c("specific", "metric", "scalar"))

# Does ENET or itemENET keep items with MI violations or specific effects more often 
#   than without manipulation? 
coef_pDiffSub <- coefMatrix[coefMatrix$cAlpha == 0.8 & # for ENET
                        coefMatrix$model %in% c("ENET", "itemENET") & # for ENET
                        coefMatrix$varT == 0.25 & # condition with larger MI effect
                        coefMatrix$R2 == 0.8, ] # due to scaling numerically bigger

coef_pDiffSub <- tidyr::pivot_longer(coef_pDiffSub, 
                               cols = !c(MIvio, miItems, P, Q, cAlpha, varT, corGO,
                                         R2, sample, N, model, predSet, miItemsPlot),
                               names_to = "coefficient", values_to = "value")

ggplot(coef_pDiffSub[coef_pDiffSub$P == 4 &
                       coef_pDiffSub$predSet == "comb" &
                       coef_pDiffSub$miItemsPlot == 1 &
                       coef_pDiffSub$coefficient %in% c("X1", "X2", "X3", "X4", "score"),], 
       aes(x = interaction(N, model), y = value, fill = coefficient, color = coefficient)) +
  coord_flip() + 
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients")+
  facet_grid(corGO ~ MIvio)
  

ggplot(testSub[testSub$MIvio == "metric" &
                 testSub$predSet == "comb" & 
                 testSub$coefficient %in% c("X1", "X2", "X3", "X4", "score"),], 
       aes(x = N, y = value, fill = coefficient, color = coefficient)) +
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients")+
  facet_wrap(~ cAlpha + corGO) +
  coord_flip() 

ggplot(testSub[testSub$MIvio == "scalar" &
                 testSub$predSet == "items" & 
                 testSub$coefficient %in% c("X1", "X2", "X3", "X4"),], 
       aes(x = N, y = value, fill = coefficient, color = coefficient)) +
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients")+
  facet_wrap(~ cAlpha + corGO) +
  coord_flip() 

ggplot(testSub[testSub$MIvio == "specific" &
                 testSub$predSet == "items" & 
                 testSub$coefficient %in% c("X1", "X2", "X3", "X4"),], 
       aes(x = N, y = value, fill = coefficient, color = coefficient)) +
  geom_violinhalf(width=2.1, alpha = 0.5) +
  xlab("") +
  ylab("Regression Coefficients")+
  facet_wrap(~ cAlpha + corGO) +
  coord_flip() 

# Item 4 is the one with the manipulation
colVec <- c( "#006400", "#FF7F50", "#B22222") 
ggplot(testSub[testSub$MIvio == "metric" &
                  testSub$predSet == "items",], 
       aes(x = X1, y = X4, group = N, 
           color = N)) +
  #geom_abline(intercept = 0, slope = 0, "black", alpha = 0.5) +
  #geom_errorbar(aes(ymin = M - 2*SD, ymax = M + 2*SD), alpha = 0.2, width = 0.1) + 
  geom_point(size = 1.5, alpha = 0.5) +
  #geom_line(alpha = 0.4, linewidth = 0.75) +
  #scale_shape_manual(values = shapeVec) +
  scale_color_manual(values = colVec) +
  #scale_linetype_manual(values = lineTypeVec) +
  #ylab(paste0("difference in ", iMT, " (items - LM score)")) +
  facet_wrap(~ cAlpha + corGO)

# maybe remove ENET and itemENET if alpha = 0, because ENET would not have been able 
#   to do variable selection!
# Is there a difference between the manipulations in item selection?
