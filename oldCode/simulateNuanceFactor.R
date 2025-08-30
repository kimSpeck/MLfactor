# mini simulation nuance factor
# Mottus et al. formulieren in ihren inhaltlichen Beschreibungen dieses dgm
# nuancen bilden zusätzlichen Faktor ab (trait), der Outcome beeinflusst
# aber: auf Basis dieses dgms ist R^2 mit Items als Prädiktoren nicht besser als 
#       nur Faktor als Prädiktor
#   -> Suche das datengenerierende Modell, für das Performanz in R² mit Einzelitems >>> Performanz 
#     mit Faktor als Prädiktor wird?
# Odd: mit lambdaN = 0.8 & betaNY = 0.4 factor + items ist plötzlich besser als "trueR2" -> Why?!
# erweitere Mini-Simulation mit diesem dgm: Nuancenfaktor simulieren

library(lavaan)

# Ladungsstruktur: hohe Ladungen factor != nuance; wie wenn nur 1 Item?
setParam <- list()
setParam$N <- 10000
setParam$nOutcome <- 1
setParam$nItems <- 6
# Abweichungen in Lambda pro Gruppe dürfen in Summe nicht 0 sein (balanced DIF; 
#   see Robitzsch & Lüdtke, 2020)
setParam$lambdaFa <- c(1, 0.8, 0.2)
setParam$lambdaFb <- c(0.8, 0.8, 0.8)
setParam$lambdaFc <- c(1, 0.6, 0.2)
setParam$lambdaF <- list(lambdaFa = setParam$lambdaFa, 
                         lambdaFb = setParam$lambdaFb, 
                         lambdaFc = setParam$lambdaFc)


# systamtisch manipulieren: 
#     - Faktorladungen in Nuancenfaktor größer machen
#     - theta kleiner machen (Residualvarianzen in Items)
#     - Effekt von Nuancenfaktor/Facettenfaktor auf Outcome variieren
condGrid <- expand.grid(lambdaF = c("lambdaFa", "lambdaFb", "lambdaFc"),
                        lambdaN = c(0.2, 0.5, 0.8),
                        error = c(0.2, 0.4),
                        betaNY = c(0.4, 0.6, 0.8),
                        betaFY = c(0.6))

str(condGrid) # ! factors are interpreted as level numbers!; only character variables are interpreted by their name!
condGrid$lambdaF <- as.character(condGrid$lambdaF)

# create seed number for parallel cluster (reproducibility of generated data)
set.seed(8967369)
seedNum <- sample(1:999999, dim(condGrid)[1], replace = FALSE) 
condGrid$sampleSeed <- seedNum[1:dim(condGrid)[1]]

################################################################################
# simulate data
################################################################################
# sample data in parallel
sampleData <- function(lambdaF, lambdaN, error, betaNY, betaFY, sampleSeed){
  # matrix of factor loadings
  lambdaSEM <- cbind(factor = c(setParam$lambdaF[[lambdaF]][1], setParam$lambdaF[[lambdaF]][2], 
                                setParam$lambdaF[[lambdaF]][3], setParam$lambdaF[[lambdaF]][1], 
                                setParam$lambdaF[[lambdaF]][2], setParam$lambdaF[[lambdaF]][3], 
                                rep(0, setParam$nOutcome)),
                     nuance1 = c(lambdaN, lambdaN, lambdaN, 0, 0, 0, rep(0, setParam$nOutcome)),
                     fY = c(rep(0, setParam$nItems), rep(1, setParam$nOutcome)))
  nFactors <- ncol(lambdaSEM)
  
  # factors are not correlated 
  corSEM <- diag(nFactors) * 0.5
  # one additional nuance factor
  corSEM[upper.tri(corSEM, diag = F)] <- c(0,
                                           0, 0) 
  corSEM <- corSEM + t(corSEM)
  
  # path coefficient matrix
  betaMat <- matrix(c(0, 0, 0,
                      0, 0, 0,
                      betaFY, betaNY, 0), ncol = nFactors, byrow = T)
  
  # error
  thetaSEM <- diag(setParam$nItems + setParam$nOutcome) * error
  
  # calculate model implied covariance matrix
  sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% corSEM %*% 
    t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM
  
  # set mean values for variables
  set.seed(sampleSeed)
  M <- c(rnorm(mean = 3.1, sd = 0.6, sum(setParam$nItems)), 
         rnorm(mean = 10, sd = 3, setParam$nOutcome))
  
  # simulate date from multivariate normal distribution
  set.seed(sampleSeed)
  data <- mvtnorm::rmvnorm(setParam$N, mean = M, sigma = sigmaSEM)
  colnames(data) <- c(paste0("X", seq_len(setParam$nItems[1])), "y")
  data <- data.frame(data)
}

dataList <- do.call(mapply, c(FUN = sampleData, condGrid, SIMPLIFY = FALSE))
names(dataList) <- paste0(condGrid$lambdaF, 
                          "_lambdaN", condGrid$lambdaN, 
                          "_error", condGrid$error, 
                          "_betaNY", condGrid$betaNY, 
                          "_betaFY", condGrid$betaFY)

fileName <- "nuanceTestData.rda"
save(dataList, file = fileName)

################################################################################
# fit data
################################################################################
source("testModels.R")
load("nuanceTestData.rda")

# all SEMs
names(testSEMs)

# data names in list of data
listNames <- paste0(condGrid$lambdaF,
                    "_lambdaN", condGrid$lambdaN, 
                    "_error", condGrid$error, 
                    "_betaNY", condGrid$betaNY, 
                    "_betaFY", condGrid$betaFY)

predictivePerformance <- sapply(listNames, function(iData) {
  tmp.pp <- sapply(names(testSEMs), function(iSEM) {
    fit <- sem(testSEMs[[iSEM]], data = dataList[[iData]])
    inspect(fit, "rsquare")[["fY"]]
  }, simplify = T)  
  
  dataList[[iData]]$sumScore <- (dataList[[iData]]$X1 +  dataList[[iData]]$X2 +  
                                   dataList[[iData]]$X3 + dataList[[iData]]$X4 +  
                                   dataList[[iData]]$X5 +  dataList[[iData]]$X6)/6 
  fit4.lm <- lm(y ~ sumScore, data = dataList[[iData]])
  fit5.lm <- lm(y ~ sumScore + X1 + X2 + X3 + X4 + X5 + X6, 
                data = dataList[[iData]])
  
  tmp.pp <- c(tmp.pp,
              "sumScore" = summary(fit4.lm)$r.squared,
              "sumScore+items" = summary(fit5.lm)$r.squared)
  tmp.pp
}, simplify = T)

pP <- t(predictivePerformance)
pP <- data.frame(pP)
colnames(pP) <- c(names(testSEMs), "sumScore", "sumScore+items")

pP$measure <- rownames(pP)
pP$measure <- stringr::str_replace(pP$measure, "\\.[:digit:]{1,}$", "")

pP <- tidyr::separate(pP, measure, c("lambdaF", "lambdaN", "error", 
                                     "betaNY", "betaFY"), sep = "_")
pP$lambdaF <- stringr::str_sub(pP$lambdaF, start = 8L)
pP$lambdaN <- stringr::str_sub(pP$lambdaN, start = 8L)
pP$error <- stringr::str_sub(pP$error, start = 6L)
pP$betaNY <- stringr::str_sub(pP$betaNY, start = 7L)
pP$betaFY <- stringr::str_sub(pP$betaFY, start = 7L)
rownames(pP) <- NULL

pP[which(pP$items > pP$factor),]

################################################################################
# analyse data
################################################################################
library(ggplot2)
pP <- tidyr::pivot_longer(pP, cols = !c(lambdaF, lambdaN, error, betaNY, betaFY),
                          names_to = "measure", values_to = "rsquared")
pP$measure <- factor(pP$measure, 
                     levels = c("dgm", "factor", "items", "factorItems", "factorItem1",
                                "sumScore", "sumScore+items")) 


colValues = c("darkgreen", "darkblue", "darkred", "darkmagenta", "magenta",
              "cyan", "darkcyan")
ggplot(pP,
       aes(x = measure, y = rsquared, 
           group = interaction(lambdaF, lambdaN), colour = lambdaN,
           linetype = lambdaF)) +
   geom_point() +
   geom_line() +
   scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
   # geom_errorbar(aes(ymin = M - SE, ymax = M + SE), width=.2) +
   scale_color_manual(values = colValues) +
   # geom_hline(aes(yintercept = 0)) +
   facet_grid(error ~ betaNY, labeller = label_both) +
   ylab("") +
   xlab("pTrash (decreasing) x N (increasing)") +
   ggtitle("R^2: Training vs. Test performance") +
   theme(axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust=1),
         axis.title.x = element_text(size = 20),
         axis.title.y = element_text(size = 20),
         strip.text.x = element_text(size = 15),
         strip.text.y = element_text(size = 15))

# preliminary results: 
#   - predictive performance items > factor for:
#     betaNY = 0.6 x lambdaN = {0.5, 0.2} x error = 0.2
#     betaNY = 0.8 x lambdaN = {0.5, 0.2} x error = 0.2
#     betaNY = 0.8 x lambdaN = {0.5} x error = 0.4
# betaNY needs to be sufficiently high
# error needs to be sufficiently small
# paradoxical lambda influence: smaller lambda paradoxically better for items
