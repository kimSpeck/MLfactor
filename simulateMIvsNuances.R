# simulate data to check the assumption that measurement invariance could account
#   for some portion of the increase in explained variance due to nuance level based
#   prediction in comparison with domain/facet based prediction
# reanalyse data: test data for measurement invariance?
# this is a freely available NEO-PI dataset on which Mottus & Rozgonjuk demonstrated
#   age prediction based on nuances
#   https://osf.io/tbmh5/

# To Do:
# plot correlations for predictors
#   - only one factor
#   - additional nuance factor
#   - measurement invariance violation

################################################################################
# what is where in matrices when simulating SEMs
################################################################################
# simulate data from common factor model and throw data into lavaan 
#   check dimensionalities and entries of matrices when simulating data from SEM
N <- 10000
nItems <- c(3)
lambdaFull <- cbind(factor = c(0.8, 0.6, 0.6),
                    nuance1 = c(0.2, 0.2, 0))

# # with both nuance factor not enough parameter to fix in SEM
# lambdaFull <- cbind(factor = c(0.8, 0.6, 0.6),
#                     nuance1 = c(0.2, 0.2, 0),
#                     nuance2 = c(0, 0, -0.2))

## set factor correlations
# factors are not correlated 
#   nuance information is independent from facet information (this is the dgm, i.e., 
#   extra information in nuances that has nothing to do with what is already explained by facet)
nFactors <- ncol(lambdaFull)
corMat <- diag(nFactors) * 0.5
# one additional nuance factor
corMat[upper.tri(corMat, diag = F)] <- c(0) 

# # two additional nuance factors
# corMat[upper.tri(corMat, diag = F)] <- c(0, 
#                                          0, 0)
corMat <- corMat + t(corMat)

# error
theta <- diag(sum(nItems)) * 0.2

# calculate covariance matrix
# for simulating from common factor model
sigmaMat <- lambdaFull %*% corMat %*% t(lambdaFull) + theta

# set mean values for variables
set.seed(289301)
M <- rnorm(mean = 3.1, sd = 0.6, sum(nItems))

# simulate date from multivariate normal distribution
set.seed(830492)
data <- mvtnorm::rmvnorm(N, mean = M, sigma = sigmaMat)
colnames(data) <- c(paste0("X", seq_len(nItems[1])))
data <- data.frame(data)

# get random outcome (purpose is only looking inside matrices!)
data$y <- rnorm(n = dim(data)[1], m = 0, sd = 1)

# estimate SEM according to model idea to look at matrices in SEM
# -> necessary to do all latent formulation of y!
#     that is, set factor loading to one and residual variance to zero
library(lavaan)

m0 <- '# measurement model
       f =~ 0.8*X1 + 0.6*X2 + 0.6*X3
       n1 =~ 0.2*X1 + 0.2*X2
       
       # uncorrelated factor and nuances
       f ~~ 0*n1
      
       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
    
       # structure model
       fY ~ f + n1'

fit <- sem(m0, data = data)
summary(fit, fit.measures = T, standardized = T)

# this is the matrix of path coefficients
# in columns: {f, n1, fY} -> factors in order of appearance within model syntax
# in rows: {f, n1, fY} -> factors in order of appearance within model syntax
fit@Model@GLIST[["beta"]]

# this is the matrix of factor loadings
# in columns: {f, n1, fY} -> factors in order of appearance within model syntax
# in rows: {X1, X2, X3, y} -> item in order of appearance within model syntax
fit@Model@GLIST[["lambda"]]

# this is the matrix of factor correlations
# in columns: {f, n1, fY} -> factors in order of appearance within model syntax
# in rows: {f, n1, fY} -> factors in order of appearance within model syntax
fit@Model@GLIST[["psi"]]

# model-implied mean matrix
fit@h1[["implied"]][["mean"]][[1]]

# # model-implied variance-covariance matrix
# # -> what is the difference between these to matrices
# #### same
# fit@Fit@Sigma.hat[[1]]
# fit@implied[["cov"]][[1]]
# 
# #### different
# fit@h1[["implied"]][["cov"]][[1]]
# 
################################################################################
# data simulation from SEM
# to do: switch to standardized factor loadings
# in order to simulate associations between factor(s) and outcome
N <- 10000
nOutcome <- 1
nItems <- 6
lambdaF1 <- 1
lambdaF2 <- 0.8
lambdaF3 <- 0.2

lambdaN <- 0.2 # {0.2, 0.5, 0.8}
# fake Belohnung für unreliable Prädiktoren (siehe Mulaik?)

# matrix of factor loadings
# factors and items in order of appearance within model syntax
# in columns: {f, n1, fY}; in rows: {X1, X2, X3, y} 

# hier: lambda factor = {1, 0.8, 0.2, 0.2, 0.2, 0.2}, lambda nuance = {0.2, 0.2, 0.2, 0, 0, 0}
#       sind factor und nuance sich zu ähnlich! dadurch factor = items
# lambdaSEM <- cbind(factor = c(lambdaF1, lambdaF2, lambdaF3, lambdaF3, lambdaF3, lambdaF3, rep(0, nOutcome)),

# hier: lambda factor = {1, 0.8, 0.2, 1, 0.8, 0.2}, lambda nuance = {0.2, 0.2, 0.2, 0, 0, 0}
#       bilden factor und nuance unterschiedliche Dinge ab, sodass factor << items (für lambdaN = 0.2)
#       -> Wie soll diese Struktur vorhanden sein, wenn immer nur 1 Item für nuance factor vorhanden ist?
lambdaSEM <- cbind(factor = c(lambdaF1, lambdaF2, lambdaF3, lambdaF1, lambdaF2, lambdaF3, rep(0, nOutcome)),
                   # nuance1 = c(lambdaN, lambdaN, lambdaN, 0, 0, 0, rep(0, nOutcome)),
                   # nuance1 = c(lambdaN, 0.5, 0.8, 0, 0, 0, rep(0, nOutcome)), # this increases advantage of items over factor
                   nuance1 = c(0.8, 0.5, 0.2, 0, 0, 0, rep(0, nOutcome)),
                   fY = c(rep(0, nItems), rep(1, nOutcome)))
nFactors <- ncol(lambdaSEM)

# factors are not correlated 
#   nuance information is independent from facet information (this is the dgm, i.e., 
#   extra information in nuances that has nothing to do with what is already explained by facet)
# variables from structure model need to be included in factor correlation matrix 
# for the SEM simulation 
corSEM <- diag(nFactors) * 0.5
# one additional nuance factor
corSEM[upper.tri(corSEM, diag = F)] <- c(0,
                                         0, 0) 
corSEM <- corSEM + t(corSEM)

# path coefficient matrix
# ! factors in order of appearance within model syntax
# in columns: {f, n1, fY}; in rows: {f, n1, fY} 
betaFY <- 0.6 # path coefficient factor -- outcome
betaNY <- 0.8 # path coefficient nuance -- outcome {0.8, 0.4}
# je kleiner betaNY, desto kleiner die Unterschiede zwischen Analysemodellen {dgm, factor, items, ...}
betaMat <- matrix(c(0, 0, 0,
                    0, 0, 0,
                    betaFY, betaNY, 0), ncol = nFactors, byrow = T)

# error
# {0.2, 0.4}
thetaSEM <- diag(nItems + nOutcome) * 0.2

# calculate model implied covariance matrix
# for simulating from SEM
# # test inverse
# solve(diag(3) - betaMat)
# MASS::ginv(diag(3) - betaMat)

sigmaSEM <- lambdaSEM %*% solve(diag(nFactors) - betaMat) %*% corSEM %*% 
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaSEM) + thetaSEM


# set mean values for variables
set.seed(289301)
M <- c(rnorm(mean = 3.1, sd = 0.6, sum(nItems)), rnorm(mean = 10, sd = 3, nOutcome))

# simulate date from multivariate normal distribution
set.seed(830492)
data <- mvtnorm::rmvnorm(N, mean = M, sigma = sigmaSEM)
colnames(data) <- c(paste0("X", seq_len(nItems[1])), "y")
data <- data.frame(data)

# warum muss ich eine Korrelation zwischen Faktoren oder eine Ladung festsetzen?
#     à la bi-factor model(?)
#     assumption of uncorrelated factor and nuances fine from a substantive point
#       of view; even if facet (and therefore domain) variance is partialed/regressed out 
#       from items/nuances, predictive power does not change significantly which 
#       indicates uncorrelated facet and nuance factor
# warum habe ich 6 (bei 5 Items) bzw. 11 (bei 6 Items) df aber ohne das fixieren 
#     ist das Modell identifiziert?
m0 <- '# measurement model
       f =~ X1 + X2 + X3 + X4 + X5 + X6
       n1 =~ X1 + X2 + X3

       # uncorrelated factor and nuances
       f ~~ 0*n1
      
       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ f + n1'

# fit <- sem(mTrue, data = dataList[["lambdaN0.2_error0.2_betaNY0.6_betaFY0.6"]])
fit <- sem(m0, data = data)
summary(fit, fit.measures = T, standardized = T)
(trueR2 <- inspect(fit, "rsquare")["fY"])
################################################################################
# fit different models and compare their Rsquared
## only factor as predictor and ignore nuance level explanation of the outcome
##    -> this should be the worst case...
m1 <- '# measurement model
       f =~ X1 + X2 + X3 + X4 + X5 + X6
      
       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ f'

fit1 <- sem(m1, data = data)
#summary(fit1, fit.measures = T, standardized = T)

inspect(fit1, "rsquare")["fY"]
(deltaR2.ignore <- inspect(fit1, "rsquare")["fY"] - trueR2)

# with lambdaN = 0.2: (R^2 decreases as a result of ignoring nuance information)
#   R² = 0.2686896; deltaR2.ignore = -0.187608
# with lambdaN = 0.8: (R^2 only shows minor change in R²)
#   R² = 0.4482701; deltaR2.ignore = -0.01008881 

## use items with nuance level information as additional manifest predictors of the outcome
## == combined analysis in literature (factor + items)
m2 <- '# measurement model
       f =~ X1 + X2 + X3 + X4 + X5 + X6

       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ f + X1 + X2 + X3'

fit2 <- sem(m2, data = data)
#summary(fit2, fit.measures = T, standardized = T)
inspect(fit2, "rsquare")["fY"]
(deltaR2.addNuances <- inspect(fit2, "rsquare")["fY"] - trueR2)

# with lambdaN = 0.2: (R^2 decreases as a result of adding nuance information)
#   R² = 0.2825819; deltaR2.addNuances = -0.1737157 
# with lambdaN = 0.8: (R^2 only shows minor change in R²)
#   R² = 0.4313762; deltaR2.ignore = -0.02698269 

### add only one additional manifest nuance level predictor 
### -> nuance level predictors are highly correlated and therefore can damage performance
###     again    
m2a <- '# measurement model
       f =~ X1 + X2 + X3 + X4 + X5 + X6
      
       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ f + X1'

fit2a <- sem(m2a, data = data)
#summary(fit2, fit.measures = T, standardized = T)
inspect(fit2a, "rsquare")["fY"]
(deltaR2.addOneNuance <- inspect(fit2a, "rsquare")["fY"] - trueR2)

# with lambdaN = 0.2: (R^2 decreases as a result of adding nuance information)
#   R² = 0.2874407 ; deltaR2.addNuances = -0.1688569  
# with lambdaN = 0.8: (R^2 only shows minor change in R²)
#   R² = 0.4491565 ; deltaR2.ignore = -0.009202325  

## only use items to predict outcome (basically regression type analysis) 
## == nuance level analysis in literature (only items) without adding factor
m3 <- '# all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ X1 + X2 + X3 + X4 + X5 + X6'

fit3 <- sem(m3, data = data)
#summary(fit3, fit.measures = T, standardized = T)
inspect(fit3, "rsquare")["fY"]
(deltaR2.items <- inspect(fit3, "rsquare")["fY"] - trueR2)

# with lambdaN = 0.2: (R^2 decreases as a result of using only nuance information)
#   R² = 0.2667204; deltaR2.addNuances = -0.1895772 
# with lambdaN = 0.8: (R^2 only shows minor change in R²)
#   R² = 0.4311547 ; deltaR2.ignore = -0.0272042  

# # Rsquared in this SEM is equivalent to multiple linear regression 
# fit3.lm <- lm(y ~ X1 + X2 + X3 + X4 + X5 + X6, data = data)
# summary(fit3.lm)$r.squared

data$sumScore <- (data$X1 + data$X2 + data$X3 + data$X4 + data$X5 + data$X6)/6 
fit4.lm <- lm(y ~ sumScore, data = data)
summary(fit4.lm)$r.squared
(deltaR2.sumScore <- summary(fit4.lm)$r.squared - trueR2)

fit5.lm <- lm(y ~ sumScore + X1 + X2 + X3 + X4 + X5 + X6, data = data)
summary(fit5.lm)$r.squared
(deltaR2.sumScoreItems <- summary(fit5.lm)$r.squared - trueR2)


predictivePerformance <- c("dgm" = inspect(fit, "rsquare")[["fY"]],  # dgm
                           "factor" = inspect(fit1, "rsquare")[["fY"]], # only factor
                           "items" = inspect(fit3, "rsquare")[["fY"]], # only items
                           "factor+items" = inspect(fit2, "rsquare")[["fY"]], # factor + items 
                           "sumScore" = summary(fit4.lm)$r.squared, # only sum Score
                           "sumScore+items" = summary(fit5.lm)$r.squared) # sum score + items
predictivePerformance

# Mottus et al. formulieren in ihren inhaltlichen Beschreibungen dieses dgm
# nuancen bilden zusätzlichen Faktor ab (trait), der Outcome beeinflusst
# aber: auf Basis dieses dgms ist R^2 mit Items als Prädiktoren nicht besser als 
#       nur Faktor als Prädiktor
#   -> Suche das datengenerierende Modell, für das Performanz in R² mit Einzelitems >>> Performanz 
#     mit Faktor als Prädiktor wird?
# Odd: mit lambdaN = 0.8 & betaNY = 0.4 factor + items ist plötzlich besser als "trueR2" -> Why?!
# erweitere Mini-Simulation mit diesem dgm: Nuancenfaktor simulieren
# systamtisch manipulieren: 
#     - Faktorladungen in Nuancenfaktor größer machen
#     - theta kleiner machen (Residualvarianzen in Items)
#     - Effekt von Nuancenfaktor/Facettenfaktor auf Outcome variieren

# Ladungsstruktur: hohe Ladungen factor != nuance; wie wenn nur 1 Item?
# lambdaN = {0.2, 0.5, 0.8}
# error = {0.2, 0.4}
# betaNY = {0.4, 0.6, 0.8} # 0.6 == betaFY

# anderes dgm testen: use residuals of each item as predictor for outcome
# andere Modellierung testen: Summenscore/Mittelwertscore als Faktor/Facettenscore 
#     (Prädiktor) für Outcome benutzen (= manifester Prädiktor)
#     - Skrondal & Laaker?  

################################################################################
# dgm: residuals of each item as predictor for outcome
################################################################################


################################################################################
# simulate MI violations
################################################################################
# simulated measurement invariance
nOutcome <- 1
nItems <- 6
nGroups <- 2 # {2, 4}
N <- 10000
N <- nGroups * N

# simulate two factors with item that either loads on one or the other factor

# matrix of factor loadings
# factors and items in order of appearance within model syntax
# in columns: {f, n1, fY}; in rows: {X1, X2, X3, y} 
##### configural
lambdaConfigMI.1 <- cbind(factor = c(1, 0.8, 0.2, 1, 0.8, 0.2, rep(0, nOutcome)),
                          fY = c(rep(0, nItems), rep(1, nOutcome)))
lambdaConfigMI.2 <- cbind(factor = c(0, 0, 0.2, 1, 0, 0.2, rep(0, nOutcome)),
                          fY = c(rep(0, nItems), rep(1, nOutcome)))
# lambdaConfigMI.3 <- cbind(factor = c(0, 0, 0.2, 1, 0.4, 0, rep(0, nOutcome)),
#                           fY = c(rep(0, nItems), rep(1, nOutcome)))
# lambdaConfigMI.4 <- cbind(factor = c(0, 0, 0, 1, 0.4, 0.2, rep(0, nOutcome)),
#                           fY = c(rep(0, nItems), rep(1, nOutcome)))

# höhere Ladungen in zweiter Gruppe! sonst hat Faktor weniger Varianz -> R^2 sinkt

##### metric (weak)
# more extreme: change direction of effect between groups
lambdaConfigMI.1 <- cbind(factor = c(1, 0.8, 0.2, 1, 0.8, 0.2, rep(0, nOutcome)),
                          fY = c(rep(0, nItems), rep(1, nOutcome)))
lambdaConfigMI.2 <- cbind(factor = c(0.2, -0.2, -0.4, 1, 0.2, -0.4, rep(0, nOutcome)),
                          fY = c(rep(0, nItems), rep(1, nOutcome)))
# lambdaConfigMI.3 <- cbind(factor = c(1, 0.2, 0.8, 1, -0.2, 0.8, rep(0, nOutcome)),
#                           fY = c(rep(0, nItems), rep(1, nOutcome)))
# lambdaConfigMI.4 <- cbind(factor = c(0.2, -0.8, 1, 1, 0.8, 0.2, rep(0, nOutcome)),
#                           fY = c(rep(0, nItems), rep(1, nOutcome)))
nFactors <- ncol(lambdaConfigMI.1)

# factors are not correlated 
#   nuance information is independent from facet information (this is the dgm, i.e., 
#   extra information in nuances that has nothing to do with what is already explained by facet)
# variables from structure model need to be included in factor correlation matrix 
# for the SEM simulation 
corConfigMI <- diag(nFactors) * 0.5
# one additional nuance factor
corConfigMI[upper.tri(corConfigMI, diag = F)] <- c(0) 
corConfigMI <- corConfigMI + t(corConfigMI)

# path coefficient matrix
# ! factors in order of appearance within model syntax
# in columns: {f, n1, fY}; in rows: {f, n1, fY} 
betaFY <- 0.8 # path coefficient factor -- outcome {0.4, 0.8}
betaMat <- matrix(c(0, 0,
                    betaFY, 0), ncol = nFactors, byrow = T)

# error
thetaMI <- diag(nItems + nOutcome) * 0.2

# calculate model implied covariance matrix
# for simulating from SEM
# # test inverse
# solve(diag(3) - betaMat)
# MASS::ginv(diag(3) - betaMat)

sigmaConfigMI.1 <- lambdaConfigMI.1 %*% solve(diag(nFactors) - betaMat) %*% corConfigMI %*% 
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaConfigMI.1) + thetaMI
sigmaConfigMI.2 <- lambdaConfigMI.2 %*% solve(diag(nFactors) - betaMat) %*% corConfigMI %*% 
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaConfigMI.2) + thetaMI
# sigmaConfigMI.3 <- lambdaConfigMI.3 %*% solve(diag(nFactors) - betaMat) %*% corConfigMI %*%
#   t(solve(diag(nFactors) - betaMat)) %*% t(lambdaConfigMI.3) + thetaMI
# sigmaConfigMI.4 <- lambdaConfigMI.4 %*% solve(diag(nFactors) - betaMat) %*% corConfigMI %*%
#   t(solve(diag(nFactors) - betaMat)) %*% t(lambdaConfigMI.4) + thetaMI

# set mean values for variables
set.seed(289301)
M <- c(rnorm(mean = 3.1, sd = 0.6, sum(nItems)), rnorm(mean = 10, sd = 3, nOutcome))

# simulate date from multivariate normal distribution
set.seed(830492)
dataConfigMI.1 <- mvtnorm::rmvnorm(N/nGroups, mean = M, sigma = sigmaConfigMI.1)
dataConfigMI.2 <- mvtnorm::rmvnorm(N/nGroups, mean = M, sigma = sigmaConfigMI.2)
# dataConfigMI.3 <- mvtnorm::rmvnorm(N/nGroups, mean = M, sigma = sigmaConfigMI.3)
# dataConfigMI.4 <- mvtnorm::rmvnorm(N/nGroups, mean = M, sigma = sigmaConfigMI.4)
colnames(dataConfigMI.1) <- c(paste0("X", seq_len(nItems[1])), "y")
colnames(dataConfigMI.2) <- c(paste0("X", seq_len(nItems[1])), "y")
# colnames(dataConfigMI.3) <- c(paste0("X", seq_len(nItems[1])), "y")
# colnames(dataConfigMI.4) <- c(paste0("X", seq_len(nItems[1])), "y")
dataConfigMI.1 <- data.frame(dataConfigMI.1)
dataConfigMI.2 <- data.frame(dataConfigMI.2)
# dataConfigMI.3 <- data.frame(dataConfigMI.3)
# dataConfigMI.4 <- data.frame(dataConfigMI.4)

dataConfigMI.1$groupe <- 0
dataConfigMI.2$groupe <- 1
# dataConfigMI.3$groupe <- 2
# dataConfigMI.4$groupe <- 3
dataConfigMI <- rbind(dataConfigMI.1, dataConfigMI.2)
# dataConfigMI <- rbind(dataConfigMI.1, dataConfigMI.2, dataConfigMI.3, dataConfigMI.4)
str(dataConfigMI)
dataConfigMI$groupe <- as.factor(dataConfigMI$groupe)

# "sanity check" 
model <- 'f =~ NA*X1 + X2 + X3 + 1*X4 + X5 + X6'

# without measurement invariance violations
fit0 <- cfa(model, data = dataConfigMI, meanstructure = T)

# configural invariance
fit1 <- cfa(model, data = dataConfigMI, group = "groupe")

# weak invariance
fit2 <- cfa(model, data = dataConfigMI, group = "groupe",
            group.equal = "loadings")

# strong invariance
fit3 <- cfa(model, data = dataConfigMI, group = "groupe",
            group.equal = c("intercepts", "loadings"))

# model comparison tests
lavTestLRT(fit0, fit1, fit2, fit3)

# oracle R^2 is multiple group SEM
mMI.mg <- '# measurement model
          f =~ NA*X1 + X2 + X3 + 1*X4 + X5 + X6
          
          # all latent representation of manifest outcome
          fY =~ 1 * y
          y ~~ 0*y
          
          # structure model
          fY ~ f'
fitMI.mg <- sem(mMI.mg, data = dataConfigMI, group = "groupe",
                group.equal = "loadings")
summary(fitMI.mg, fit.measures = T, standardized = T)
inspect(fitMI.mg, "rsquare")

# was ist das "Gesamt-R²" in einer multiple group Analyse?
# researcher is ignorant of multi group structure
#   -> R² for comparing model performances is single group model

# fit model to data
## only factor as predictor and ignore explanation of the outcome by groupe
mMI.ignore <- '# measurement model
               f =~ NA*X1 + X2 + X3 + 1*X4 + X5 + X6
                     
               # all latent representation of manifest outcome
               fY =~ 1 * y
               y ~~ 0*y
               
               # structure model
               fY ~ f'

fitMI.ignore <- sem(mMI.ignore, data = dataConfigMI)
#summary(fitMI.ignore, fit.measure = T, standardized = T)
inspect(fitMI.ignore, "rsquare")["fY"]

mMI.itemMI <- '# measurement model
               f =~ NA*X1 + X2 + X3 + 1*X4 + X5 + X6
                     
               # all latent representation of manifest outcome
               fY =~ 1 * y
               y ~~ 0*y
               
               # structure model
               fY ~ f + X1 + X2 + X3'

fitMI.itemMI <- sem(mMI.itemMI, data = dataConfigMI)
inspect(fitMI.itemMI, "rsquare")["fY"]

mMI.item12 <- '# measurement model
             f =~ NA*X1 + X2 + X3 + 1*X4 + X5 + X6
                     
             # all latent representation of manifest outcome
             fY =~ 1 * y
             y ~~ 0*y
               
             # structure model
             fY ~ f + X1 + X2'

fitMI.item12 <- sem(mMI.item12, data = dataConfigMI)
inspect(fitMI.item12, "rsquare")["fY"]

mMI.item1 <- '# measurement model
             f =~ X1 + X2 + X3 + X4 + X5 + X6
                     
             # all latent representation of manifest outcome
             fY =~ 1 * y
             y ~~ 0*y
               
             # structure model
             fY ~ f + X2'

fitMI.item1 <- sem(mMI.item1, data = dataConfigMI)
inspect(fitMI.item1, "rsquare")["fY"]

mMI.item2 <- '# measurement model
             f =~ X1 + X2 + X3 + X4 + X5 + X6
                     
             # all latent representation of manifest outcome
             fY =~ 1 * y
             y ~~ 0*y
               
             # structure model
             fY ~ f + X2'

fitMI.item2 <- sem(mMI.item2, data = dataConfigMI)
inspect(fitMI.item2, "rsquare")["fY"]

## only use items to predict outcome without groupe (basically regression type analysis) 
mMI.itemsOnly <- '
       # all latent representation of manifest outcome
       fY =~ 1 * y
       y ~~ 0*y
       
       # structure model
       fY ~ X1 + X2 + X3 + X4 + X5 + X6'

fitMI.itemsOnly <- sem(mMI.itemsOnly, data = dataConfigMI)
inspect(fitMI.itemsOnly, "rsquare")["fY"]

predictivePerformance <- c("factor" = inspect(fitMI.ignore, "rsquare")[["fY"]], # only factor
                           "items" = inspect(fitMI.itemsOnly, "rsquare")[["fY"]], # only items
                           "factor+items" = inspect(fitMI.itemMI, "rsquare")[["fY"]], # factor + all items
                           "factor+item12" = inspect(fitMI.item12, "rsquare")[["fY"]], # factor + items 1 & 2
                           "factor+item1" = inspect(fitMI.item1, "rsquare")[["fY"]],  # factor + item 1 
                           "factor+item2" = inspect(fitMI.item2, "rsquare")[["fY"]])  # factor + item 2 

predictivePerformance

# ## only use items to predict outcome with groupe (basically regression type analysis) 
# mMI.items <- '# all latent representation of manifest outcome
#               fY =~ 1 * y
#               y ~~ 0*y
#                      
#               # structure model
#               fY ~ X1 + X2 + X3 + X4 + X5 + X6 + groupe'
# 
# fitMI.items <- sem(mMI.items, data = dataConfigMI)
# inspect(fitMI.items, "rsquare")["fY"]

# ## use items with groupe information as additional manifest predictor of the outcome
# mMI.groupe <- '# measurement model
#                f =~ X1 + X2 + X3 + X4 + X5 + X6
#                      
#                # all latent representation of manifest outcome
#                fY =~ 1 * y
#                y ~~ 0*y
#                
#                # structure model
#                fY ~ f + groupe'
# 
# fitMI.groupe <- sem(mMI.groupe, data = dataConfigMI)
# inspect(fitMI.groupe, "rsquare")["fY"]

# Reliabilitätseffekt? 
# je größer Reliabilität der Items, desto näher sind Faktor  

fitLM <- lm(y ~ X1 + X2 + X3 + X4 + X5 + X6, data = dataConfigMI)
summary(fitLM)$"r.squared"

fitLM <- lm(y ~ X1 + X2 + X3 + X4 + X5 + X6 + groupe +
     X1:groupe + X2:groupe + X3:groupe, data = dataConfigMI)
summary(fitLM)$"r.squared"

mMI.nuance <- '# measurement model
             f =~ X1 + X2 + X3 + X4 + X5 + X6
             n1 =~ X1 + X2 + X3
            
             f ~~ 0*n1

             # all latent representation of manifest outcome
             fY =~ 1 * y
             y ~~ 0*y
               
             # structure model
             fY ~ f + n1'

fitMI.nuance <- sem(mMI.nuance, data = dataConfigMI)
inspect(fitMI.nuance, "rsquare")["fY"]


################################################################################
library(ggcorrplot)

# simulate data without additional nuance factor or measurement invariance violation
N <- 10000
nOutcome <- 1
nItems <- 6
lambdaF1 <- 1
lambdaF2 <- 0.8
lambdaF3 <- 0.2

# matrix of factor loadings
lambdaZero <- cbind(factor = c(lambdaF1, lambdaF2, lambdaF3, lambdaF1, lambdaF2, lambdaF3, rep(0, nOutcome)),
                   fY = c(rep(0, nItems), rep(1, nOutcome)))
nFactors <- ncol(lambdaZero)

# factors are not correlated 
#   nuance information is independent from facet information (this is the dgm, i.e., 
#   extra information in nuances that has nothing to do with what is already explained by facet)
# variables from structure model need to be included in factor correlation matrix 
# for the SEM simulation 
corZero <- diag(nFactors) * 0.5
# one additional nuance factor
corZero[upper.tri(corZero, diag = F)] <- c(0) 
corZero <- corZero + t(corZero)

# path coefficient matrix
betaFY <- 0.6 # path coefficient factor -- outcome
betaMat <- matrix(c(0, 0,
                    betaFY, 0), ncol = nFactors, byrow = T)

# error
# {0.2, 0.4}
thetaZero <- diag(nItems + nOutcome) * 0.2

# calculate model implied covariance matrix
# for simulating from SEM
# # test inverse
# solve(diag(3) - betaMat)
# MASS::ginv(diag(3) - betaMat)

sigmaZero <- lambdaZero %*% solve(diag(nFactors) - betaMat) %*% corZero %*% 
  t(solve(diag(nFactors) - betaMat)) %*% t(lambdaZero) + thetaZero


# set mean values for variables
set.seed(289301)
M <- c(rnorm(mean = 3.1, sd = 0.6, sum(nItems)), rnorm(mean = 10, sd = 3, nOutcome))

# simulate date from multivariate normal distribution
set.seed(830492)
dataZero <- mvtnorm::rmvnorm(N, mean = M, sigma = sigmaZero)
colnames(dataZero) <- c(paste0("X", seq_len(nItems[1])), "y")
dataZero <- data.frame(dataZero)

pZero <- ggcorrplot(cor(dataZero[, 1:6]))
pNuance <- ggcorrplot(cor(data[, 1:6]))
pMI <- ggcorrplot(cor(dataConfigMI[, 1:6]))
pMI.1 <- ggcorrplot(cor(dataConfigMI[which(dataConfigMI$groupe == 0), 1:6]))
pMI.2 <- ggcorrplot(cor(dataConfigMI[which(dataConfigMI$groupe == 1), 1:6]))

library(patchwork)
pZero | pNuance | pMI | pMI.1 | pMI.2
