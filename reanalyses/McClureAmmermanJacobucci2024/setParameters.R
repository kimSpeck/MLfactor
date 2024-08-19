params <- list()

params$nRep <- 500
params$M = 1 # number of composite variables
params$holdout = 1000 # hold out sample size

params$N <- c(50,100,250,500,1000) #this is training sample size; test/holdout sample is set to 1000 inside function
params$J <- c(9,21) # number of items in the composite 
if (!(all(params$J %% 3 == 0))) { # assumption in setBetaS; see condition 3 in item-specific effects
  stop("every item number has to be divisible by 3!")
}
params$betaT <- c(0, 0.25, 0.5, 0.75) # true score effect
params$betaS <- c(0, 0.1, 0.3, 0.5) # no effect, small effect, moderate effect large effect

##### item specific effects conditions #####
# see table 1 in paper
params$betaS.cond <- 1:5

##### variance components conditions #####
# see table 1 in paper
# use actual values by indexing matrix sigma
params$measVar.cond <- 1:5

# T, S, and e variance in columns; conditions to index with measVar.cond in rows
# half of non T variance is signal
#     variance for specific effect is 50% signal (predictive for outcome) and 50% noise
# most non T variance (90%) is signal
#     variance for specific effect is 90% signal (predictive for outcome) and 10% noise
# equal impact of T, S, and e on X
#     variance for specific effect is 33% signal (predictive for outcome) and 33% noise
# variance for specific effect is 10% signal (predictive for outcome) and 90% noise
# more S than T
params$sigmaS <- c(0.25, 0.45, 0.33, 0.05, 0.5)                  # item specific variance
params$sigma <- cbind(sigmaT = c(0.5, 0.5, 0.33, 0.5, 0.25),     # true score variance
                      sigmaS = params$sigmaS,                    # item specific variance
                      sigmaE = c(0.25, 0.05, 0.33, 0.45, 0.25))  # error variance 

# manipulation measurement invariance violations
params$nGroups = 2  # number of groups (MI due to categorical variable)
params$miItems.cond = c(1, 2, 3)  # number of items with deviation in factor loadings
# configural MI by using 0 as a weight
# metric MI by using -0.8 as a weight
# scalar MI by using unit weight (= identical weights) but instead change mean of specific item effect
params$miWeight = c(0, -0.8, 1) # weight of measurement variant items

params$mS <- 0
params$mMI <- c(1, 2, 5, 10) 

#irreducible error variance (= error in the outcome)
# for the "formative" model; R² is .9 and higher
params$sigmaErrY = 1 
# for the reflective model; R² is .9, .65, .5, .33
params$sigmaErrY_ref <- c(0.01, 0.05, 0.1, 0.2) 

# to simulate from kind of reflective measurement model (R^2 not dependent on number of items)
# with the error variance from Jacobucci on the outcome (sigma^2_e = 1): 
#   ... R^2 is around 0.1 without specific item effect
#   ... R^2 is arount 0.5 with specific item effect
# now, the specific item effect is really dominating the prediction
# only simulating specific item effects on the outcome causes superiority of items 
#     over the composite; no measurement invariance violation is able to cause item 
#     superiority in terms of prediction
# therefore, try increasing R² to the order of magnitude it had beforehand and compare 
#     results again.
# --> by reducing the error in Y, we can increase R²
# # explained variance with varying error variance in the outcome for reflective model
# sigmaErrY <- 0.01 # 0.899537 1.002019
# sigmaErrY <- 0.2 # 0.3389238 0.8858312
# sigmaErrY <- 0.1 # 0.4995846 0.9484514
# sigmaErrY <- 0.05 # 0.6581100 0.9803542
################################################################################
# get all simulation condition (as in original paper)
# N.train x J x betaT x item specific effect x variance components 
# 1000 conditions
# refSub <- rbind(cbind(sigmaErrY = params$sigmaErrY, reflective = FALSE), 
#                 expand.grid(sigmaErrY = params$sigmaErrY_ref, reflective = TRUE))

grid_noMIvio_form = expand.grid(N.train = params$N,
                                J = params$J,
                                betaT = params$betaT,
                                betaS.cond = params$betaS.cond,
                                measVar.cond = params$measVar.cond,
                                # these are all parameters to turn MI violations off
                                MI = TRUE,
                                nGroups = 1,
                                miItems.cond = 0, 
                                miWeight = NA,
                                mSpecific = params$mS, 
                                sigmaErrY = params$sigmaErrY, 
                                # "formative" model, therefore only sigmaErrY = 1
                                reflective = FALSE) 

grid_noMIvio_ref <- expand.grid(N.train = params$N,
                                J = params$J,
                                betaT = params$betaT,
                                betaS.cond = params$betaS.cond,
                                measVar.cond = params$measVar.cond,
                                # these are all parameters to turn MI violations off
                                MI = TRUE,
                                nGroups = 1,
                                miItems.cond = 0, 
                                miWeight = NA,
                                mSpecific = params$mS, 
                                # reflectiv model, therefore sigmaErrY varies to manipulate R²
                                sigmaErrY = params$sigmaErrY_ref, 
                                reflective = TRUE)

# add MI simulation conditions
grid_MIvio_form <- expand.grid(N.train = params$N,
                               J = params$J,
                               betaT = params$betaT,
                               betaS.cond = params$betaS.cond,
                               measVar.cond = params$measVar.cond,
                               # these are all parameters to turn MI violations on
                               MI = FALSE,
                               nGroups = params$nGroups, # two groups
                               miItems.cond = params$miItems.cond, 
                               miWeight = params$miWeight, # configural, metric, scalar
                               mSpecific = params$mMI, # different degrees of scalar MI vio
                               # "formative" model, therefore sigmaErrY varies to manipulate R²
                               sigmaErrY = params$sigmaErrY,
                               reflective = FALSE)

grid_MIvio_ref <- expand.grid(N.train = params$N,
                              J = params$J,
                              betaT = params$betaT,
                              betaS.cond = params$betaS.cond,
                              measVar.cond = params$measVar.cond,
                              # these are all parameters to turn MI violations on
                              MI = FALSE,
                              nGroups = params$nGroups, # two groups
                              miItems.cond = params$miItems.cond, 
                              miWeight = params$miWeight, # configural, metric, scalar
                              mSpecific = params$mMI, # different degrees of scalar MI vio
                              # reflectiv model, therefore sigmaErrY varies to manipulate R²
                              sigmaErrY = params$sigmaErrY_ref, 
                              reflective = TRUE)

params$gridFull <- rbind(grid_noMIvio_form, grid_noMIvio_ref, 
                         grid_MIvio_form, grid_MIvio_ref)
dim(params$gridFull) # 185 000 conditions in total; we will only be abel to run a subset!!!

# seed for every condition
set.seed(784359) # seed for reproducibility of seeds
params$gridFull$conditionSeed <- sample(843786, dim(params$gridFull)[1], replace = FALSE)
# # check (has to be true)
# length(unique(seed)) == dim(gridFull)[1]

rm(grid_noMIvio_form, grid_noMIvio_ref, 
   grid_MIvio_form, grid_MIvio_ref)
