params <- list()

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

params$sigmaErrY = 1 #irreducible error variance

################################################################################
# get all simulation condition (as in original paper)
# N.train x J x betaT x item specific effect x variance components 
# 1000 conditions
grid = expand.grid(N.train = params$N,
                   J = params$J,
                   betaT = params$betaT,
                   betaS.cond = params$betaS.cond,
                   measVar.cond = params$measVar.cond,
                   MI = TRUE,
                   nGroups = 1,
                   miItems.cond = 0, 
                   miWeight = NA,
                   mSpecific = params$mS)

# add MI simulation conditions
gridMI <- expand.grid(N.train = params$N,
                      J = params$J,
                      betaT = params$betaT,
                      betaS.cond = params$betaS.cond,
                      measVar.cond = params$measVar.cond,
                      MI = FALSE,
                      nGroups = params$nGroups,
                      miItems.cond = params$miItems.cond,
                      miWeight = params$miWeight,
                      mSpecific = params$mMI)

params$gridFull <- rbind(grid, gridMI)

# seed for every condition
set.seed(784359) # seed for reproducibility of seeds
params$gridFull$seed <- sample(843786, dim(params$gridFull)[1], replace = FALSE)
# # check (has to be true)
# length(unique(seed)) == dim(gridFull)[1]

rm(grid, gridMI)
