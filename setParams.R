params <- list()

params$sim$nRep <- 500
params$dgp$holdout = 1000 # hold out sample size
params$dgp$trueScore <- 1 # number of composite variables
# params$dgp$M = 1 # number of composite variables

params$dgp$trueScoreLoad <- 1 # factor loadings of the true score 
params$dgp$specificLoad <- 1 # factor loading specific item latent variable 
params$dgp$nOutcome <- 1 # number of outcome measures
params$dgp$mY_g1 <- 0 # mean of the outcome for group 1 
params$dgp$itemTargetVar <- 1 # constant variance of all manifest variables (items; != outcome)
params$dgp$meanItem <- 0 # mean value of all manifest variables (items; != outcome)

# manipulation measurement invariance violations
# we simulate a balanced design 
params$dgp$nGroups = 2  # number of groups (MI due to binary variable)
params$dgp$varGroup <- 0.25 # due to binary variable and balanced design

# params$mS <- 0

params$dgp$miItems.cond = c(1, 2, 3)  # number of items with deviation in factor loadings
# configural MI by using 0 as a weight
# metric MI by using -0.8 as a weight
# scalar MI by using unit weight (= identical weights) but instead change mean of specific item effect
params$dgp$miWeight = c(0, -0.8, 1) # weight of measurement variant items
params$dgp$rho <- c(0, 0.2, 0.5)
params$dgp$MIvio <- c("none", "config", "metric", "scalar")

params$dgp$N <- c(50,100,250,500,1000) #this is training sample size; test/holdout sample is set to 1000 inside function
params$dgp$J <- c(9,21) # number of items in the composite 
if (!(all(params$dgp$J %% 3 == 0))) { # assumption in setBetaS; see condition 3 in item-specific effects
  stop("every item number has to be divisible by 3!")
}

# generate grid of simulation conditions
# ! for MIvio = "none" miItems.cond are always 0
# ! miWeight is matches to MIvio; not every combination is simulated!
# params$simGrid <- expand.grid()

# generate seed for grid of simulation conditions
# seed for every condition
set.seed(784359) # seed for reproducibility of seeds
params$simGrid$conditionSeed <- sample(843786, dim(params$simGrid)[1], replace = FALSE)
# # check (has to be true)
# length(unique(seed)) == dim(gridFull)[1]