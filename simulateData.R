# generate grid for all simulation conditions
# simulation parameters:
# ... dgp: simulate from SEM with/without MI violation (vs. simulate nuance factor)
# ... type of measurement invariance violation {none, config, metric, scalar}
#     (simulate none from two groups without difference in parameters between groups)
#     (always on the last items)
# ... correlation between grouping variable and outcome {0, 0.2, 0.5}
# ... (with and) without specific item effect (always on first item)

# Kontrollbedingung: Messinvarianz-Annahmen sind erfüllt + spezifischer Item Effekt ist dabei
#   -> R^2 hier ungefähr konstant halten

# To Do: 
# Korrelationen zw. Gruppe und Outcome in Cohens d umrechnen
# standardisierte Faktoreladungen anschauen und Ladungen des True Scores evtl noch runter setzen
# Ladungen in metrischer MI Verletzung ändern
# Reliabilität manipulieren (Reliabilität wirkt auf Regression über Attenuationsbias)
#   -> recherchiere Attenuationsbias: was sind Konsequenzen von Attenuationsbias auf Regressionsmodell
#   -> Power bei kleinem N? Summenscore besser als Faktorenmodell?

# Ziele: 
# - alternatives datengenerierendes Modell vorschlagen
# - Differenz im composite und item basierten Modell 

# within specific dgp
# ... number of observations N = {50, 100, 250, 500, 1000}
# ... number of items with MI violation {1, J/3, J/2} 

# ... explained variance/simulated R² = {0.1, 0.3, 0.5, 0.8} 
#   -> this will either be manipulated by manipulating beta coefficients or error variance of the outcome
#   -> how does correlation between group and outcome manifest in simualted R^2 (konstant halten oder variieren lassen)

# modell-implizirte Kovarianzmatrix zeigen!
# in Intro an Literatur zur Äquivalenz unterschiedlicher modell-implizierter Kov-Matrizen für kausale Modelle anknüpfen

# Varianz der manifesten Variablen konstant halten und Manipulation von R² über Koeffizienten manipulieren
# -> standardisierte Koeffizienten sind immer dasselbe! 

# iterate across grid of simulation conditions

# simulate data by simulating the following matrices step by step

# simulate lambda
# -> lambda matrix differs between groups for configural and metric MI violation 
#     (manipulation via miWeights)
#     (control this by miItems = 0 & miWeight = list(c(NA)) for none;
#       by miItems = miItems & miWeight = params$dgp$trueScoreLoad for scalar;
#       otherwise miItems = miItems & miWeight = miWeight for config (0) or metric)
lambdaSEM <- defineLambda(nGroups, nItems = J, miItems = 0, miWeight = list(c(NA)), 
                          params$dgp$nOutcome)
nFactors <- ncol(lambdaSEM)
nVariables <- nrow(lambdaSEM)

# simulate phi
# -> phi matrix changes for scalar MI violation 
#   (variance in specific item effect is 0 in phi but is simulated via mean differences)
#   (control this by setting scalarMI = T/F; no additional parameter since variance 
#     is always set to 0 for scalarMI = T)
covSEM <- definePhi(nFactors, Sigma_T, Sigma_S[1, 1], nItems = J, scalarMI = F, 
                    varNames = colnames(lambdaSEM))

# simulate Beta = path coefficient matrix (effect of factors on the target/criterion)
# -> changes depending on whether we simulate specific item effects or not
# ! factors in order of appearance within model syntax
# in columns: {T, S1,..., SJ, Y}; in rows: {T, S1,..., SJ, Y}
betaMat <- defineBeta(nFactors, betaT, rep(0, J), J, colnames(lambdaSEM))

# simulate theta
# -> theta matrix differs between groups for configural and metric MI violation 
#     (due to reduced loadings from the factor on the item); differing values in simValuesSigmaE!
#     (due to all latent formulation of the outcome, the measurement error for the 
#     outcome is always zero)
thetaSEM <- defineTheta(nVariables, sigma_E, cItemVar = T, 
                        rownames(lambdaSEM), lambdaSEM, covSEM)

# generate the model implied variance-covariance matrix
sigmaSEM <- getSigma(lambdaSEM, betaMat, covSEM, thetaSEM)

# generate the vector of expected values for the items 
# -> M vector differs between groups ...
#     ... for scalar MI violations  
#     ... to control correlation between grouping variable and outcome
M <- defineMeanVec()

# sample data from multivariate normal distribution
# if we have this block matrix structure we only need N/nGroups observations
#   as variables for each group are simulated separately and otherwise we would
#   end up with N * nGroups observations
set.seed(seed)
data_multipleGroup <- mvtnorm::rmvnorm((N.train + params$holdout)/nGroups, mean = M, sigma = sigmaSEM)


# save data
# save stats of data 
#   ... variance of the items, correlations, ... (see checks)



