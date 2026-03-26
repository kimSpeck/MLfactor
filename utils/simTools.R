################################################################################
# parameter estimation
################################################################################

################################################################################
# sanity check functions (to control generated data)
################################################################################
getCronbachsAlpha <- function(data) {
  # this function takes simulated data and returns an estimate of Cronbach's alpha
  # ... data should be a numeric matrix (or data.frame)
  # ... data should include one dimensional data (only one factor)
  # ... data should be structured with observations in rows and items in columns
  
  # Number of items
  p <- ncol(data)  
  
  # Sum of item variances 
  itemVars <- apply(data, 2, var, na.rm = TRUE)
  
  # Variance of total (composite) scores
  scoreVar <- var(rowSums(data, na.rm = TRUE), na.rm = TRUE)
  
  # calculate alpha
  alpha <- (p / (p - 1)) * (1 - (sum(itemVars) / scoreVar))
  
  return(alpha)
}

# generate model for sem-function in lavaan
genModel <- function(itemNames, varT, nFactors = 1) {
  
  if (nFactors > 1) {
    stop("Kim, you still need to implement this!")
  }
  
  defineFactor <- paste0("T =~ NA * ", itemNames[1],
                         if (length(itemNames) > 1) {
                           paste0(" + ", paste(itemNames[-1], collapse = " + "))
                         } else {""})
  
  model <- paste0(defineFactor, "\n", 
                  "fY =~ 1 * Y\n", 
                  "Y ~~ 0*Y\n\n",
                  "T ~~ ", varT, "*T\n\n",
                  "fY ~ T")
  return(model)
}

# generate model for sem-function in lavaan
genModel_0varE <- function(itemNames, varT, nFactors = 1) {
  
  if (nFactors > 1) {
    stop("Kim, you still need to implement this!")
  }
  
  defineFactor <- paste0("T =~ NA * ", itemNames[1],
                         if (length(itemNames) > 1) {
                           paste0(" + ", paste(itemNames[-1], collapse = " + "))
                         } else {""})
  
  model <- paste0(defineFactor, "\n", 
                  "fY =~ 1 * Y\n", 
                  "Y ~~ 0*Y\n\n",
                  "T ~~ ", varT, "*T\n\n",
                  paste0(itemNames, " ~~ 0*", itemNames, "\n", collapse = ""),
                  "fY ~ T")
  return(model)
}

################################################################################
# utility functions to simulate data
################################################################################
defineLambda <- function(flT, flS = 1, nGroups, nItems, 
                         miItems, miWeight, nOutcome, nFactors = 1) {
  ###
  # define lambda matrix 
  ## input:
  # flT       - [scalar] factor loadings on the construct (essential tau-equivalent model)
  # flS       - [scalar] factor loadings on the specific item effect (1)
  # nGroups   - [scalar] number of groups with different factor loading patterns to join in block matrix manner
  # nItems    - [scalar] number of items (equal across groups); each group gets this number of items
  # miItems   - [scalar] number of items with measurement invariance violations (identical for every group)
  # miWeight  - [scalar] factor loadings for the non-MI items; groups in list, weights in vector within list 
  # nOutcome  - [scalar] number of outcome variables; !needs to be multiplied by group
  # nFactors  - [scalar] number of factors 
  ## output:
  # lambda     - [matrix] factor loading matrix with factors in columns and items in rows [factors x items]
  ###
  
  # # test it
  # condVec <- which(params$dgp$grid$P == 4 &
  #                 params$dgp$grid$cAlpha == 0.6 &
  #                 params$dgp$grid$varT == 0.5 &
  #                 params$dgp$grid$miItems == 1 &
  #                 params$dgp$grid$MIvio == "metric")
  # iCond <- condVec[1]
  # nItems <- params$dgp$grid$P[iCond]
  # flT <- params$dgp$grid$flT[iCond]
  # flS <- params$dgp$grid$flS[iCond]
  # nGroups <- params$dgp$nGroups
  # miItems <- params$dgp$grid$miItems[iCond]
  # miWeight <- params$dgp$grid$flT[iCond] + params$dgp$grid$b1[iCond]
  # nOutcome <- params$dgp$nOutcome
  # nFactors = 1
  
  # check if arguments are supplied correctly
  # tbd
  if (nFactors != 1) {
    stop("Kim, you still need to implement this!")
  }
  
  # factor loading matrix for reference group 
  lambdaRef <- cbind(TS = rep(flT, nItems),
                     specificEffects = diag(nItems) * 
                       rep(c(params$dgp$flStable, flS), c(nItems-miItems, miItems)))
  
  # multiple group case: factor loading matrices for MI violation group(s)
  #   we always generate multiple group matrices but with identical parameters if no MI violation
  if (nGroups > 1) {
    lambdaGroups <- sapply(seq_len(nGroups - 1), function(iGroup) {
      cbind(TS = rep(c(flT, miWeight), # miWeight is simply flT for "specific" or "scalar"
                     c(nItems - miItems, miItems)),
            specificEffects = diag(nItems) * 
              rep(c(params$dgp$flStable, flS), c(nItems - miItems, miItems)))
    }, simplify = F)
    
    # concatenate factor loading matrices in block matrix style
    lambda <- adiag(lambdaRef, do.call(adiag, lambdaGroups))
  }
  
  # add outcome
  lambda <- adiag(lambda, diag(nOutcome*nGroups))
  
  # name rows and columns 
  colnames(lambda) <- c(paste0(rep(c("T", 
                                     sapply(seq(nItems), 
                                            function(iItem) {paste0("S", iItem)})), times = nGroups), 
                               "_g", rep(seq_len(nGroups), each = nFactors + nItems)), 
                        paste0("Y", seq_len(nOutcome), "_g", 
                               rep(seq_len(nGroups), each = nOutcome)))
  
  rownames(lambda) <- c(paste0(rep(c(sapply(seq(nItems), 
                                            function(iItem) {paste0("X", iItem)})), times = nGroups), 
                               "_g", rep(seq_len(nGroups), each = nItems)), 
                        paste0("Y", seq_len(nOutcome), "_g", 
                               rep(seq_len(nGroups), each = nOutcome)))
  
  return(lambda)
} 

definePhi <- function(nLatentSEM, varT, varStable, varS, nItems, MIvio = F, miItems = 0, 
                      covNuanceFactor = 0, varNames = NULL, nFactors = 1) {
  ###
  # define factor correlation matrix
  ## input:
  # nLatentSEM  - [scalar] number of latent variables in SEM notation
  # varT        - [scalar] variance of the true score composite
  # varStable   - [scalar] variance of the specific item effects (usually identical for all items)
  # varS        - [scalar] variance of the specific item effects (usually identical for all items)
  # nItems      - [scalar] number of items and hence, specific item effect factors
  # MIvio       - [boolean] if MI violation is present variance of the specific item 
  #                         effect is generated by mean differences between groups
  #                         in manifest variables or by moderation of factor loadings
  #                         thus, the specific item variance is 0! 
  # miItems     - [scalar] number of items with MI violation
  # nFactors    - [scalar] number of factors 
  ## output:
  # phi         - [matrix] factor covariance matrix with factors in rows and columns [factor x factor]
  ###
  
  # maintenance!
  if (nFactors != 1) {
    stop("Kim, you still need to implement this!")
  }
  # warning("Kim, implement multiple group variables to simulate correlations with the outcome!")
  
  # # actual input checks
  # if (MIvio && miItems == 0) {
  #   warning("There are no items violating MI assumptions. Is this what you want?!")
  # }
  
  
  # if (!MIvio) { # specific or nuance
  if (MIvio %in% c("specific", paste0("NF_", params$dgp$corNuanceFactor))) {
    phi <- diag(nLatentSEM) * c(rep(c(varT, rep(varS, nItems)), times = params$dgp$nGroups),  
                                # rep(params$dgp$sigmaEy, params$dgp$nOutcome))
                                rep(params$dgp$sigmaEy, params$dgp$nOutcome * params$dgp$nGroups))
                                
  } else { # scalar or metric MI violations 
    phi <- diag(nLatentSEM) * c(
      # true factor with true factor variance, 
      rep(rep(c(varT, varStable, varS), c(nFactors, nItems-miItems, miItems)), times = params$dgp$nGroups),  
      # rep(params$dgp$sigmaEy, params$dgp$nOutcome)) 
      rep(params$dgp$sigmaEy, params$dgp$nOutcome * params$dgp$nGroups)) 
  }
  
  if (!is.null(varNames)){
    colnames(phi) <- rownames(phi) <- varNames
  }
  
  if (MIvio %in% c(paste0("NF_", params$dgp$corNuanceFactor))) {
    miVars_g1 <- paste0("S", (nItems-miItems+1):nItems, "_g1")
    miVars_g2 <- paste0("S", (nItems-miItems+1):nItems, "_g2")
    
    pairs <- rbind(expand.grid(v1 = miVars_g1, v2 = miVars_g1, stringsAsFactors = FALSE),
                  expand.grid(v1 = miVars_g2, v2 = miVars_g2, stringsAsFactors = FALSE))
    pairs <- pairs[pairs$v1 != pairs$v2, ]  # drop self-pairs
    
    # pairs <- rbind(t(combn(miVars_g1, 2)), t(combn(miVars_g2, 2))) 
    for (iRow in seq_len(dim(pairs)[1])) {
      phi[pairs[iRow,1], pairs[iRow,2]] <- covNuanceFactor
    }
  }
  
  phi
}



# Creates a diagonal residual matrix of size nItems*nGroups + nOutcome
#   ... residuals of the outcome are zero as we use all latent formulation here 
#       irreducible error of the outcome is already modeled elsewhere
defineTheta <- function(varE, nItems, nGroups, nOutcome, nFactors = 1) {
  ## input:
  # varE      - [scalar] residual variance of the manifest variables
  # nItems    - [scalar] number of items per group
  # nGroups   - [scalar] number of groups
  # nOutcome  - [scalar] number of outcome variables        
  ## output:
  # theta     - [matrix] residual variance matrix of the items with items in rows and columns [items x items]
  ###
  if (nFactors != 1) {
    stop("Kim, you still need to implement this!")
  }
  return(diag(c(rep(varE, nItems * nGroups), rep(0, nOutcome * nGroups))))
  
}

defineBeta <- function(nLatentSEM, betaT, betaS, miItems, nItems, varNames, nFactors = 1) {
  ###
  # define matrix of path coefficients
  ## input:
  # nLatentSEM  - [scalar] number of factors (SEM view) to generate empty path coefficient matrix
  # betaT       - [scalar] coefficient of the true score composite on the outcome (in this case scalar)
  # betaS       - [scalar] coefficient of the specific item factors on the outcome 
  # miItems     - [scalar] how many items with specific item effect
  # nItems      - [scalar] number of items, i.e., specific item factors
  # varNames    - [vector] vector of factor names to set values in matrix (readable and easy)
  ## output:
  # betaMat     - [matrix] matrix of path coefficients with factors in rows and columns
  ###
  
  # # test it
  # condVec <- which(params$dgp$grid$P == 4 &
  #                    params$dgp$grid$cAlpha == 0.6 &
  #                    params$dgp$grid$varT == 0.5 &
  #                    params$dgp$grid$miItems == 0 &
  #                    params$dgp$grid$MIvio == "specific")
  # iCond <- condVec[1]
  # nLatentSEM <- (params$dgp$grid$P[iCond] + params$dgp$grid$Q[iCond] + 
  #                  params$dgp$nOutcome)*params$dgp$nGroups
  # betaT <- params$dgp$betaT
  # betaS <- params$dgp$grid$betaS[iCond]
  # nItems <- params$dgp$grid$P[iCond] 
  # miItems <- params$dgp$grid$miItems[iCond]
  # varNames <- colnames(lambdaSEM)
  
  # maintenance
  if (nFactors != 1) {
    stop("Kim, you still need to implement this!")
  }
  
  if (params$dgp$nGroups != 2) {
    stop("This function is only tested for two groups!")
  }
  
  # beta matrix without any effect (empty matrix)
  betaMat <- matrix(rep(0, nLatentSEM + nLatentSEM*(nLatentSEM-1)), ncol = nLatentSEM, byrow = T)
  # povide row and colnames to select cells of the matrix
  colnames(betaMat) <- rownames(betaMat) <- varNames 
  
  # outcome variables are sorted in the following order: Y1_g1 Y2_g1 Y2_g2 Y2_g2
  # either only the true factors load on the outcome or for MIvio == "specific", 
  #   miItem other specific item variances load on the outcome, too
  outVars <- varNames[grep("^Y", varNames)]
  for (iOut in seq_along(outVars)) {
    if (grepl("g1", outVars[iOut])) { # this is group 1; rowvector of loadings
      betaMat[outVars[iOut],] <- c(betaT, rep(c(0, betaS), c(nItems-miItems, miItems)), 
                                   rep(0, nItems + nFactors), 
                                   rep(0, params$dgp$nOutcome*params$dgp$nGroups))
    } else { # this is group 2
      # T & S 
      betaMat[outVars[iOut],] <-c(rep(0, nItems + nFactors), 
                                  betaT, rep(c(0, betaS), c(nItems-miItems, miItems)), 
                                  rep(0, params$dgp$nOutcome*params$dgp$nGroups))
    }
  }
  
  betaMat
}

getSigma <- function(lambda, beta, phi, theta) {
  ###
  # calculate model implied variance-covariance matrix
  ## input:
  # lambda    - [matrix] factor loading matrix with items in rows and factors in columns [items x factor]
  # beta      - [matrix] path coefficients with factors in rows and columns [factor x factor]
  # phi       - [matrix] variance-covariance matrix with factors in rows and columns [factor x factor]
  # theta     - [matrix] residual variance matrix (often diagonal matrix) with items in rows and columns [items x items]
  ## output:
  # sigma     - [matrix] variance-covariance matrix of the items with items in rows and columns [items x items]
  ###
  
  # check dimensionality of matrices
  nF <- ncol(lambda)
  if (!all(c(dim(beta), dim(phi)) == nF) && all(dim(theta) == nrow(lambda))) {
    stop("dimensionality of matrices does not fit. Check your input matrices!")
  } else if (!all(sapply(list(lambda, beta, phi, theta), is.numeric))) {
    stop("non numeric input in matrices. Check your input matrices!")
  }
  
  sigma <- lambda %*% solve(diag(nF) - beta) %*% phi %*%
    t(solve(diag(nF) - beta)) %*% t(lambda) + theta
  
  colnames(sigma) <- rownames(sigma) <- rownames(lambda)
  
  sigma
}

getMeanYg2 <- function(targetRho, varGroup, varWithinY) {
  ###
  # calculate mean for the second group of the outcome (assuming mean of group 1 is zero)
  #   -> mean differences will generate correlation between outcome and MI violation 
  #       inducing group
  ## input:
  # rho         - [scalar] target correlation (by simulating different means for subgroups)
  # varGroup    - [scalar] variance of the grouping variable (fixed at 0.25 due to balanced 
  #                        design and dichotomous grouping variable)
  # varWithinY  - [scalar] within group variance of the outcome 
  #                        var(y) = var(within group) + var(between group); 
  ## output:
  #             - [scalar] mean of the second group in the outcome 
  ###  
  if (varGroup != 0.25) {
    warning("Is your grouping variable binary and balanced? If not: check the formula!")
  }
  
  # warning("Variance within Y is only an estimate! Can I derive this values and precalculate the mean?")
  
  sqrt((targetRho**2 * varGroup * varWithinY)/(0.25**2 - targetRho**2 * varGroup * 0.25))
}

# # test it
# getMeanYg2(0, params$dgp$varGroup, varWithinY)
# getMeanYg2(0.2, params$dgp$varGroup, varWithinY)
# getMeanYg2(0.5, params$dgp$varGroup, varWithinY)

# this function needs to enable ...
#     ... mean differences for scalar MI violations and 
#     ... mean differences between groups in the outcome to simulate their correlation
# instead of switching between different kinds of mean vector structures, just
#   supply identical values if scalar MI is not violated or correlation between y and group is 0
defineMeanVec <- function(mItems, mOut, nItems, nGroups, nOutcome, 
                          mXg2 = mItems, mYg2 = mOut, miItems = 0, nFactors = 1) {
  # maintenance
  if (nFactors != 1) {
    stop("Kim, you still need to implement this!")
  }
  
  # user warning for probably unintential input
  if (mXg2 != miItems & miItems == 0) {
    warning("There are no items violating MI assumptions. Is this what you want?!")
  }
  
  M <- c(rep(mItems, nItems), # manifest variables group 1
         rep(c(mItems, mXg2), c(nItems-miItems, miItems)), # manifest variables group 2
         rep(mOut, nOutcome), rep(mYg2, nOutcome)) # outcome variable means
  
  return(M)
}






