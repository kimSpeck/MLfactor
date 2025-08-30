# sample data for different levels of measurement invariance violations
# ... specific item effects, (configural), metric, scalar, (residual)

# function gets full environment of the 01_simulateData script; thus no arguments
#   besides different process arguments
sampleData <- function() {
  
  # generate folder structure for this condition
  conditionFolder <- paste0("/simData_MIvio", MIvio, "_miItems", miItems,
                            "_P", P, "_Q", Q, "_cAlpha", cAlpha, "_varT", varT, 
                            "_corGO", corGO, "_R2", R2)
  createFolder(paste0(dataFolder, conditionFolder))
  
  # sample data in parallel: 
  #   generate samples in parallel as samples are drawn as random & independent
  data <- parLapply(cl, seq_len(params$sim$nSamples), function(iSample) {
    
    if (iSample > params$sim$nTrain) {
      # test sample with fixed sample size across all simulated conditions (independent of N)
      N <- params$sim$testN # N = 1000 for holdout test data 
    } else {
      N <- max(params$dgp$trainN)
    }
    
    # simulate data + save parameter checks
    checkFlag <- T 
    if (checkFlag) {N <- 20000}
    
    # simulate date from multivariate normal distribution
    # if we have this block matrix structure we only need N/nGroups observations
    #   as variables for each group are simulated separately and otherwise we would
    #   end up with N * nGroups observations
    
    dataWide <- mvtnorm::rmvnorm(N/params$dgp$nGroups, mean = M, sigma = sigmaSEM)
  
    # name columns of the data
    colnames(dataWide) <- c(paste0(rep(c(sapply(seq(P), 
                                                function(iP) {paste0("X", iP)})), times = params$dgp$nGroups), 
                                   "_g", rep(seq_len(params$dgp$nGroups), each = P)), 
                            "Y_g1", "Y_g2")
    
    # wide to long
    data <- data.frame(dataWide)
    data <- tidyr::pivot_longer(data, cols = !c(Y_g1, Y_g2), names_to = "item", values_to = "value")
    data <- tidyr::separate(data, col = item, into = c("item", "group"), sep = "_")
    data <- tidyr::pivot_wider(data, names_from = item, values_from = value)
    data$Y <- ifelse(data$group == "g1", data$Y_g1, data$Y_g2)
    data$group <- ifelse(data$group == "g1", 0, 1)
    
    # calculate score and average depending on P 
    itemNames <- sapply(seq_len(P), function(x) paste0("X", x))
    data$score <- rowSums(data[,itemNames])
    data$avg <- data$score/length(itemNames)
    
    
    # what to save form data simulation?
    # ... what do I need to fit the models on the data?
    # ... what do I need to compute outcome measures?
    # save ...
    #     ... data 
    #     ... 
    #     ... any checks?
    if (!checkFlag) {
      dataList <- list()   
    }
    
    
    # here: maybe save checkList parameters differently?
    if (checkFlag) {
      # check variances of the manifest variables
      var_hat <- apply(data, 2, var)
      
      # check correlations between group and outcome
      corGO_hat <- cor(data["Y"], data["group"])
      
      # check SEM data
      # function that creates lavaan syntax automatically with respective 
      #       number of items and variance of the true score
      
      m0 <- genModel(itemNames = itemNames, varT = varT)  
      
      # wenn ich die Varianz im Faktor auf die simulierte Varianz festsetze und die Markiervariable befreie (NA*X1), dann ...
      # ... entsprechen die geschätzten Ladungen im SEM den simulierten Faktorladungen
      # ... sind die Residualvarianzen der Items genau varS + varE (also orange + rot)
      # ... sind die regressionskoeffizienten der Items λ^2 ∗var(η_T)
      # ... die Residualvarianz im Outcome ist auch 1 (wie reinsimuliert)
      outSEM <- tryCatch({
        # with identical factor models, modelling both groups at once works as well
        # lavaan::lavOptions("optim.force.converged")
        fit_singleGroup <- lavaan::sem(m0, data = data, optim.force.converged = TRUE)
        # summary(fit_singleGroup, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
        
        # vector of estimated factor loadings (= flT)
        lambda_hat <- fit_singleGroup@Model@GLIST$lambda[seq_len(length(itemNames)),1]
        
        # residual variances of the items (= varS + varE)
        theta_hat <- diag(fit_singleGroup@Model@GLIST$theta)[seq_len(length(itemNames))]
        
        # residual variance of the outcome (= sigmaEy)
        sigmaEy_hat <- fit_singleGroup@Model@GLIST$psi[Q+1, Q+1]
        
        # R² factor -- items (= λ^2 ∗var(η_T))
        R2_itemsT <- lavInspect(fit_singleGroup, what = "rsquare")[itemNames]
        
        # # structure model
        # regression coefficient: factor -- outcome (= betaT)
        betaT_hat <- fit_singleGroup@Model@GLIST$beta[Q+1, Q]
        
        # get R² from lavaan to check betaT manipulation (= R2)
        R2_hat <- lavInspect(fit_singleGroup, what = "rsquare")[["fY"]]
        
        # check Cronbach's Alpha
        cAlpha_semTools <- semTools::reliability(fit_singleGroup)["alpha",]
        
        list(lambda_hat = lambda_hat, theta_hat = theta_hat, sigmaEy_hat = sigmaEy_hat, 
             R2_itemsT = R2_itemsT, betaT_hat = betaT_hat, R2_hat = R2_hat, 
             cAlpha_semTools = cAlpha_semTools)
      }, error = function(e) {
        list(lambda_hat = rep(NA, P), theta_hat = rep(NA, P), sigmaEy_hat = NA, 
             R2_itemsT = rep(NA, P), betaT_hat = NA, R2_hat = NA, 
             cAlpha_semTools = NA)
      })
      
      
      cAlpha_myFUN <- getCronbachsAlpha(data[,itemNames])
      
      
      outSEMgroups <- tryCatch({
        # for the multiple group analysis, the data looks fine in both groups (see first test above)
        fit_multipleGroups <- lavaan::sem(m0, data = data[,-c(1, 2)], group = "group")
        # fit_multipleGroups@Model@GLIST
        # summary(fit_multipleGroups, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
        
        # regression coefficient: factor -- outcome (= betaT)
        betaIdx <- which(names(fit_multipleGroups@Model@GLIST) %in% "beta")
        betaT_hat_g1 <- fit_multipleGroups@Model@GLIST[[betaIdx[1]]][Q+1, Q]
        betaT_hat_g2 <- fit_multipleGroups@Model@GLIST[[betaIdx[2]]][Q+1, Q]
        
        # get R² from lavaan to check betaT manipulation (= R2)
        R2_hat_g1 <- lavInspect(fit_multipleGroups, what = "rsquare")$"0"[["fY"]]
        R2_hat_g2 <- lavInspect(fit_multipleGroups, what = "rsquare")$"1"[["fY"]]
        
        c(betaT_hat_g1 = betaT_hat_g1, betaT_hat_g2 = betaT_hat_g2, 
          R2_hat_g1 = R2_hat_g1, R2_hat_g2 = R2_hat_g2)
      }, error = function(e) {
        c(betaT_hat_g1 = NA, betaT_hat_g2 = NA, R2_hat_g1 = NA, R2_hat_g2 = NA)
      })
      
      
      # save data and checks {variances, CronbachsAlpha}
      lm_modelSyntax <- paste0("Y ~ (1 +", paste(itemNames, sep = "", collapse = "+"), ")") # all items
      
      fitItems <- lm(lm_modelSyntax, data = data)
      R2_items <- summary(fitItems)$r.squared
      
      fitScore <- lm(Y ~ score, data = data)
      R2_score <- summary(fitScore)$r.squared
      
      fitAvg <- lm(Y ~ avg, data = data)
      R2_avg <- summary(fitAvg)$r.squared
      
      checkList <- list(condition = c(P = P, Q = Q, cAlpha = cAlpha, varT = varT,
                                      miItems = miItems, corGO = corGO, R2 = R2, 
                                      MIvio = MIvio),
                        # vectors
                        var_hat = var_hat,
                        lambda_hat = outSEM[["lambda_hat"]],
                        theta_hat = outSEM[["theta_hat"]], 
                        R2_itemsT = outSEM[["R2_itemsT"]],
                        # scalars
                        # correlation between outcome and group (corGO)
                        trueVShat = c(corGO_sim = corGO, corGO_hat = corGO_hat,
                        # single group SEM
                        # residual variance of the outcome (sigmaEy)
                        sigmaEy = params$dgp$sigmaEy, sigmaEy_hat = outSEM[["sigmaEy_hat"]],
                        # regression coefficient factor - outcome (betaT)
                        betaT = betaT, betaT_hat = outSEM[["betaT_hat"]], 
                        # R2 factor - outcome (R2)
                        R2_sim = R2, R2_hat = outSEM[["R2_hat"]], 
                        # Cronbach's Alpha (cAlpha) 
                        cAlpha_sim = cAlpha, 
                        cAlpha_myFUN = cAlpha_myFUN, cAlpha_semTools = outSEM[["cAlpha_semTools"]],
                        # multiple group SEM estimates
                        # regression coefficients factor -- outcome
                        betaT_hat_g1 = outSEMgroups["betaT_hat_g1"], betaT_hat_g2 = outSEMgroups["betaT_hat_g2"], 
                        # rsquared factor -- outcome
                        R2_hat_g1 = outSEMgroups["R2_hat_g1"], R2_hat_g2 = outSEMgroups["R2_hat_g2"]),
                        R2_fit = c(R2_items = R2_items, R2_score = R2_score, R2_avg = R2_avg))
    }
    
    # save files
    if (!checkFlag){
      sampleNames <- c(seq_len(params$sim$nTrain), paste0("test", seq_len(params$sim$nTest)))
      fileName <- paste0("sample_", sampleNames[iSample], ".rda")
      save(dataList, file = paste0(dataFolder, conditionFolder, "/", fileName))
    } else {
      checkList
    }
    
  })
  
  # only save check measures (this does not save the generated data)
  if (checkFlag) {
    fileName <- paste0("checkData_MIvio", MIvio, "_miItems", miItems, "_P", P, 
                       "_Q", Q, "_cAlpha", cAlpha, "_varT", varT, "_corGO", corGO, 
                       "_R2", R2, ".rda")
    save(data, file = paste0(dataFolder, conditionFolder, "/", fileName))
  }
  
}



