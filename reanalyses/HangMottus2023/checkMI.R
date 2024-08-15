# check MI violations
# scalar measurement invariance might be able to provoke superiority of item based 
#   prediction over factor based prediction
# it is widely known that personality inventories often violate scalar measurement
#   invariance assumptions (Dong & Dumas, 2020)
# check MI violations in the personality data used for prediction of general and 
#   specific psychopathology factors
# 
# is superiority of item based prediction due to measurement invariance violations?
#   -> cannot unambigiously answer this question but we can test if measurement invariance
#       violations are present and therefore can provide an alterantive explanation of the results
library(lavaan)
library(tidyverse)

################################################################################
##### data preprocessing #####
################################################################################
# data preprocessing as in original code (identical subset of cases and missing values)
# load data
adata_full = read.csv("Dataset/AAAData (1).csv")

#clean data----
## drop data with more than 40 missing personality item responses
adata_full = adata_full %>% filter( (adata_full %>% select(210:457) %>% is.na %>% rowSums) <= 40 )

## quality check 
adata_full = adata_full %>% select(-contains("SS_"))

#substance abuse scores ASSIST questionnaire (never used is coded as 0)
adata_full$tobacco <- select(adata_full, c("assist_q2_tobacco", "assist_q3_tobacco", "assist_q4_tobacco", "assist_q5_tobacco", "assist_q6_tobacco", "assist_q7_tobacco")) %>% rowMeans
adata_full$alcohol <- select(adata_full, c("assist_q2_alcohol", "assist_q3_alcohol", "assist_q4_alcohol", "assist_q5_alcohol", "assist_q6_alcohol", "assist_q7_alcohol")) %>% rowMeans
adata_full$cannabis <- select(adata_full, c("assist_q2_cannabis", "assist_q3_cannabis", "assist_q4_cannabis", "assist_q5_cannabis", "assist_q6_cannabis", "assist_q7_cannabis")) %>% rowMeans

adata_full[,451:453] = adata_full %>% select(451:453) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), 0)))
colnames(adata_full)

## restrict age to 18 plus 
adata = adata_full %>% filter(agegroup == 2 & age > 17) %>% select(-contains("SS_"))

## median replace
adata[,210:449] = adata %>% select(210:449) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))

## reverse personality items, calculate alpha
p = select(adata, 210:449) 

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

# alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
#   map(~psych::alpha(p %>% select(grep(.x, names(p)))))
# 
# map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

##### check factor structure (measurement model) #####
## items, domains and facets
items = p
facet.names <- paste(rep(c("N","E","O","A","C"), each = 6), rep(1:6, times = 5), sep= "")
domain.names <- c("N","E","O","A","C")
facets <- as.data.frame(matrix(ncol = 30, nrow = nrow(p))) %>% `names<-`(facet.names)
domains <- as.data.frame(matrix(ncol = 5, nrow = nrow(p))) %>% `names<-`(domain.names)
for(i in facet.names) 
  facets[,i] = p %>% select(contains(i)) %>% rowMeans(na.rm = TRUE)
for (i in domain.names)
  domains[,i] = p %>% select(contains(i)) %>% rowMeans(na.rm = TRUE)

##### find possible group variables for MI testing #####
# those are variables that are used as an outcome
# therefore do not search for measurement invariance due to these variables 
data <- adata[,c(27:46, 451:453)]
colnames(data)

# age as a grouping variable 
# only two age groups; maybe use different age group according to MI literature on
#   MI violations due to age as a covariate
unique(adata_full$agegroup)
unique(adata_full$age)

# basically two levels
unique(adata_full$sex)

# length of education in years? 
unique(adata_full$edu_length) 

unique(adata_full$job) # don't know what numbers mean

p = select(adata, c("age", "edu_length", "sex", "job", 210:449))
table(p$sex) # one group is more than twice the size of the other group
table(p$edu_length) # maybe split in the middle and simply generate 2 artificial edu_length groups
table(p$job)# weird
hist(p$age, breaks = seq(15, 800, 1), xlim = range(15, 60))
################################################################################
# fit general factor model
################################################################################
##### first order factor model #####
# (check if factor loadings of items are high enough)
# build second order factors on top
facet.names
colnames(items)

genModel <- function(factorNames, items) {
  varsPerFac <- sapply(factorNames, function(iFacet) {
    tmp = colnames(p %>% select(contains(iFacet)))
    paste0(sapply(tmp, function(iVar) {
      paste0(iVar, " +")}, simplify = "array"), collapse = " ")
  })
  cfaModel <- c()
  for (iFactor in factorNames) {
    (cfaModel <- paste0(cfaModel, "\n", 
                        iFactor, " =~ ", 
                        substr(varsPerFac[[iFactor]], 1, nchar(varsPerFac[[iFactor]])-2)))
  }
  return(cfaModel)
}
facetModel <- genModel(facet.names, p)

# this is too huge for my RAM if standard errors are estimated in the model, too
# fitFacets <- lavaan::cfa(facetModel, data = items, meanstructure = T)
covItems <- cov(items)
sampleSize <- dim(items)[1]
fitFacets <- lavaan::cfa(facetModel, # estimation based on cov matrix
                         sample.cov = covItems, 
                         sample.nobs = sampleSize,
                         se = "none")
# warning: not positive definite due to highly correlated facet factors since 
#   higher order factors are missing so far

# next: implement higher order factors
summary(fitFacets, fit.measures = T, stand = T)

# fitFacets@Model@GLIST
# diag(fitFacets@Model@GLIST$theta)

##### second order factor model #####
factorNames <- list(first = facet.names, 
                    second = domain.names)

genModel <- function(factorNames, items) {
  varsPerFac <- sapply(factorNames[['first']], function(iFacet) {
    tmp = colnames(items %>% select(contains(iFacet)))
    paste0(sapply(tmp, function(iVar) {
      paste0(iVar, " +")}, simplify = "array"), collapse = " ")
  })
  cfaModel <- c()
  for (iFactor in factorNames[['first']]) {
    (cfaModel <- paste0(cfaModel, "\n", 
                        iFactor, " =~ ", 
                        substr(varsPerFac[[iFactor]], 1, nchar(varsPerFac[[iFactor]])-2)))
  }
  facPerDom <- sapply(factorNames[['second']], function(iDomain) {
    tmp = stringr::str_subset(factorNames[['first']], iDomain)
    paste0(sapply(tmp, function(iFac) {
      paste0(iFac, " +")}, simplify = "array"), collapse = " ")
  })
  for (iDomain in factorNames[['second']]) {
    (cfaModel <- paste0(cfaModel, "\n", 
                        iDomain, " =~ ", 
                        substr(facPerDom[[iDomain]], 1, nchar(facPerDom[[iDomain]])-2)))
  }
  return(cfaModel)
}
facetModel <- genModel(factorNames, p)

# fit CFA with second order factors based on covariance matrix 
fitDomains <- lavaan::cfa(facetModel, # estimation based on cov matrix
                         sample.cov = covItems, 
                         sample.nobs = sampleSize,
                         se = "none")
summary(fitDomains, fit.measures = T, stand = T)

# Comparative Fit Index (CFI)                    0.601
# Tucker-Lewis Index (TLI)                       0.597
# RMSEA                                          0.042
# SRMR                                           0.101

# fit again based on actual data
fitDomainsData <- lavaan::cfa(facetModel, # estimation based on cov matrix
                          data = items,
                          meanstructure = T,
                          se = "none")
summary(fitDomainsData, fit.measures = T, stand = T)

################################################################################
# MI for sex
################################################################################
## reverse personality items, calculate alpha
# p = select(adata, c("age", "sex", "edu_length", 210:449)) 
p = select(adata, c("sex", 210:449)) 

# for every MI check model the following warning is displayed
# Warnmeldung:
#   lavaan->lav_object_post_check():  
#   some estimated lv variances are negative 

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

# alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
#   map(~psych::alpha(p %>% select(grep(.x, names(p)))))
# 
# map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

##### check factor structure (measurement model) 
## items, domains and facets
items = p

items <- items[which(!is.na(items$sex)), ]

# # in every model for separate groups:
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# configural invariance
fit.config <- lavaan::cfa(facetModel, data = items, meanstructure = T, se = "none",
                    group = "sex")
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative
summary(fit.config, fit.measures = T, stand = T)
# fitMeasures(fit1, se = "none") # this kills RAMS

# metric (weak) invariance
fit.metric <- cfa(facetModel, data = items, meanstructure = T, se = "none",
            group = "sex",
            group.equal = "loadings")
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# scalar (strong) invariance
fit.scalar <- cfa(facetModel, data = items, meanstructure = T, se = "none",
            group = "sex",
            group.equal = c("intercepts", "loadings"))
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# residual invariance
fit.residual <- cfa(facetModel, data = items, meanstructure = T, se = "none",
            group = "sex",
            group.equal = c("residuals", "intercepts", "loadings"))
# Warnmeldung:
#   lavaan->lav_object_post_check():  
#   some estimated lv variances are negative 

# model comparison tests
lavTestLRT(fitDomainsData, fit.config, fit.metric, fit.scalar, fit.residual)
lavTestLRT(fit.config, fit.metric, fit.scalar, fit.residual)

# where do negative latent variable variances come from?
# check missing values or impossible values in the data
any(is.na(items))
any(items > 4)
any(items < 0)
checkData <- psych::describe(items)

##### fit models for single higher order factors 
# generate cfa models for every domain
# factor model contains one domain factor which is loading on all facettes of the domain
#   and all items of the respective domain are on the lowest level
highOrderFactors <- 5
nFacPerDom <- 6
factorNames <- vector(mode = "list", length = highOrderFactors)
names(factorNames) <- domain.names

for (i in seq(factorNames)) {
  factorNames[[i]] <- vector(mode = "list", length = nFacPerDom)
  facetteNames <- stringr::str_subset(facet.names, names(factorNames)[i]) 
  names(factorNames[[i]]) <- facetteNames
  for (iFacet in facetteNames){
    tmp = colnames(items %>% select(contains(iFacet)))
    tmpFacette <- paste0(sapply(tmp, function(iVar) {
      paste0(iVar, " +")},  simplify = "array"), collapse = " ")
    factorNames[[i]][[iFacet]] <- substr(tmpFacette, 1, nchar(tmpFacette)-2)
  }
}


factorModelSyntax <- vector(mode = "list", length = highOrderFactors)
for (iDomain in seq(factorNames)) {
  oneDomModel <- c()
  for (iFacette in seq_len(nFacPerDom)) {
    oneDomModel <- paste0(oneDomModel, "\n", 
                          names(factorNames[[iDomain]])[iFacette], " =~ ", 
                          factorNames[[iDomain]][[iFacette]])
  }
  allFacettes <- paste0(names(factorNames[[iDomain]]), sep = " + ", collapse = " ")
  oneDomModel <- paste0(oneDomModel, "\n", 
                        names(factorNames)[iDomain], " =~ ", 
                        substr(allFacettes, 1, nchar(allFacettes)-3))
  factorModelSyntax[[iDomain]] <- oneDomModel
}
names(factorModelSyntax) <- domain.names

# test measurement invariance violations for each of these models
results <- vector(mode = "list", length(factorModelSyntax))
for (iModel in seq(factorModelSyntax)){
  # generate subset of items with only respective variables in the data
  tmpDom <- names(factorModelSyntax)[iModel]
  itemsSub <- items %>% select("sex", contains(tmpDom))
  
  # fit cfa 
  tmpModel <- factorModelSyntax[[iModel]]
  
  # without groups
  fitDomainsData <- lavaan::cfa(tmpModel, data = itemsSub,
                                meanstructure = T,
                                se = "none")
  
  # configural invariance
  fit.config <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                            group = "sex")
  # metric (weak) invariance
  fit.metric <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                    group = "sex",
                    group.equal = "loadings")
  
  # scalar (strong) invariance
  fit.scalar <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                    group = "sex",
                    group.equal = c("intercepts", "loadings"))
  
  # residual invariance
  fit.residual <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                      group = "sex",
                      group.equal = c("residuals", "intercepts", "loadings"))
  
  # model comparison tests
  results[[iModel]] <- lavTestLRT(fitDomainsData, fit.config, fit.metric, fit.scalar, fit.residual)
}
# every MI test becomes significant for every domain in the data

################################################################################
# MI for education length
################################################################################
p = select(adata, c("edu_length", 210:449)) 

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

# alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
#   map(~psych::alpha(p %>% select(grep(.x, names(p)))))
# 
# map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

##### check factor structure (measurement model) 
## items, domains and facets
items = p
items <- items[which(!is.na(items$edu_length)), ]

table(items$edu_length)
items$eduGroup <- ifelse(items$edu_length <= 5, 1, 2)
items$edu_length <- NULL
table(items$eduGroup)

# # in every model for separate groups:
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# configural invariance
fit.config <- lavaan::cfa(facetModel, data = items, meanstructure = T, se = "none",
                          group = "eduGroup")
# fitMeasures(fit1, se = "none") # this kills RAMS

# metric (weak) invariance
fit.metric <- cfa(facetModel, data = items, meanstructure = T, se = "none",
                  group = "eduGroup",
                  group.equal = "loadings")
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# scalar (strong) invariance
fit.scalar <- cfa(facetModel, data = items, meanstructure = T, se = "none",
                  group = "eduGroup",
                  group.equal = c("intercepts", "loadings"))
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# residual invariance
fit.residual <- cfa(facetModel, data = items, meanstructure = T, se = "none",
                    group = "eduGroup",
                    group.equal = c("residuals", "intercepts", "loadings"))
# Warnmeldung:
#   In lav_object_post_check(object) :
#   lavaan WARNING: some estimated lv variances are negative

# model comparison tests
lavTestLRT(fitDomainsData, fit.config, fit.metric, fit.scalar, fit.residual)

##### fit models for single higher order factors 
# generate cfa models for every domain
# factor model contains one domain factor which is loading on all facettes of the domain
#   and all items of the respective domain are on the lowest level
highOrderFactors <- 5
nFacPerDom <- 6
factorNames <- vector(mode = "list", length = highOrderFactors)
names(factorNames) <- domain.names

for (i in seq(factorNames)) {
  factorNames[[i]] <- vector(mode = "list", length = nFacPerDom)
  facetteNames <- stringr::str_subset(facet.names, names(factorNames)[i]) 
  names(factorNames[[i]]) <- facetteNames
  for (iFacet in facetteNames){
    tmp = colnames(items %>% select(contains(iFacet)))
    tmpFacette <- paste0(sapply(tmp, function(iVar) {
      paste0(iVar, " +")},  simplify = "array"), collapse = " ")
    factorNames[[i]][[iFacet]] <- substr(tmpFacette, 1, nchar(tmpFacette)-2)
  }
}


factorModelSyntax <- vector(mode = "list", length = highOrderFactors)
for (iDomain in seq(factorNames)) {
  oneDomModel <- c()
  for (iFacette in seq_len(nFacPerDom)) {
    oneDomModel <- paste0(oneDomModel, "\n", 
                          names(factorNames[[iDomain]])[iFacette], " =~ ", 
                          factorNames[[iDomain]][[iFacette]])
  }
  allFacettes <- paste0(names(factorNames[[iDomain]]), sep = " + ", collapse = " ")
  oneDomModel <- paste0(oneDomModel, "\n", 
                        names(factorNames)[iDomain], " =~ ", 
                        substr(allFacettes, 1, nchar(allFacettes)-3))
  factorModelSyntax[[iDomain]] <- oneDomModel
}
names(factorModelSyntax) <- domain.names

# test measurement invariance violations for each of these models
resultsEdu <- vector(mode = "list", length(factorModelSyntax))
for (iModel in seq(factorModelSyntax)){
  # generate subset of items with only respective variables in the data
  tmpDom <- names(factorModelSyntax)[iModel]
  itemsSub <- items %>% select("eduGroup", contains(tmpDom))
  
  # fit cfa 
  tmpModel <- factorModelSyntax[[iModel]]
  
  # without groups
  fitDomainsData <- lavaan::cfa(tmpModel, data = itemsSub,
                                meanstructure = T,
                                se = "none")
  
  # configural invariance
  fit.config <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                            group = "eduGroup")
  # metric (weak) invariance
  fit.metric <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                            group = "eduGroup",
                            group.equal = "loadings")
  
  # scalar (strong) invariance
  fit.scalar <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                            group = "eduGroup",
                            group.equal = c("intercepts", "loadings"))
  
  # residual invariance
  fit.residual <- lavaan::cfa(tmpModel, data = itemsSub, meanstructure = T, se = "none",
                              group = "eduGroup",
                              group.equal = c("residuals", "intercepts", "loadings"))
  
  # model comparison tests
  resultsEdu[[iModel]] <- lavTestLRT(fitDomainsData, fit.config, fit.metric, fit.scalar, fit.residual)
}

# this is the only domain which shows at least configural MI
resultsEdu[[3]]

################################################################################
# replicate predictive model
################################################################################
# replicate predictive models
# elastic net regression models with 10-fold cross-validation via the glmnet package
#   -> what are alpha parameters?
# selecting the regularisation parameter (‘lambda min’) 
# 75–25% training–validation sample split to avoid overfitting 
# repeated the training-validation procedure for 100 random 75–25 sample splits 
#   and computed the mean and standard deviation for predictive accuracies across these splits.

# analyse the associations between personality traits and the three specific psychopathology factors.
# correlated predicted values of each psychopathology factor from the elastic net 
#     regressions (on each level of the personality hierarchy) with the observed values 
#     of all other psychopathology factors

##### prediction from domains #####

##### prediction from facets #####

##### prediction from items #####




################################################################################
# improve prediction by adding covariates 
################################################################################
