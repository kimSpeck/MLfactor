# 1. check analysis (more comments)
# 2. replicate analysis
# 3. compare results with different predictor sets (add age, add sex, add both)
#     -> compare performance of domains, facets and items against each other for 
#         different predictor sets

# nItems = 240
# nFacets = 30
# nDomains = 5

source("analysisTools.R")

#setup----
library(tidyverse)
library(lavaan)
library(psych)
library(careless)
library(EFAtools)
library(glmnet)
# 
library(bootnet)
library(qgraph)
#library(sjPlot)
library(matrixStats)

adata_full = read.csv("Dataset/AAAData (1).csv")

#clean data----
## drop data with more than 40 missing personality item responses
adata_full = adata_full %>% filter( (adata_full %>% select(210:457) %>% is.na %>% rowSums) <= 40 )

## quality check 
adata_full = adata_full %>% select(-contains("SS_"))

# substance abuse scores ASSIST questionnaire (never used is coded as 0)
# substance abuse scores build the third factor in the bifactor model later that 
#     will be used as a criterion/target predicted by personality items/features 
adata_full$tobacco <- select(adata_full, c("assist_q2_tobacco", "assist_q3_tobacco", "assist_q4_tobacco", "assist_q5_tobacco", "assist_q6_tobacco", "assist_q7_tobacco")) %>% rowMeans
adata_full$alcohol <- select(adata_full, c("assist_q2_alcohol", "assist_q3_alcohol", "assist_q4_alcohol", "assist_q5_alcohol", "assist_q6_alcohol", "assist_q7_alcohol")) %>% rowMeans
adata_full$cannabis <- select(adata_full, c("assist_q2_cannabis", "assist_q3_cannabis", "assist_q4_cannabis", "assist_q5_cannabis", "assist_q6_cannabis", "assist_q7_cannabis")) %>% rowMeans

# recode missing values as never used = 0
adata_full[,451:453] = adata_full %>% select(451:453) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), 0))) 

## restrict age to 18 plus 
adata = adata_full %>% filter(agegroup == 2 & age > 17) %>% select(-contains("SS_"))

## median replace (impute missing values with median values)! 
colnames(adata[,210:449]) # these are the items! imputation for missing values in items
sum(is.na(adata[,210:449])) # 406 missing values in the data

##### !!! this is leakage!!!
# due to more leakage on item level compared to facet or domain level, item based
#   prediction might show more overfitting of the data which might increase the test performance
adata[,210:449] = adata %>% select(210:449) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))

## reverse personality items, calculate alpha
p = select(adata, 210:449) 

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
  map(~psych::alpha(p %>% select(grep(.x, names(p)))))

map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

# calculate factes and domains as sumscores----
## items, domains and facets
items = p
facet.names <- paste(rep(c("N","E","O","A","C"), each = 6), rep(1:6, times = 5), sep= "")
domain.names <- c("N","E","O","A","C")
facets <- as.data.frame(matrix(ncol = 30, nrow = nrow(p))) %>% `names<-`(facet.names)
domains <- as.data.frame(matrix(ncol = 5, nrow = nrow(p))) %>% `names<-`(domain.names)

# unweighted sum scores to calculate facets and domains
varInfac <- list()
for(i in seq_along(facet.names)) {
  #varInfac[[i]] <- names(p %>% select(contains(facet.names[i])))
  facets[,facet.names[i]] = p %>% select(contains(facet.names[i])) %>% rowMeans(na.rm = TRUE)
  # print(names(p %>% select(contains(facet.names[i])))) # check
}

for (i in domain.names) {
  domains[,i] = p %>% select(contains(i)) %>% rowMeans(na.rm = TRUE)
  # print(names(p %>% select(contains(i)))) # check
}


#EFA----
# factors from EFA will become the criterion values to predict with items, facets and domains
#subset DSM data
data <- adata[,c(27:46, 451:453)]
colnames(data)

#test suitablity for factor analysis
KMO(data)
BARTLETT(data)

#select number of factors
PARALLEL(data, cor_method = "spearman") # PCA-determined eigenvalues:  3
VSS(data)
eigenvalues<-eigen(cor(data))$values
plot(eigenvalues, type = 'b', pch = 16,
     main = "Scree Plot", xlab="", ylab="Eigenvalues")
axis(1, at = 1:10, labels = 1:10)

#build bifactor model with 3 specific factors
efa_dsm <- psych::omega(data, fm="minres", nfactors = 3, poly = F)

#extract factors scores for further analysis
efa_scores <- as.data.frame(psych::factor.scores(data,efa_dsm$schmid$sl)$scores)

data$p_factor <- efa_scores$g # general factor from bifactor model
data$intern <- efa_scores$`F1*` 
data$psych <- efa_scores$`F2*` 
data$substance <- efa_scores$`F3*` # tobacco, alcohol, cannabis

#predict p_factor using personality----
sex = adata$sex
agegroup = adata$agegroup
age = adata$age
adata <- cbind(data, domains, facets, items, sex, agegroup, age)

table(adata$age %>% round(0))


Data_all <- cbind(data, domains, facets, items)
model <- list()
coef <- list()

# predict substance by domains with ENET
# alpha is not cross-validated
for (i in 1:100){
  model[[i]]<-cv.glmnet(as.matrix(Data_all[,c(28:32)]), Data_all$substance, alpha = .05)
  coef[[i]] <- round(coef(model[[i]], s = "lambda.min"),3)
} #repeat for all factors across all levels of the personality hierarchy

# investigate coefficients from ENET with domains as predictors
coef <- lapply(coef, as.matrix) 
coef <- do.call("cbind", coef)
estimate <- as.data.frame(rowMeans(coef)) # average coef values across 100 repetitions of the ENET
## results: 
#   more neuroticism, agreeableness, conscientiousness tend to reduce substance abuse
#   more openness and extraversion tend to increase substance abuse
#           (Intercept)        1.21854
#           N                 -0.15973
#           E                  0.22764
#           O                  0.12628
#           A                 -0.19547
#           C                 -0.50016
sd <- as.data.frame(rowSds(coef)) # standard error of coefficients
nuance <- cbind(estimate, sd)
nuance <- nuance %>% 
  arrange(desc(abs(rowMeans(estimate)))) %>% round(.,3) # averaged coefficients and SE in one table
domain_substance <- as.data.frame(cbind(nuance, CIlower = nuance$`rowMeans(coef)`-nuance$`rowSds(coef)`, Est = nuance$`rowMeans(coef)`, CIupper = nuance$`rowMeans(coef)` + nuance$`rowSds(coef)`))[,3:5]

#create and validate elastic net prediction models----
# train:test as 75:25 percent of the data
#   ... but missing values are already imputed for the full data set -> leakage
# 10-fold cross validation but ...
#   ... alpha is set to alpha = 0.05 (close to ridge)
#   ... no lambda vector explicitly given to the model
#   ... min instead of se as criterion 
foo1 = function(x1,y,x2=NULL,p=.75) {
  
  s = sample(nrow(x1), p*nrow(x1))
  if(is.null(x2)) x2 = x1
  
  x1 = sapply(x1, scale)
  x2 = sapply(x2, scale)
  
  tmp <- cv.glmnet(x1[s,], y[s], alpha = .05) %>% 
    predict(x2[-s,], s = "lambda.min") %>%  
    cor(y[-s])
  tmp**2
}

# for 4 different criterion values
# for only personality data, personality + sex, personality + age, personality + age + sex
# Run the models through domains, facets and items 
# Replicate 100 times

# ! sex has missing values
sum(is.na(sex)) # there is exactly one missing values
if (any(is.na(sex))) {# there is exactly one missing values)
  idxNA <- which(is.na(sex)) # remove missing value from domains, facets and items (and sex + age)
  domains <- domains[-idxNA,]
  facets <- facets[-idxNA,]
  items <- items[-idxNA,]
  adata <- adata[-idxNA,]
}
sex = ifelse(adata$sex == 1, 0, 1) # sex as dummy; standardising the dummy is fine
agegroup = adata$agegroup
age = adata$age

# criterion to predict with ENET
criterionVec <- list(adata$p_factor, adata$intern, adata$psych, adata$substance)
# which data matrix will we use for prediction
#   which columns will be removed from data in this specific iteration?
#   1. sex and age, 2. only sex, 3. only age, 4. no extra variables
dataVec <- list(NULL, "age", "sex", c("age", "sex")) # thus reduce data set further in every iteration

# # check if removal of variables works
# domains <- cbind(domains, sex, age)
# domains[,c(dataVec[[1]])] <- list(NULL)
# colnames(domains)
# domains[,c(dataVec[[4]])] <- list(NULL)
# colnames(domains)

colnames(domains)
colnames(facets)
colnames(items)

results <- lapply(seq_along(dataVec), function(iRM) {
  # add sex and age to the data
  domains <- cbind(domains, sex, age)
  facets <- cbind(facets, sex, age)
  items <- cbind(items, sex, age)
  
  # remove extra predictors for each iteration (depending on condition)
  domains[,c(dataVec[[iRM]])] <- list(NULL)
  facets[,c(dataVec[[iRM]])] <- list(NULL)
  items[,c(dataVec[[iRM]])] <- list(NULL)
  
  print(colnames(domains))
  
  TMP <- lapply(criterionVec, function(iCrit) {
    tmp <- list(domains, facets, items) %>% 
      map(~replicate(100, foo1(., iCrit, p=.75))) %>%
      bind_cols %>%
      describe %>% as.data.frame() %>% select(mean, sd) %>% t
    colnames(tmp) <- c("domains", "facets", "items")
    tmp
  })
  names(TMP) <- c("pFactor", "intern", "psych", "substance")
  
  # remove all extra variables before the next iteration for different data set starts
  domains[,c(dataVec[[4]])] <- list(NULL)
  facets[,c(dataVec[[4]])] <- list(NULL)
  items[,c(dataVec[[4]])] <- list(NULL)
  
  TMP
})

names(results) <- c("sexANDage", "onlySex", "onlyAge", "none")

# compare results with different predictor sets
pFactor <- rbindResults(results, "pFactor")
intern <- rbindResults(results, "intern")
psych <- rbindResults(results, "psych")
substance <- rbindResults(results, "substance")

predRes <- rbind(pFactor, intern, psych, substance)

idxCondVec <- c("sexANDage", "onlySex", "onlyAge", "none")

predRes <- rowNames2col(predRes)
predRes$rowNames <- stringr::str_replace(predRes$rowNames,  
                                         "[:digit:]{1,}$", "") # nachbessern
predRes$idxCond <- plyr::mapvalues(predRes$idxCond, 
                                   from = unique(predRes$idxCond), 
                                   to = idxCondVec)

# add original values
original <- cbind(domains = c(0.50, 0.03, 0.61, 0.03, 0.11, 0.04, 0.32, 0.03),
                  facets = c(0.51, 0.03, 0.67, 0.02, 0.22, 0.04, 0.41, 0.03), 
                  items = c(0.53, 0.02, 0.69, 0.02, 0.27, 0.04, 0.48, 0.03), 
                  idxCond = rep("original", 4*2),
                  measure = rep(c("pFactor", "intern", "psych", "substance"), each = 2),
                  rowNames = rep(c("mean", "sd"), times = length(criterionVec)))

predRes <- rbind(predRes, original)

str(predRes)

col2num <- c("domains", "facets", "items")
predRes[col2num] <- lapply(predRes[col2num], as.numeric)

col2fac <- c("idxCond", "measure", "rowNames")
predRes[col2fac] <- lapply(predRes[col2fac], factor)

# wide to long for domains, facets, items
predRes <- tidyr::pivot_longer(predRes, cols = c("domains", "facets", "items"),
                               names_to = "levels", values_to = "R2")

# long to wide for rowNames, i.e., mean & se
predRes <- tidyr::pivot_wider(predRes, names_from = rowNames, 
                              values_from = R2)

predRes$idxCond <- factor(predRes$idxCond, 
                          levels = c("original", "none", "onlyAge", "onlySex", "sexANDage"),
                          labels = c("original", "replicate", "+Age", "+Sex", "+Age+Sex"))

colValVec <- c("black", '#339900', '#0000FF', '#9900CC', '#FF0000')

# plot results
(pCompare <- ggplot(predRes, aes(x = levels, y = mean, 
                                 group = idxCond, colour = idxCond, linetype = idxCond)) +
    geom_point(size = 2) +
    geom_line() +
    scale_color_manual(values = colValVec) +
    ylab("average test R²") +
    geom_errorbar(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), width = 0.2) +
    facet_wrap(~ measure, scales = "free") + 
    theme(panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid', color = "grey"), 
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid', color = "grey"),
          panel.background = element_rect(color = "white", fill = "white"),
          axis.text.y = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15)))

save(predRes, file = "Dataset/replicateResultsData.rda")

plotFolder <- "plots"
if (!file.exists(plotFolder)){
  dir.create(plotFolder)
}

ggplot2::ggsave(filename = paste0(plotFolder, "/replicateResults.png"),
                plot = pCompare,
                width = 10.5,
                height = 9.41,
                units = "in")

# age zusätzlich zu domains oder facets aufzunehmen macht Prediction schlechter!
# mit {domains, facets, items} + age
# für p_factor:
#       mean 0.49760988 0.5035403 0.52357323
#       sd   0.03000685 0.0322614 0.03386979

# für intern:
#       mean 0.59878570 0.65475812 0.68837846
#       sd   0.03126554 0.03303537 0.02257346

# für psych:
#       mean 0.07215395 0.17665516 0.2516669
#       sd   0.04061207 0.04711811 0.0438922

# für substance:
#       mean 0.30232471 0.40385241 0.47000532
#       sd   0.04438687 0.03543191 0.02983622


# mit {domains, facets, items} + sex
# für p-factor
#       mean 0.49611923 0.50857456 0.52509687
#       sd   0.03297089 0.03002387 0.02992118

# für intern
#       mean 0.62885804 0.68392415 0.70217615
#       sd   0.02270342 0.02116176 0.02234138

# für psych
#       mean 0.09818933 0.22397509 0.27430425
#       sd   0.04561470 0.03946214 0.03830359

# für substance 
#       mean 0.35333274 0.43227869 0.47963990
#       sd   0.03652258 0.03307682 0.03261478

list(domains, facets, items) %>% 
  map(~replicate(100, foo1(., adata$p_factor, p=.75))) %>%
  bind_cols %>%
  describe %>% as.data.frame() %>% select(mean, sd) %>% t

list(domains, facets, items) %>% 
  map(~replicate(100, foo1(., adata$intern, p=.75))) %>%
  bind_cols %>%
  describe %>% as.data.frame() %>% select(mean, sd) %>% t

list(domains, facets, items) %>% 
  map(~replicate(100, foo1(., adata$psych, p=.75))) %>%
  bind_cols %>%
  describe %>% as.data.frame() %>% select(mean, sd) %>% t

list(domains, facets, items) %>% 
  map(~replicate(100, foo1(., adata$substance, p=.75))) %>%
  bind_cols %>%
  describe %>% as.data.frame() %>% select(mean, sd) %>% t

#test discriminant validity----

# function for cross prediction 
foo2 = function(x1,y,x2=NULL,p=.67, otheroutcomes) {
  val = !is.na(y)
  y = y[val]
  x1 = x1[val,]
  x2 = x2[val,]
  otheroutcomes = otheroutcomes[val,]
  s = sample(nrow(x1), p*nrow(x1))
  if(is.null(x2)) x2 = x1
  x1 = sapply(x1,scale)
  x2 = sapply(x2,scale)
  cv.glmnet(x1[s,], y[s], alpha = .05) %>%
    predict(x2[-s,], s = "lambda.min") %>%
    cor(otheroutcomes[-s,], use="pairwise") %>% as.numeric
}

outcomes1 <- adata[,24:27]
crossdom1 = crossfacet1 = crossitem1 = matrix(ncol=ncol(outcomes1),nrow=ncol(outcomes1))

for(i in seq_along(outcomes1)) {
  crossdom1[,i] =  list(domains) %>%
    map(~replicate(100, foo2(., outcomes1[,i], otheroutcomes = outcomes1[,]))) %>%
    map(~rowMeans(.x))  %>% unlist
  
  crossfacet1[,i] =  list(facets) %>%
    map(~replicate(100, foo2(., outcomes1[,i], otheroutcomes = outcomes1))) %>%
    map(~rowMeans(.x))  %>% unlist
  
  crossitem1[,i] =  list(items) %>%
    map(~replicate(100, foo2(., outcomes1[,i], otheroutcomes = outcomes1))) %>%
    map(~rowMeans(.x))  %>% unlist
}


lt = function(x) x[lower.tri(x)]

harmonize = function(X){
  m = rowMeans(cbind(lt(X), lt(t(X))))
  X[lower.tri(X)] = m
  X = t(X)
  X[lower.tri(X)] = m
  X
}

summary(lt(harmonize(cov2cor(crossdom1)))^2)
summary(lt(harmonize(cov2cor(crossfacet1)))^2)
summary(lt(harmonize(cov2cor(crossitem1)))^2)


#manhatten plots----
keys <- names(items)

keys_names <- names(items)
keys <-
  list(N1 = keys_names[grepl("^N1.*", keys_names)], N2 = keys_names[grepl("^N2.*", keys_names)], N3 = keys_names[grepl("^N3.*", keys_names)], N4 = keys_names[grepl("^N4.*", keys_names)], N5 = keys_names[grepl("^N5.*", keys_names)], N6 = keys_names[grepl("^N6.*", keys_names)], E1 = keys_names[grepl("^E1.*", keys_names)], E2 = keys_names[grepl("^E2.*", keys_names)], E3 = keys_names[grepl("^E3.*", keys_names)], E4 = keys_names[grepl("^E4.*", keys_names)], E5 = keys_names[grepl("^E5.*", keys_names)], E6 = keys_names[grepl("^E6.*", keys_names)], O1 = keys_names[grepl("^O1.*", keys_names)], O2 = keys_names[grepl("^O2.*", keys_names)], O3 = keys_names[grepl("^O3.*", keys_names)], O4 = keys_names[grepl("^O4.*", keys_names)], O5 = keys_names[grepl("^O5.*", keys_names)], O6 = keys_names[grepl("^O6.*", keys_names)], A1 = keys_names[grepl("^A1.*", keys_names)], A2 = keys_names[grepl("^A2.*", keys_names)], A3 = keys_names[grepl("^A3.*", keys_names)], A4 = keys_names[grepl("^A4.*", keys_names)], A5 = keys_names[grepl("^A5.*", keys_names)], A6 = keys_names[grepl("^A6.*", keys_names)], C1 = keys_names[grepl("^C1.*", keys_names)], C2 = keys_names[grepl("^C2.*", keys_names)], C3 = keys_names[grepl("^C3.*", keys_names)], C4 = keys_names[grepl("^C4.*", keys_names)], C5 = keys_names[grepl("^C5.*", keys_names)], C6 = keys_names[grepl("^C6.*", keys_names)])
keys_facet <- c("N1", "N2", "N3", "N4", "N5", "N6", "E1", "E2", "E3", "E4", "E5", "E6", "O1", "O2", "O3", "O4", "O5", "O6", "A1", "A2", "A3", "A4", "A5", "A6", "C1", "C2", "C3", "C4", "C5", "C6") 

#p-factor
par(mfrow = c(1,1), mar=c(4,4,2,0), oma=c(3,1,0,0))
manhattan(items, data.frame(P_Factor = adata$p_factor), keys = keys, labels = keys_facet, abs=F, adjust="holm", ylim = c(-.35, .35), ylab = "Correlation with", axis=F, main="Item Level Personality Profile of") #repeat for all factors




