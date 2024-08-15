#setup----
library(tidyverse)
library(lavaan)
library(psych)
library(careless)
library(EFAtools)
library(glmnet)
library(bootnet)
library(qgraph)
library(sjPlot)
library(matrixStats)

adata_full = read.csv("AAAData (1).csv")

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

## restrict age to 18 plus 
adata = adata_full %>% filter(agegroup == 2 & age > 17) %>% select(-contains("SS_"))

## median replace
adata[,210:449] = adata %>% select(210:449) %>% 
  mutate_if(is.numeric, list(~ replace(., is.na(.), median(., na.rm = TRUE))))

## reverse personality items, calculate alpha
p = select(adata, 210:449) 

p[grep("P",names(p))] = 4 - p[grep("P",names(p))]

alphas = names(p) %>% strtrim(2) %>% unique %>% as.list %>%
  map(~psych::alpha(p %>% select(grep(.x, names(p)))))

map_dbl(alphas, function(x) x[["total"]][1] %>% unlist)

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

#EFA----
#subset DSM data
data <- adata[,c(27:46, 451:453)]

#test suitablity for factor analysis
KMO(data)
BARTLETT(data)

#select number of factors
PARALLEL(data, cor_method = "spearman")
VSS(data)
eigenvalues<-eigen(cor(data))$values
plot(eigenvalues, type = 'b', pch = 16,
     main = "Scree Plot", xlab="", ylab="Eigenvalues")
axis(1, at = 1:10, labels = 1:10)

#build bifactor model with 3 specific factors
efa_dsm <- omega(data, fm="minres", nfactors = 3, poly = F)

#extract factors scores for further analysis
efa_scores <- as.data.frame(factor.scores(data,efa_dsm$schmid$sl)$scores)

data$p_factor <- efa_scores$g
data$intern <- efa_scores$`F1*`
data$psych <- efa_scores$`F2*`
data$substance <- efa_scores$`F3*`


#predict p_factor using personality----
adata <- cbind(data, domains, facets, items)

table(adata$age %>% round(0))


Data_all <- cbind(data, domains, facets, items)
model <- list()
coef <- list()

for (i in 1:100){
 model[[i]]<-cv.glmnet(as.matrix(Data_all[,c(28:32)]), Data_all$substance, alpha = .05)
 coef[[i]] <- round(coef(model[[i]], s = "lambda.min"),3)
} #repeat for all factors across all levels of the personality hierarchy


coef <- lapply(coef, as.matrix)
coef <- do.call("cbind", coef)
estimate <- as.data.frame(rowMeans(coef))
sd <- as.data.frame(rowSds(coef))
nuance <- cbind(estimate, sd)
nuance <- nuance %>% 
  arrange(desc(abs(rowMeans(estimate)))) %>% round(.,3)
domain_substance <- as.data.frame(cbind(nuance, CIlower = nuance$`rowMeans(coef)`-nuance$`rowSds(coef)`, Est = nuance$`rowMeans(coef)`, CIupper = nuance$`rowMeans(coef)` + nuance$`rowSds(coef)`))[,3:5]

#create and validate elastic net prediction models----
foo1 = function(x1,y,x2=NULL,p=.75) {
  s = sample(nrow(x1), p*nrow(x1))
  if(is.null(x2)) x2 = x1
  x1 = sapply(x1,scale)
  x2 = sapply(x2,scale)
  cv.glmnet(x1[s,], y[s], alpha = .05) %>% 
    predict(x2[-s,], s = "lambda.min") %>%  
    cor(y[-s])
}

# Run the models through domains, facets and items 
# Replicate 100 times


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




