
rm(list = ls())

library(tidyverse)
library(survival)
library(survminer)
library(CrIMMix)
library(iClusterPlus)
library(NbClust)
library(car)

load("Data/Mergedata.Rdata")

Mergedata$Methy <- logit(Mergedata$Methy) %>% 
  t %>% scale 

Mergedata$FPKM <- Mergedata$FPKM %>% 
  t %>% scale 

Mergedata$RPPA <- Mergedata$RPPA %>% t 

mad <- lapply(Mergedata[-1], function(dd) {
  apply(dd, 2, mad)
})

Clinical <- Mergedata$Clinical
CRC_mad <- Mergedata[-1]

CRC_mad$Methy <- CRC_mad$Methy[, order(mad$Methy, decreasing = T)[1:1000]]
CRC_mad$FPKM <- CRC_mad$FPKM[, order(mad$FPKM, decreasing = T)[1:1000]]
CRC_mad$RPPA <- CRC_mad$RPPA[, order(mad$RPPA, decreasing = T)[1:100]]

pathOutput <- R.utils::Arguments$getWritablePath("Output/CRC")
print(pathOutput)

##### LRAcluster

print("LRAcluster")
LRAcluster_results <- IntMultiOmics(CRC_mad, method = "LRAcluster", K = 3, 
                                    type = c("gaussian", "gaussian", "gaussian"))

scr <- LRAcluster_results$fit$coordinate %>% t
res.nbclust <- NbClust(scr, distance = "euclidean",
                       min.nc = 2, max.nc = 6, 
                       method = "ward.D2", index ="all")
K_LRA <- res.nbclust$Best.nc[1, ] %>% table %>% which.max %>% names %>% as.numeric()


LRAcluster_results <- IntMultiOmics(CRC_mad, method = "LRAcluster", K = K_LRA, 
                                    type = c("gaussian", "gaussian", "gaussian"))
saveRDS(LRAcluster_results, file = file.path(pathOutput, sprintf("LRAcluster_res.rds")))

##### PINSPlus

print("PINSPlus")
PINSPlus_results <- IntMultiOmics(CRC_mad, method = "PINSPlus", K = 6)
K_PINSPlus <- PINSPlus_results$clust %>% table %>% length %>% as.numeric()

doPINSPlus_V2 <- function (data, K) {
  require(PINSPlus)
  fit <- SubtypingOmicsData(data, k = K)
  fit$cluster2 <- as.integer(factor(fit$cluster2, levels = unique(fit$cluster2), 
                                    labels = 1:length(unique(fit$cluster2))))
  names(fit$cluster2) <- names(fit$cluster1)
  res <- list(clust = fit$cluster2, fit = fit)
  return(res)
}

PINSPlus_results <-  doPINSPlus_V2(CRC_mad, K = K_PINSPlus)
saveRDS(PINSPlus_results, file = file.path(pathOutput, sprintf("PINSPlus_res.rds")))

##### SNF

print("SNF")
SNF <- IntMultiOmics(CRC_mad, method = "SNF", K = 3, K_n = nrow(CRC_mad[[1]])/10)
estimationResult <- SNFtool::estimateNumberOfClustersGivenGraph(SNF$fit, 2:6) %>% unlist %>% table
K_SNF <- estimationResult %>% which.max %>% names() %>% as.numeric()

SNFresults <- IntMultiOmics(CRC_mad, method = "SNF", K = K_SNF, K_n = nrow(CRC_mad[[1]])/10)
saveRDS(SNFresults, file = file.path(pathOutput, sprintf("SNF_res.rds")))

### CIMLR

print("CIMLR")
NUMC <- 2:6
CIMLR_results <- lapply(CRC_mad, t) %>% 
  CIMLR::CIMLR_Estimate_Number_of_Clusters(NUMC = NUMC, cores.ratio = 0)
K_CIMLR <- NUMC[which.min(CIMLR_results$K1)] %>% as.numeric()

CIMLR_results <-  IntMultiOmics(CRC_mad,  method = "CIMLR", K = K_CIMLR)
saveRDS(CIMLR_results, file = file.path(pathOutput, sprintf("CIMLR_res.rds")))

##### MCIA

print("MCIA")
MCIAresults <- IntMultiOmics(CRC_mad, method = "MCIA", ncomp = 2, K = 3)
res.nbclust <- NbClust(MCIAresults$fit$mcoa$SynVar, distance = "euclidean",
                       min.nc = 2, max.nc = 6, 
                       method = "ward.D2", index ="all")
K_MCIA <- res.nbclust$Best.nc[1, ] %>% table %>% which.max %>% names %>% as.numeric()

MCIAresults <- IntMultiOmics(CRC_mad, method = "MCIA", ncomp = 2, K = K_MCIA)
saveRDS(MCIAresults, file = file.path(pathOutput, sprintf("MCIA_res.rds")))

##### Mocluster
print("Mocluster")

Moaresults <- IntMultiOmics(CRC_mad, method = "Mocluster", K = 3, k = c(0.1, 0.1, 0.1), ncomp = 4)
scr <- mogsa::moaScore(Moaresults$fit)
res.nbclust <- NbClust(scr, distance = "euclidean",
                       min.nc = 2, max.nc = 6, 
                       method = "ward.D2", index ="all")
K_Moa <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names %>% as.numeric()


c_1 <- (10/sqrt((CRC_mad %>% sapply(dim))[2,])) %>% round(2)
SGCCAresults <- IntMultiOmics(CRC_mad, method = "SGCCA", K = K_Moa, c1 = c_1)
a <- SGCCAresults$fit$a
posDat <- lapply(a, function(aa) which(aa %>% rowSums != 0))
selectVars_sgcca <- lapply(posDat, names)
k <- selectVars_sgcca %>% sapply(length)/(CRC_mad %>% sapply(dim))[2,]

Moaresults <- IntMultiOmics(CRC_mad, method = "Mocluster", K = K_Moa, k = k, ncomp = 4)
saveRDS(Moaresults, file = file.path(pathOutput, sprintf("Mocluster_res.rds")))

##### intNMF

print("intNMF")

dat <- lapply(CRC_mad, function (dd){
  if (!all(dd>=0)) dd <- pmax(dd + abs(min(dd)), .Machine$double.eps)
  dd <- dd/max(dd)
  return(dd %>% as.matrix)
})

NMFresults <- IntNMF::nmf.opt.k(dat, k.range = 2:6, n.runs = 5, progress = TRUE)
K_intNMF <- gsub("k", "", NMFresults %>% rowMeans %>% which.max %>% names %>% as.numeric())

NMFresults <- IntMultiOmics(CRC_mad, method = "intNMF", K = K_intNMF)
saveRDS(NMFresults, file = file.path(pathOutput, sprintf("NMF_res.rds")))

##### iClusterPlus

print("iClusterPlus")

icluster_best_k <- lapply(1:5, function(kk){
  cv.fit <- tune.iClusterPlus(cpus = 15, 
                              dt1 = as.matrix(CRC_mad[[1]]), 
                              dt2 = as.matrix(CRC_mad[[2]]), 
                              dt3 = as.matrix(CRC_mad[[3]]), 
                              type = c("gaussian", "gaussian", "gaussian"), 
                              K = 2, scale.lambda = c(1, 1, 1))
})
saveRDS(icluster_best_k, file = file.path(pathOutput, sprintf("icluster_best_k.rds")))
icluster_kbest <- readRDS(file.path(pathOutput, sprintf("icluster_best_k.rds")))

nK <- length(icluster_kbest)
BIC <- getBIC(icluster_kbest)
devR <- getDevR(icluster_kbest)

minBICid <- apply(BIC, 2, which.min)
devRatMinBIC = rep(NA, nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i], i]
}

### Here icluster extraction
plot(1:(nK + 1), c(0, devRatMinBIC), type = "b", xlab = "Number of clusters (K+1)", ylab = "%Explained Variation")
K_icluster <- c(K_best=diff(c(0,devRatMinBIC)) %>% abs %>% which.min , method = "iClusterPlus")

if (length(K_icluster) == 0){
  K_iclusterplus <- 3
} else {
  K_iclusterplus <- K_icluster
}

fit_lambda <- tune.iClusterPlus(cpus = 40,
                                dt1 = CRC_mad[[1]],
                                dt2 = CRC_mad[[2]],
                                dt3 = CRC_mad[[3]],
                                type = c("gaussian","gaussian","gaussian"),
                                K = 2, scale.lambda = c(1, 1, 1))

nLambda <- nrow(fit_lambda$lambda)
nK <- length(fit_lambda)
BIC <- do.call(c, lapply(1:nLambda, function(nn){fit_lambda[["fit"]][[nn]][["BIC"]]}))
para_lambda <- fit_lambda$lambda[which.min(BIC), ]

iCluster_results <- IntMultiOmics(lapply(CRC_mad, as.matrix), 
                                  method = "iCluster", K = 2,
                                  #lambda = c(0.3, 0.3, 0.3), 
                                  type = c("gaussian", "gaussian", "gaussian"))
saveRDS(iCluster_results, file = file.path(pathOutput, sprintf("iCluster_res.rds")))

