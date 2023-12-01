
load("Data/CORE_intersim_data.Rdata")
source("Script/Simulated study/1.Simulated function.R")
pathDat <- R.utils::Arguments$getWritablePath("Data/SimDat_Proportion")
pathMeth <- R.utils::Arguments$getWritablePath("Output/SimDat_Proportion")

S <- 50

n.sample <- c(rep(150, 6))
cluster.sample.prop = list(B1 = c(0.1, 0.3, 0.6), B2 = c(0.2, 0.3, 0.5), B3 = c(1/3, 1/3, 1/3), 
                           B4 = c(0.1, 0.1, 0.2, 0.6), B5 = c(0.1, 0.2, 0.3, 0.4), B5 = c(1/4, 1/4, 1/4, 1/4))
n.clust <- sapply(1:length(n.sample), function(x) length(cluster.sample.prop[[x]]))
p.DMP <- c(rep(0.1, 6))
delta.methyl <- c(rep(2, 6))
delta.expr <- c(rep(2, 6))
delta.protein <- c(rep(2, 6))
noise.methyl <- c(rep(FALSE, 6))
noise.expr <- c(rep(FALSE, 6))
noise.protein <- c(rep(FALSE, 6))

set.seed(6666)

lapply(1:S, function(ss) {
  pathDat_sim <- lapply(1:length(n.sample), function(bb) R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, bb)))
  file <- file.path(pathDat_sim, sprintf("simu%s.rds", ss))
  print(file)
  l <- list(n.sample, cluster.sample.prop , p.DMP, delta.methyl, delta.expr, delta.protein, noise.methyl, noise.expr, noise.protein, file)
  t <- purrr::pmap(l, CORE_intersim)
})

library(CrIMMix)
library(parallel)
listBenchmark <- gtools::mixedsort(list.files(pathDat))
nbCPU <- 15
time_start <- Sys.time()

for(ii in seq(along = listBenchmark)){
  b <- listBenchmark[ii]
  K <- n.clust[ii]
  N <- n.sample[ii]
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  data <- lapply(list.sim, function (ll) ll$data)
  
  ### Early integration
  print("LRAcluster")
  LRAcluster_results <- lapply(data, IntMultiOmics, method = "LRAcluster", K = K, type = c("gaussian", "gaussian", "gaussian"))
  saveRDS(LRAcluster_results, file = file.path(pathMeth_sub, sprintf("LRAcluster_res.rds")))
  
  ### Late integration PINSPlus
  print("PINSPlus")
  doPINSPlus_V2 <- function (data, K) {
    require(PINSPlus)
    fit <- SubtypingOmicsData(data, k = K)
    fit$cluster2 <- as.integer(factor(fit$cluster2, levels = unique(fit$cluster2), 
                                      labels = 1:length(unique(fit$cluster2))))
    names(fit$cluster2) <- names(fit$cluster1)
    res <- list(clust = fit$cluster2, fit = fit)
    return(res)
  }
  
  PINSPlusresults <- mclapply(data, doPINSPlus_V2, K = K, mc.cores = nbCPU)
  saveRDS(PINSPlusresults, file = file.path(pathMeth_sub, sprintf("PINSPlus_res.rds")))
  
  ### Similarity network
  print("SNF")
  SNFresults <- mclapply(data, IntMultiOmics, method = "SNF", K = K, K_n = N/10, mc.cores = nbCPU)
  saveRDS(SNFresults, file = file.path(pathMeth_sub, sprintf("SNF_res.rds")))
  
  print("CIMLR")
  CIMLR_results <- lapply(data, IntMultiOmics, method = "CIMLR", K = K)
  saveRDS(CIMLR_results, file = file.path(pathMeth_sub, sprintf("CIMLR_res.rds")))
  
  ### Dimensional reduction
  print("MCIA")
  MCIAresults <- mclapply(data, IntMultiOmics, method = "MCIA", K = K, ncomp = 2, mc.cores = nbCPU)
  saveRDS(MCIAresults, file = file.path(pathMeth_sub, sprintf("MCIA_res.rds")))
  
  print("Mocluster")
  Moaresults <- mclapply(data, IntMultiOmics, method = "Mocluster", K = K, ncomp = 4,
                         k = c(0.2, 0.4, 0.2), mc.cores = nbCPU)
  saveRDS(Moaresults, file = file.path(pathMeth_sub, sprintf("Mocluster_res.rds")))
  
  print("NMF")
  NMFresults <- mclapply(data, IntMultiOmics, method = "intNMF", K = K, mc.cores = nbCPU)
  saveRDS(NMFresults, file = file.path(pathMeth_sub, sprintf("NMF_res.rds")))
  
  ### Statistical methods
  print("iCluster")
  
  para_lambda <- do.call(rbind, lapply(data, function (dd) {
    fit_lambda <- tune.iClusterPlus(cpus = 12, dt1 = dd$dat.methyl, dt2 = dd$dat.expr, dt3 = dd$dat.protein,
                                    type = c("gaussian","gaussian","gaussian"), K = K-1, scale.lambda = c(1,1,1))
    nLambda = nrow(fit_lambda$lambda)
    nK = length(fit_lambda)
    BIC <- do.call(c, lapply(1:nLambda, function(nn){fit_lambda[["fit"]][[nn]][["BIC"]]}))
    lambda <- fit_lambda$lambda[which.min(BIC), ]
  }))
  
  iCluster_results <- mclapply(1:S, function (ss) {
    IntMultiOmics(data[[ss]], method = "iCluster", K = K-1, lambda = para_lambda[ss, ],
                  type = c("gaussian", "gaussian", "gaussian"))
  })
  saveRDS(iCluster_results, file = file.path(pathMeth_sub, sprintf("iCluster_res.rds")))
  
  time_end <- Sys.time()
  print(time_end)
}


library(tidyverse)
methods <- c("LRAcluster", "PINSPlus", "SNF", "CIMLR", "MCIA", "Mocluster", "NMF", "iCluster")

listBenchmark <- gtools::mixedsort(list.files(pathDat))

true.clusters <- list()
for (ii in seq(along = listBenchmark)) {
  b <- listBenchmark[ii]
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  true.clusters[[ii]] <- lapply(list.sim, function(ii) ii[["clustering.assignment"]]$cluster.id)
}

########## Computing ARI

No.clust <- paste0("Clust = ", n.clust)
Inequality <- rep(c("Extreme", "Moderate", "Equality"), 2)

ARI_dat <- do.call(rbind, lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  
  adjustedRIs <- do.call(rbind, lapply(methods, function(mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    
    adjRI <- sapply(1:S, function(sss) {r[[sss]]$clust %>% mclust::adjustedRandIndex(true.clusters[[ii]][[sss]])})
    
    data.frame(method = mm, ARI = adjRI)
  }))
  adjustedRIs$benchmark <- b
  adjustedRIs$clust <- No.clust[ii]
  adjustedRIs$equality <- Inequality[ii]
  
  return(adjustedRIs)
}))

ARI_dat$method <- gsub("NMF", "intNMF", ARI_dat$method)
ARI_dat$method <- gsub("iCluster", "iClusterPlus", ARI_dat$method)
ARI_dat$method <- ARI_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                       "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
ARI_dat$benchmark <- ARI_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:6))
ARI_dat$clust <- ARI_dat$clust %>% factor(levels = c("Clust = 3", "Clust = 4"))
ARI_dat$equality <- ARI_dat$equality %>% factor(levels = c("Extreme", "Moderate", "Equality"))

ARI_dat_mean <- ARI_dat %>% group_by(clust, equality, method) %>% 
  summarise(mean = mean(ARI)) %>% data.frame() %>% mutate(Index = "ARI")

##########  Computing NMI

NMI_dat <- do.call(rbind, lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  
  NormalizedMIs <- do.call(rbind, lapply(methods, function(mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    
    NMI <- sapply(1:S, function(sss) {r[[sss]]$clust %>% aricode::NMI(true.clusters[[ii]][[sss]])})
    
    data.frame(method = mm, NMI = NMI)
  }))
  NormalizedMIs$benchmark <- b
  NormalizedMIs$clust <- No.clust[ii]
  NormalizedMIs$equality <- Inequality[ii]
  
  return(NormalizedMIs)
}))

NMI_dat$method <- gsub("NMF", "intNMF", NMI_dat$method)
NMI_dat$method <- gsub("iCluster", "iClusterPlus", NMI_dat$method)
NMI_dat$method <- NMI_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                       "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
NMI_dat$benchmark <- NMI_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:6))
NMI_dat$clust <- NMI_dat$clust %>% factor(levels = c("Clust = 3", "Clust = 4"))
NMI_dat$equality <- NMI_dat$equality %>% factor(levels = c("Extreme", "Moderate", "Equality"))

NMI_dat_mean <- NMI_dat %>% group_by(clust, equality, method) %>% 
  summarise(mean = mean(NMI)) %>% data.frame() %>% mutate(Index = "NMI")

##########  Computing F1-score

F1_dat <- do.call(rbind, lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  
  F_measures <- do.call(rbind, lapply(methods, function(mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    
    F1_score <- sapply(1:S, function(sss) {r[[sss]]$clust %>% FlowSOM::FMeasure(true.clusters[[ii]][[sss]])})
    
    data.frame(method = mm, F1 = F1_score)
  }))
  F_measures$benchmark = b
  F_measures$clust <- No.clust[ii]
  F_measures$equality <- Inequality[ii]
  
  return(F_measures)
}))

F1_dat$method <- gsub("NMF", "intNMF", F1_dat$method)
F1_dat$method <- gsub("iCluster", "iClusterPlus", F1_dat$method)
F1_dat$method <- F1_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                     "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
F1_dat$benchmark <- F1_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:6))
F1_dat$clust <- F1_dat$clust %>% factor(levels = c("Clust = 3", "Clust = 4"))
F1_dat$equality <- F1_dat$equality %>% factor(levels = c("Extreme", "Moderate", "Equality"))

F1_dat_mean <- F1_dat %>% group_by(clust, equality, method) %>% 
  summarise(mean = mean(F1)) %>% data.frame() %>% mutate(Index = "F1-score")

library(ggplot2)
library(ggalt)

Equal_data <- rbind(ARI_dat_mean, NMI_dat_mean, F1_dat_mean) %>% 
  mutate(Index = factor(Index, levels = c("ARI", "NMI", "F1-score")))

Equal_data_wide <- Equal_data %>% mutate(mean = sprintf("%.3f", mean) %>% as.numeric)
Equal_data_wide <- Equal_data_wide %>% pivot_wider(names_from = method, values_from = mean)

g_equal_line <- Equal_data %>% ggplot(aes(x = equality, y = mean, group = method, color = method, shape = method)) +
  geom_point(size = 2.5) + 
  geom_line(linewidth = 1) + 
  scale_color_manual(values  = MetBrewer::met.brewer("VanGogh2", 8)) +
  scale_shape_manual(values = c(1:8)) + 
  facet_grid(clust ~ Index, labeller = label_wrap_gen(multi_line = T)) +
  theme_bw() + xlab("") + ylab("") + 
  theme(panel.border = element_rect(linewidth = 1),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18),
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 18), 
        legend.position = "top",
        legend.spacing.x = unit(0.5, 'cm'),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
