
load("Data/CORE_intersim_data.Rdata")
source("Script/Simulated study/1.Simulated function.R")
pathDat <- R.utils::Arguments$getWritablePath("Data/SimDat_Percent")
pathMeth <- R.utils::Arguments$getWritablePath("Output/SimDat_Percent")

S <- 50

p.DMP <- c(0.01, seq(0.05, 0.25, by = 0.05))
n <- length(p.DMP)

n.sample <- c(rep(150, n))
cluster.sample.prop = list(B1 = c(0.2, 0.3, 0.5), B2 = c(0.2, 0.3, 0.5), B3 = c(0.2, 0.3, 0.5), 
                           B4 = c(0.2, 0.3, 0.5), B5 = c(0.2, 0.3, 0.5), B6 = c(0.2, 0.3, 0.5))
n.clust <- sapply(1:n, function(x) length(cluster.sample.prop[[x]]))

delta.methyl <- c(rep(2, n))
delta.expr <- c(rep(2, n))
delta.protein <- c(rep(2, n))
noise.methyl <- c(rep(FALSE, n))
noise.expr <- c(rep(FALSE, n))
noise.protein <- c(rep(FALSE, n))


set.seed(1111)

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
  adjustedRIs$benchmark = b
 
  return(adjustedRIs)
}))

ARI_dat$method <- gsub("NMF", "intNMF", ARI_dat$method)
ARI_dat$method <- gsub("iCluster", "iClusterPlus", ARI_dat$method)
ARI_dat$method <- ARI_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                       "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
ARI_dat$benchmark <- ARI_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:10))
ARI_dat_mean <- ARI_dat %>% 
  group_by(benchmark, method) %>% summarise(mean = mean(ARI)) %>% data.frame() %>% 
  mutate(Index = "ARI")

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
  NormalizedMIs$benchmark = b
 
  return(NormalizedMIs)
}))

NMI_dat$method <- gsub("NMF", "intNMF", NMI_dat$method)
NMI_dat$method <- gsub("iCluster", "iClusterPlus", NMI_dat$method)
NMI_dat$method <- NMI_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                       "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
NMI_dat$benchmark <- NMI_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:10))
NMI_dat_mean <- NMI_dat %>% 
  group_by(benchmark, method) %>% summarise(mean = mean(NMI)) %>% data.frame() %>% 
  mutate(Index = "NMI")

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
 
  return(F_measures)
}))

F1_dat$method <- gsub("NMF", "intNMF", F1_dat$method)
F1_dat$method <- gsub("iCluster", "iClusterPlus", F1_dat$method)
F1_dat$method <- F1_dat$method %>% factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                     "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))
F1_dat$benchmark <- F1_dat$benchmark %>% factor(levels = paste0("Benchmark", 1:10))
F1_dat_mean <- F1_dat %>% 
  group_by(benchmark, method) %>% summarise(mean = mean(F1)) %>% data.frame() %>% 
  mutate(Index = "F1-score")


Sim_percent_res <- rbind(ARI_dat_mean, NMI_dat_mean, F1_dat_mean)
Sim_percent_res$Index <- factor(Sim_percent_res$Index, levels = c("ARI", "NMI", "F1-score"))

Sim_percent_res_wide <- Sim_percent_res %>% mutate(mean = sprintf("%.3f", mean) %>% as.numeric)
Sim_percent_res_wide <- Sim_percent_res_wide %>% pivot_wider(names_from = method, values_from = mean)

library(ggplot2)
library(ggalt)

g_percent_line <- Sim_percent_res %>% 
  ggplot(aes(x = benchmark, y = mean, group = method, color = method, shape = method)) +
  geom_point(size = 2.5) + 
  geom_line(linewidth = 1) + 
  scale_color_manual(values  = MetBrewer::met.brewer("VanGogh2", 8)) +
  scale_shape_manual(values = c(1:8)) + 
  facet_wrap(. ~ Index, ncol = 3, labeller = label_wrap_gen(multi_line = T)) +
  theme_bw() + xlab("Proportion") + ylab("") +
  scale_x_discrete(labels = as.character(paste0(p.DMP*100, "%"))) + 
  theme(panel.border = element_rect(linewidth = 1),
        strip.text.x = element_text(size = 18), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 18), 
        legend.position = "top",
        legend.spacing.x = unit(0.5, 'cm'),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

