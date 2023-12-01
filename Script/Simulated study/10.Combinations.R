
pathDat <- R.utils::Arguments$getWritablePath("Data/SimDat")
pathMeth <- R.utils::Arguments$getWritablePath("Output/SimDat")

methods <- c("LRAcluster", "PINSPlus", "SNF", "CIMLR", "MCIA", "Mocluster", "NMF", "iCluster")
listBenchmark <- gtools::mixedsort(list.files(pathDat))
nbCPU <- 50
ii <- 1
b <- listBenchmark[ii]
K <- n.clust[ii]
N <- n.sample[ii]

pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
pathMeth_comb <- R.utils::Arguments$getWritablePath("Output/Comb")
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth_comb, b))
print(pathMeth_sub)

######### Benchmark1
list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
data <- lapply(list.sim, function (ll) ll$data)

comb <- list(c(1:3), 1:2, 2:3, c(1,3))
time_start <- Sys.time()

lapply(1:length(data), function (rr){
  data1 <- data[[rr]]
  
  fit_lambda <- tune.iClusterPlus(cpus = 15, dt1 = data1$dat.methyl, dt2 = data1$dat.expr, dt3 = data1$dat.protein,
                                  type = c("gaussian","gaussian","gaussian"), K = K-1, scale.lambda = c(1,1,1))
  nLambda <- nrow(fit_lambda$lambda)
  nK <- length(fit_lambda)
  BIC <- do.call(c, lapply(1:nLambda, function(nn) {fit_lambda[["fit"]][[nn]][["BIC"]]}))
  para_lambda <- fit_lambda$lambda[which.min(BIC), ]
  
  res <- lapply(1:length(comb), function (cc){
    c <- comb[[cc]]
    
    print("LRAcluster")
    LRAcluster_results <- IntMultiOmics(data1[c],  method = "LRAcluster", K = K,
                                        type = c("gaussian", "gaussian", "gaussian")[c])
    saveRDS(LRAcluster_results, file = file.path(pathMeth_sub, sprintf("LRAcluster_res_comb%s.rds", cc)))

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

    PINSPlus_results <-  doPINSPlus_V2(data1[c], K = K)
    saveRDS(PINSPlus_results, file = file.path(pathMeth_sub, sprintf("PINSPlus_res_comb%s.rds", cc)))

    print("SNF")
    SNFresults <- IntMultiOmics(data1[c], method = "SNF", K = K, K_n = N/10)
    saveRDS(SNFresults, file = file.path(pathMeth_sub, sprintf("SNF_res_comb%s.rds", cc)))

    print("CIMLR")
    CIMLR_results <-  IntMultiOmics(data1[c],  method = "CIMLR", K = K)
    saveRDS(CIMLR_results, file = file.path(pathMeth_sub, sprintf("CIMLR_res_comb%s.rds", cc)))

    print("MCIA")
    MCIAresults <- IntMultiOmics(data1[c],  method = "MCIA", ncomp = 2, K = K)
    saveRDS(MCIAresults, file = file.path(pathMeth_sub, sprintf("MCIA_res_comb%s.rds", cc)))

    print("Mocluster")
    Moaresults <-  IntMultiOmics(data1[c],  method = "Mocluster", K = K, ncomp = 4,
                                 k = c(0.2 ,0.4, 0.2)[c])
    saveRDS(Moaresults, file = file.path(pathMeth_sub, sprintf("Mocluster_res_comb%s.rds", cc)))

    print("NMF")
    NMFresults <- IntMultiOmics(data1[c],  method = "intNMF", K = K)
    saveRDS(NMFresults, file = file.path(pathMeth_sub, sprintf("NMF_res_comb%s.rds", cc)))

    print("iCluster")
    iCluster_results <-  IntMultiOmics(data1[c],  method = "iCluster", K = K-1, lambda = para_lambda[c],
                                       type = c("gaussian", "gaussian", "gaussian")[c])
    saveRDS(iCluster_results, file = file.path(pathMeth_sub, sprintf("iCluster_res_comb%s.rds", cc)))
    
    time_end <- Sys.time()
    print(time_end)

    res_tot <- list(LRAcluster_results$clust, PINSPlus_results$clust, SNFresults$clust, CIMLR_results$clust, 
                    MCIAresults$clust, Moaresults$clust, NMFresults$clust, iCluster_results$clust)
    
    names(res_tot) <- c("LRAcluster", "PINSPlus", "SNF", "CIMLR", "MCIA", "Mocluster", "NMF", "iCluster")
    return(res_tot)
  })
  names(res) <- sprintf("com%s", 1:4)

  saveRDS(res, file = file.path(pathMeth_sub, sprintf("res_sim%s.rds", rr)))
})


true.clusters <- list()

pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
print(pathMeth_sub)
list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
true.clusters <- lapply(list.sim, function(ii) ii[["clustering.assignment"]]$cluster.id)

S <- 50
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth_comb, b))
print(pathMeth_sub)

ii <- 1

### ARI
ARI_dat_comb <- do.call(rbind, lapply(1:S, function(ii){
  files <- list.files(pathMeth_sub, pattern = sprintf("res_sim%s.rds", ii), full.names = TRUE)
  fil1 <- readRDS(files)
  dat <- do.call(cbind, lapply(fil1, function (r){
    adjRI <- r  %>% sapply(mclust::adjustedRandIndex, true.clusters[[ii]])
  }))
  dat <- dat %>% data.frame %>% mutate(method = rownames(dat), sim = sprintf("sim%s", ii), index = rep('ARI', length(methods))) %>% 
    gather(com1:com4, value = value, key = comb)
  return(dat)
}))

### NMI
NMI_dat_comb <- do.call(rbind, lapply(1:S, function(ii){
  files <- list.files(pathMeth_sub, pattern = sprintf("res_sim%s.rds", ii), full.names = TRUE)
  fil1 <- readRDS(files)
  dat <- do.call(cbind, lapply(fil1, function (r){
    NormalMI <- r  %>% sapply(aricode::NMI, true.clusters[[ii]])
  }))
  dat <- dat %>% data.frame %>% mutate(method = rownames(dat), sim = sprintf("sim%s", ii), index = rep('NMI', length(methods))) %>% 
    gather(com1:com4, value = value, key = comb)
  return(dat)
}))

### F1
F1_dat_comb <- do.call(rbind, lapply(1:S, function(ii){
  files <- list.files(pathMeth_sub, pattern = sprintf("res_sim%s.rds", ii), full.names = TRUE)
  fil1 <- readRDS(files)
  dat <- do.call(cbind, lapply(fil1, function (r){
    Fscore <- r  %>% sapply(FlowSOM::FMeasure, true.clusters[[ii]])
  }))
  dat <- dat %>% data.frame %>% mutate(method = rownames(dat), sim = sprintf("sim%s", ii), index = rep('F1', length(methods))) %>% 
    gather(com1:com4, value = value, key = comb)
  return(dat)
}))
  
dat_omics_comb <- rbind(ARI_dat_comb, NMI_dat_comb, F1_dat_comb)
  
comb_list <- list(c("M", "G", "P"), c("M", "G"), c("G", "P"), c("M", "P"))
dat_omics_comb <- dat_omics_comb %>% mutate(comb = as.numeric(as.factor(comb)))
combination <- sapply(dat_omics_comb$comb, function (cc) paste(comb_list[[cc]], collapse = "+"))
dat_omics_comb <-  dat_omics_comb %>% mutate(combination = factor(combination, levels = sapply(comb_list, paste, collapse = "+")))

dat_omics_comb$method <- gsub("NMF", "intNMF", dat_omics_comb$method)
dat_omics_comb$method <- gsub("iCluster", "iClusterPlus", dat_omics_comb$method)
dat_omics_comb$method <- dat_omics_comb$method %>% 
  factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster", "intNMF", "MCIA", "iClusterPlus", "PINSPlus"))

g_omics_comb <- dat_omics_comb %>% ggplot(aes(y = value, x = index, fill = combination)) +  
  geom_boxplot() + 
  facet_wrap(. ~ method, ncol = 4) + 
  scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2", 4)) + 
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1),
        strip.text.x = element_text(size = 18), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 18), 
        legend.position = "top",
        legend.spacing.x = unit(0.5, 'cm'),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
g_omics_comb


library(ggradar)
library(tibble)

omics_comb_radar <- dat_omics_comb %>% group_by(method, combination, index) %>% summarise(value_mean = mean(value))
omics_comb_method <- omics_comb_radar %>% spread(key = combination, value = value_mean)

levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster", "intNMF", "MCIA", "iClusterPlus", "PINSPlus")

radarplot <- lapply(1:length(methods), function(mm) {
  ggradar(data.frame(omics_comb_method[which(omics_comb_method == levels[mm]), 2:6], check.names = F),
          values.radar = c("0", "0.5", "1"),
          background.circle.colour = "white",
          grid.min = 0,
          grid.mid = 0.5,
          grid.max = 1,
          grid.line.width = 1,
          gridline.min.linetype = 4,
          gridline.mid.linetype = 4,
          gridline.max.linetype = 4,
          gridline.mid.colour = "purple",
          grid.label.size = 6,
          group.point.size = 6,
          group.line.width = 1,
          plot.title = levels[mm], 
          legend.position = "bottom")
})

library(patchwork)

omics_radarplot <- radarplot[[1]] + radarplot[[2]] + radarplot[[3]] + radarplot[[4]] | 
  radarplot[[5]] + radarplot[[6]] + radarplot[[7]] + radarplot[[8]]
omics_radarplot
