
load("Data/CORE_intersim_data.Rdata")
source("Script/Simulated study/1.Simulated function.R")
pathDat <- R.utils::Arguments$getWritablePath("Data/SimDat")
pathMeth <- R.utils::Arguments$getWritablePath("Output/SimDat")

################################# 0.parameters setting #############################
S <- 50

n.sample <- c(rep(150, 10), 100, 300)
cluster.sample.prop = list(B1 = c(0.2, 0.3, 0.5), B2 = c(0.2, 0.3, 0.5), B3 = c(0.2, 0.3, 0.5), 
                           B4 = c(0.2, 0.3, 0.5), B5 = c(0.2, 0.3, 0.5), B6 = c(0.2, 0.3, 0.5), 
                           B7 = c(0.2, 0.3, 0.5), B8 = c(1/3, 1/3, 1/3), B9 = c(0.4, 0.6), 
                           B10 = c(0.2, 0.2, 0.3, 0.3), B11 = c(0.2, 0.3, 0.5), B12 = c(0.2, 0.3, 0.5))
n.clust <- sapply(1:length(n.sample), function(x) length(cluster.sample.prop[[x]]))
p.DMP <- c(rep(0.1, 5), 0.2, 0.05, rep(0.1, 5))
delta.methyl <- c(rep(2, 3), 1.5, 2.5, rep(2, 7))
delta.expr <- c(rep(2, 3), 1.5, 2.5, rep(2, 7))
delta.protein <- c(rep(2, 3), 1.5, 2.5, rep(2, 7))
noise.methyl <- c(1, 0.5, 2, rep(1, 9))
noise.expr <- c(1, 0.5, 2, rep(1, 9))
noise.protein <- c(1, 0.5, 2, rep(1, 9))


################################# 1.Simulated data generation #############################
set.seed(0000)

lapply(1:S, function(ss) {
  pathDat_sim <- lapply(1:length(n.sample), function(bb) R.utils::Arguments$getWritablePath(sprintf("%s/Benchmark%s", pathDat, bb)))
  file <- file.path(pathDat_sim, sprintf("simu%s.rds", ss))
  print(file)
  l <- list(n.sample, cluster.sample.prop , p.DMP, delta.methyl, delta.expr, delta.protein, noise.methyl, noise.expr, noise.protein, file)
  t <- purrr::pmap(l, CORE_intersim)
})


################################# 2.Integration analysis #############################

library(CrIMMix)
library(parallel)
listBenchmark <- gtools::mixedsort(list.files(pathDat))
nbCPU <- 20
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
  LRAcluster_results <- lapply(data, IntMultiOmics, method = "LRAcluster", K = K,
                               type = c("gaussian", "gaussian", "gaussian"))
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

