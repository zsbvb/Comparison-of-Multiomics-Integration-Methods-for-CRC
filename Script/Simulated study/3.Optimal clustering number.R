
nbCPU <- 50
ii = 1 ### based on benchmark1

b <- listBenchmark[ii]
K <- n.clust[ii]
pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
print(pathMeth_sub)

list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
data <- lapply(list.sim, function (ll) ll$data)

S <- 50

Kbest <- lapply(1:S, function(ii){
  print("LRAcluster")
  LRAcluster_results <- data[[ii]] %>% IntMultiOmics(method = "LRAcluster", K = K, type = c("gaussian", "gaussian", "gaussian"))
  scr <- LRAcluster_results$fit$coordinate %>% t
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 6, 
                         method = "ward.D2", index ="all")
  K_LRA <- res.nbclust$Best.nc[1, ] %>% table %>% which.max %>% names
  
  
  print("PINSPlus")
  PINSPlus_results <- data[[ii]] %>% IntMultiOmics(method = "PINSPlus", K = 6)
  K_PINSPlus <- PINSPlus_results$clust %>% table %>% length
  
  
  print("SNF")
  SNF <- IntMultiOmics(data[[ii]], method = "SNF", K = K)
  estimationResult <- SNFtool::estimateNumberOfClustersGivenGraph(SNF$fit, 2:6) %>% unlist %>% table
  K_SNF <- estimationResult %>% which.max %>% names()
  
  
  print("CIMLR")
  NUMC <- 2:6
  CIMLR_results <- lapply(data[[ii]], t) %>% CIMLR::CIMLR_Estimate_Number_of_Clusters(NUMC = NUMC, cores.ratio = 0)
  K_CIMLR <- NUMC[which.min(CIMLR_results$K1)]
  
  
  library(NbClust)
  print("MCIA")
  MCIAresults <- data[[ii]] %>% IntMultiOmics(method = "MCIA", K = K)
  res.nbclust <- NbClust(MCIAresults$fit$mcoa$SynVar, distance = "euclidean",
                         min.nc = 2, max.nc = 6, 
                         method = "ward.D2", index ="all")
  K_MCIA <- res.nbclust$Best.nc[1, ] %>% table %>% which.max %>% names
  
  
  print("Mocluster")
  Moaresults <- data[[ii]] %>% IntMultiOmics(method = "Mocluster", K = K, ncomp = 3, k = c(0.1, 0.1, 0.1))
  scr <- mogsa::moaScore(Moaresults$fit)
  res.nbclust <- NbClust(scr, distance = "euclidean",
                         min.nc = 2, max.nc = 6, 
                         method = "ward.D2", index ="all")
  K_Moa <- res.nbclust$Best.nc[1,] %>% table %>% which.max %>% names
  
  
  print("NMF")
  dat <- lapply(data[[ii]], function (dd){
    if (!all(dd>=0)) dd <- pmax(dd + abs(min(dd)), .Machine$double.eps)
    dd <- dd/max(dd)
    return(dd %>% as.matrix)
  })
  NMFresults <- IntNMF::nmf.opt.k(dat[c(1,3)], k.range = 2:6, n.runs = 5, progress = TRUE)
  K_intNMF <- gsub("k", "", NMFresults %>% rowMeans %>% which.max %>% names)
  
  K_best <- c(K_LRA, K_PINSPlus, K_SNF, K_CIMLR, K_MCIA, K_Moa, K_intNMF)
  df <- cbind(K_best, method= c("LRAcluster", "PINSPlus", "SNF", "CIMLR", "MCIA", "Mocluster", "NMF"))
})

save(Kbest, file = "Kbest.rds")

path_icluster_k <- R.utils::Arguments$getWritablePath(paste(pathMeth_sub, "icluster_best_k", sep = "/"))
library(iClusterPlus)
lapply(1:S, function(ii) {
  input <- data[[ii]]
  icluster_best_k <- lapply(1:5, function(kk){
    cv.fit <- tune.iClusterPlus(cpus = 12, dt1 = input$dat.methyl, dt2 = input$dat.expr, dt3 = input$dat.protein, 
                                type = c("gaussian", "gaussian", "gaussian"), K = kk)
  })
  saveRDS(icluster_best_k, file = file.path(path_icluster_k, sprintf("simu%s.rds", ii)))
})


icluster_kbest <- list.files(path_icluster_k, full.names = T) %>% lapply(readRDS)

Kbest_all <- lapply(1:S, function (s){
  nK <- length(icluster_kbest[[s]])
  BIC <- getBIC(icluster_kbest[[s]])
  devR <- getDevR(icluster_kbest[[s]])
  
  minBICid <- apply(BIC, 2, which.min)
  devRatMinBIC = rep(NA, nK)
  for(i in 1:nK){
    devRatMinBIC[i] = devR[minBICid[i], i]
  }
  
  ### Here icluster extraction
  plot(1:(nK + 1), c(0, devRatMinBIC), type = "b", xlab = "Number of clusters (K+1)", ylab = "%Explained Variation")
  K_icluster <- c(K_best=diff(c(0,devRatMinBIC)) %>% abs %>% which.min ,method="iClusterPlus")
  rbind(Kbest[[s]], K_icluster)
})


cluster_k <- do.call(rbind, Kbest_all) %>% as.data.frame()
cluster_k <- cluster_k %>% mutate(K_best= K_best %>% as.character %>% as.numeric) 
cluster_k <- cluster_k %>% mutate(method = as.character(method))
cluster_k$method <- gsub("NMF", "intNMF", cluster_k$method)
cluster_k <- cluster_k %>% mutate(method= factor(method, levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster",
                                                                    "intNMF", "MCIA", "iClusterPlus", "PINSPlus")))



g_kbest <- cluster_k %>% 
  ggplot(aes(x = method, y = K_best, color = method)) + 
  geom_boxplot() + 
  theme_bw() +
  geom_hline(yintercept = 3, col = "black", lty = "dashed") + 
  scale_color_manual(values  = MetBrewer::met.brewer("VanGogh2", 8)) + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) + 
  xlab("") + ylab("Number of Best Clusters")

