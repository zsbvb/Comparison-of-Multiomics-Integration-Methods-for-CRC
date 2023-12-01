
######### CIMLR
dofilter_CIMLR <- function(truth, fit){
  selectVars_1 <- fit$selectfeatures$names
  k_grid <- stringr::str_extract(pattern = "_dat*.", selectVars_1) %>% unlist %>% unique %>% sort()
  filter <- lapply(k_grid, function (kk){
    idx <- grep(kk, selectVars_1)
    gsub(kk, "", selectVars_1[idx])
  })
  
  N <- sapply(filter, length)
  
  selectVars_1 <- fit$selectfeatures$names[fit$selectfeatures$pval < 0.05]
  k_grid <- stringr::str_extract(pattern = "_dat*.", selectVars_1) %>% unlist %>% unique %>% sort()
  filter <- lapply(k_grid, function (kk){
    idx <- grep(kk, selectVars_1)
    gsub(kk, "", selectVars_1[idx])
  })
  
  filter_tp <- sapply(filter, length)
  truth_tp <- sapply(truth, length) %>% as.numeric()
  
  TPR <- sapply(1:length(truth_tp), function (ii) {
    filter[[ii]] %>% intersect(truth[[ii]]) %>% length()/truth_tp[ii]
  })
  
  FPR <- sapply(1:length(truth_tp), function (ii) {
    filter[[ii]] %>% setdiff(truth[[ii]]) %>% length()/(N[ii]-truth_tp[ii])
  })
  
  ACC <- sapply(1:length(truth_tp), function (ii) {
    ((filter[[ii]] %>% intersect(truth[[ii]]) %>% length()) + 
       (filter[[ii]] %>% setdiff(truth[[ii]]) %>% length()))/(N[ii])
  })
  
  Rank <- sapply(1:length(truth_tp), function (ii) {
    match <- filter[[ii]] %>% intersect(truth[[ii]]) %>% match(filter[[ii]])
    rank <- match %>% paste(length(filter[[ii]]), sep = "-")
  })
  return(list(TPR = TPR, FPR = FPR, ACC = ACC, Rank = Rank))
}

########## Mocluster 

dofilter_Moa <- function(truth, fit){
  K <- fit@data %>% length
  a <- fit@loading

  selectVars_1 <- which(a %>% rowSums != 0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  test <- rowSums(abs(a))
  idx <- which(test != 0)
  filter <- sort(test[idx], decreasing = TRUE) %>% names
  
  filter <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), filter)
    gsub(sprintf("_dat%s", kk), "", filter[idx])})
  
  N <- sapply(fit@data, nrow)
  filter_tp <- sapply(filter, length)
  truth_tp <- sapply(truth, length) %>% as.numeric()
  
  TPR <- sapply(1:length(truth_tp), function (ii) {
    filter[[ii]] %>% intersect(truth[[ii]]) %>% length()/truth_tp[ii]
  })
  
  FPR <- sapply(1:length(truth_tp), function (ii) {
    filter[[ii]] %>% setdiff(truth[[ii]]) %>% length()/(N[ii]-truth_tp[ii])
  })
  
  ACC <- sapply(1:length(truth_tp), function (ii) {
    ((filter[[ii]] %>% intersect(truth[[ii]]) %>% length()) + 
       (filter[[ii]] %>% setdiff(truth[[ii]]) %>% length()))/(N[ii])
  })
  
  Rank <- sapply(1:length(truth_tp), function (ii) {
    match <- filter[[ii]] %>% intersect(truth[[ii]]) %>% match(filter[[ii]])
    rank <- match %>% paste(length(filter[[ii]]), sep = "-")
  })
  return(list(TPR = TPR, FPR = FPR, ACC = ACC, Rank = Rank))
}

######## MCIA

dofilter_MCIA <- function(truth, fit){
  coas <- fit$coa
  a <- lapply(coas, function (cc) cc$li)
  
  test <- lapply(a, function(aa) rowSums(abs(aa)))
  filter <- lapply(test, function (tt){
    sort(tt, decreasing = TRUE) %>% names 
  })
  
  N <- sapply(a, nrow)
  filter_tp <- sapply(filter, length)
  truth_tp <- sapply(truth, length) %>% as.numeric()
  
  Rank <- sapply(1:length(truth_tp), function (ii) {
    match <- filter[[ii]] %>% intersect(truth[[ii]]) %>% match(filter[[ii]])
    rank <- match %>% paste(length(filter[[ii]]), sep = "-")
  })
  return(list(Rank = Rank))
}

########## intNMF

dofilter_intNMF <- function(truth, fit){
  a <- fit$H
  
  test <- lapply(a, function(aa) colSums((aa)))
  filter <- lapply(test, function (tt){
    sort(tt, decreasing = TRUE) %>% names 
  })
  
  N <- sapply(a, ncol)
  filter_tp <- sapply(filter, length)
  truth_tp <- sapply(truth, length) %>% as.numeric()
  
  Rank <- sapply(1:length(truth_tp), function (ii) {
    match <- filter[[ii]] %>% intersect(truth[[ii]]) %>% match(filter[[ii]])
    rank <- match %>% paste(length(filter[[ii]]), sep = "-")
  })
  return(list(Rank = Rank))
}


Var_importance <- function(truth, fit, method){
  
  dofilter <- switch(method,
                     "NMF" =  dofilter_intNMF,
                     "MCIA" =  dofilter_MCIA,
                     "Mocluster" = dofilter_Moa, 
                     "CIMLR" = dofilter_CIMLR)
  
  res <- dofilter(truth,fit)
  
  Eval <- do.call(c, res$Rank) %>% 
    data.frame() %>% 
    mutate(Datatype = c(rep("DNA Methylation", length(res$Rank[[1]])), 
                        rep("Gene Expression", length(res$Rank[[2]])), 
                        rep("Protein Expression", length(res$Rank[[3]]))))
  colnames(Eval) <- c("Rank", "Datatype")

  return(list(Rank = Eval))
}


listBenchmark <- gtools::mixedsort(list.files(pathDat))
S <- 50

importance_eval_dat <- lapply(1:length(listBenchmark), function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  print(pathMeth_sub)
  print(pathDat_sim)
  
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  
  trueDat1 <- sapply(list.sim, function (ss) ss$biomark$biomark_DMP %>% unlist %>% unique)
  trueDat2 <- sapply(list.sim, function (ss) ss$biomark$biomark_DEG %>% unlist %>% unique)
  trueDat3 <- sapply(list.sim, function (ss) ss$biomark$biomark_DEP %>% unlist %>% unique)
  
  truth <- lapply(1:S, function (ss) list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]]))
  
  mm <- "CIMLR"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit)
  l <- list(truth=truth, fits)
  importance_eval_CIMLR <- purrr::pmap(l, Var_importance, method = mm)
  
  
  mm <- "Mocluster"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  importance_eval_moclust <- purrr::pmap(l, Var_importance, method = mm)
  
  
  mm <- "NMF"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  importance_eval_nmf <- purrr::pmap(l, Var_importance, method = mm)
  
  
  mm <- "MCIA"
  print(mm)
  pp <- list.files(pathMeth_sub, pattern=mm, full.names = TRUE)
  ff <- readRDS(pp)
  fits <- lapply(1:S, function (ss) ff[[ss]]$fit )
  l <- list(truth=truth, fits)
  importance_eval_mcia <- purrr::pmap(l, Var_importance, method = mm)

  importance_eval_CIMLR <- lapply(1:S, function(ss) {
    importance_eval_CIMLR[[ss]]$Rank <- importance_eval_CIMLR[[ss]]$Rank %>% mutate(Simu = paste0("simu", ss))
  })

  importance_CIMLR <- data.frame(do.call(rbind, importance_eval_CIMLR))
  importance_CIMLR <- importance_CIMLR %>% 
    mutate(Method = 'CIMLR', TOP = importance_CIMLR$Rank %>% as.character() %>% strsplit("-") %>% 
             sapply("[", 1) %>% as.numeric(), noise = b)
  
  
  importance_eval_moclust <- lapply(1:S, function(ss) {
    importance_eval_moclust[[ss]]$Rank <- importance_eval_moclust[[ss]]$Rank %>% mutate(Simu = paste0("simu", ss))
  })
  importance_Mocluster <- data.frame(do.call(rbind, importance_eval_moclust))
  importance_Mocluster <- importance_Mocluster %>% 
    mutate(Method = 'Mocluster', TOP = importance_Mocluster$Rank %>% as.character() %>% strsplit("-") %>% 
             sapply("[", 1) %>% as.numeric(), noise = b)
  
  
  importance_eval_nmf <- lapply(1:S, function(ss) {
    importance_eval_nmf[[ss]]$Rank <- importance_eval_nmf[[ss]]$Rank %>% mutate(Simu = paste0("simu", ss))
  })
  importance_intNMF <- data.frame(do.call(rbind, importance_eval_nmf))
  importance_intNMF <- importance_intNMF %>% 
    mutate(Method = 'intNMF', TOP = importance_intNMF$Rank %>% as.character() %>% strsplit("-") %>% 
             sapply("[", 1) %>% as.numeric(), noise = b)
  
  
  importance_eval_mcia <- lapply(1:S, function(ss) {
    importance_eval_mcia[[ss]]$Rank <- importance_eval_mcia[[ss]]$Rank %>% mutate(Simu = paste0("simu", ss))
  })
  importance_MCIA <- data.frame(do.call(rbind, importance_eval_mcia))
  importance_MCIA <- importance_MCIA %>% 
    mutate(Method = 'MCIA', TOP = importance_MCIA$Rank %>% as.character() %>% strsplit("-") %>% 
             sapply("[", 1) %>% as.numeric(), noise = b)
  
  return(list(CIMLR = importance_CIMLR, Mocluster = importance_Mocluster, 
              intNMF = importance_intNMF, MCIA = importance_MCIA))
})


################### TOP10
TOP10 <- do.call(cbind, lapply(importance_eval_dat, function(ii) {
  do.call(rbind, lapply(ii, function(iii) {
    a <- iii %>% group_by(Datatype, Simu, Method) %>% summarise(Number = length(TOP[which(TOP <= 10)]))
    b <- a %>% group_by(Datatype, Method) %>% summarise(mean_Number = mean(Number) %>% sprintf("%.1f", .) %>% as.numeric())
    }))
}))

top_10 <- TOP10[, -c(seq(4, 34, by = 3), seq(5, 35, by = 3))]
colnames(top_10) <- c("Datatype", "Method", paste0('Benchmark', 1:12))
top_10 <- top_10 %>% mutate(Datatype = factor(Datatype, levels = c("DNA Methylation", "Gene Expression", "Protein Expression")), 
                            Method = factor(Method, levels = c('CIMLR', 'Mocluster', 'intNMF', 'MCIA')))

library(pheatmap)
top_10_heat <- top_10[, 3:14] %>% data.frame()
rownames(top_10_heat) <- paste(top_10$Method, top_10$Datatype, sep = "_")

pdf(paste0(pathFig, "/TOP10_heatmap.pdf"), width = 10, height = 10)
pheatmap(as.matrix(top_10_heat), color = colorRampPalette(c("white", "#eebe04", "#bd3106"))(50), 
         cluster_cols = F, cluster_rows = F, show_colnames = T, show_rownames = T, fontsize = 11, 
         cellwidth = 25, cellheight = 20, main = "", 
         display_numbers = T, number_format = "%.1f", number_color = "black", border_color = "white")
dev.off()


top_10_long <- gather(top_10, Benchmark, Rank, Benchmark1:Benchmark12)

g_top10_boxplot <- top_10_long %>% ggplot(aes(y = Rank, x = Method, fill = Method,)) +  
  geom_boxplot() + 
  facet_wrap(. ~ Datatype, ncol = 3) + 
  scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2", 4)) + 
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1),
        strip.text.x = element_text(size = 18), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 18), 
        legend.position = "none")


################### TOP20

TOP20 <- do.call(cbind, lapply(importance_eval_dat, function(ii) {
  do.call(rbind, lapply(ii, function(iii) {
    a <- iii %>% group_by(Datatype, Simu, Method) %>% summarise(Number = length(TOP[which(TOP <= 20)]))
    b <- a %>% group_by(Datatype, Method) %>% summarise(mean_Number = mean(Number) %>% sprintf("%.1f", .) %>% as.numeric())
  }))
}))

top_20 <- TOP20[, -c(seq(4, 34, by = 3), seq(5, 35, by = 3))]
colnames(top_20) <- c("Datatype", "Method", paste0('Benchmark', 1:12))
top_20 <- top_20 %>% mutate(Datatype = factor(Datatype, levels = c("DNA Methylation", "Gene Expression", "Protein Expression")), 
                            Method = factor(Method, levels = c('CIMLR', 'Mocluster', 'intNMF', 'MCIA')))

library(pheatmap)
top_20_heat <- top_20[, 3:14] %>% data.frame()
rownames(top_20_heat) <- paste(top_20$Method, top_20$Datatype, sep = "_")

pheatmap(as.matrix(top_20_heat), color = colorRampPalette(c("white", "#eebe04", "#bd3106"))(50), 
         cluster_cols = F, cluster_rows = F, show_colnames = T, show_rownames = T, fontsize = 11, 
         cellwidth = 25, cellheight = 20, main = "", 
         display_numbers = T, number_format = "%.1f", number_color = "black", border_color = "white")

top_20_long <- gather(top_20, Benchmark, Rank, Benchmark1:Benchmark12)

g_top20_boxplot <- top_20_long %>% ggplot(aes(y = Rank, x = Method, fill = Method,)) +  
  geom_boxplot() + 
  facet_wrap(. ~ Datatype, ncol = 3) + 
  scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2", 4)) + 
  xlab("") + ylab("") + 
  theme_bw() + 
  theme(panel.border = element_rect(linewidth = 1),
        strip.text.x = element_text(size = 18), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 18), 
        legend.position = "none")
