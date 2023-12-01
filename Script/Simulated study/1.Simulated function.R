
CORE_intersim <- function(n.sample, cluster.sample.prop, p.DMP, delta.methyl, delta.expr, delta.protein, 
                          noise.methyl, noise.expr, noise.protein, file){
  
  require(MASS)
  require(InterSIM)
  
  ### sample size and sample size among subgroups
  n.sample <- n.sample
  cluster.sample.prop <- cluster.sample.prop
  
  n.cluster <- length(cluster.sample.prop)
  n.sample.in.cluster <- c(round(cluster.sample.prop[-n.cluster] * n.sample), 
                           n.sample - sum(round(cluster.sample.prop[-n.cluster] * n.sample)))
  
  cluster.id <- do.call(c, lapply(1:n.cluster, function(x) rep(x, n.sample.in.cluster[x])))
  
  ##################################### DNA methylation #####################################
  
  CORE_cov.M <- CORE_intersim_data$CORE_cov.M
  CORE_mean.M <- CORE_intersim_data$CORE_mean.M
  
  n.CpG <- ncol(CORE_cov.M)
  sigma.methyl <- NULL
  
  if (!is.null(sigma.methyl)) {   
    if (sigma.methyl == "indep")
      cov.str <- diag(diag(CORE_cov.M))
    else cov.str <- sigma.methyl
  }else cov.str <- CORE_cov.M
  
  p.DMP <- p.DMP
  delta.methyl <- delta.methyl

  DMP <- sapply(1:n.cluster, function(x) rbinom(n.CpG, 1, prob = p.DMP))
  rownames(DMP) <- names(CORE_mean.M)
  
  d <- lapply(1:n.cluster, function(i) {
    effect <- CORE_mean.M + DMP[, i] * delta.methyl
    mvrnorm(n = n.sample.in.cluster[i], mu = effect, Sigma = cov.str)
  })
  
  biomark_DMP <- lapply(1:n.cluster, function(i) names(which(DMP[, i] == 1)))
  non_DMP <- names(rowSums(DMP))[which(rowSums(DMP) %in% 0)]
  
  if(noise.methyl == FALSE) {
    d <- d 
  } else {
    for (i in 1:n.cluster) {
      index <- match(non_DMP, colnames(d[[i]]))
      for (j in index) {
        d[[i]][, j] <- rnorm(nrow(d[[i]]), mean(d[[i]][, j]), noise.methyl*sd(d[[i]][, j]))
      }
    }
  }
  
  sim.methyl <- do.call(rbind, d)
  sim.methyl <- rev.logit(sim.methyl)
  
  #################################### gene expression ###############################
  
  CORE_cov.expr <- CORE_intersim_data$CORE_cov.expr
  CORE_mean.expr <- CORE_intersim_data$CORE_mean.expr
  CORE_rho.methyl.expr <- CORE_intersim_data$CORE_rho.methyl.expr
  CpG_to_gene <- CORE_intersim_data$CpG_to_gene
  
  n.gene <- ncol(CORE_cov.expr)
  sigma.expr <- NULL
  
  if (!is.null(sigma.expr)){
    if (sigma.expr == "indep") 
      cov.str <- diag(diag(CORE_cov.expr))
    else cov.str <- sigma.expr
  }else cov.str <- CORE_cov.expr
  
  cor.methyl.expr <- NULL
  
  if (!is.null(cor.methyl.expr)){
    CORE_rho.m.e <- cor.methyl.expr
  }else CORE_rho.m.e <- CORE_rho.methyl.expr
  
  p.DEG <- NULL
  
  if (!is.null(p.DEG)) {
    DEG <- sapply(1:n.cluster, function(x) rbinom(n.gene, 1, prob = p.DEG))
    rownames(DEG) <- names(CORE_mean.expr)
  } else {DEG <- sapply(1:n.cluster, function(x){
    cg.name <- rownames(subset(DMP, DMP[, x] == 1))
    gene.name <- as.character(CpG_to_gene[cg.name, ]$tmp.gene)
    as.numeric(names(CORE_mean.expr) %in% gene.name)})
  rownames(DEG) <- names(CORE_mean.expr)
  }
  
  delta.expr <- delta.expr
  sigma.protein <- NULL
  cor.expr.protein <- NULL
  p.DEP <- NULL
  delta.protein <- delta.protein
  
  CORE_methyl.gene.level.mean <- CORE_intersim_data$CORE_methyl.gene.level.mean
  
  d <- lapply(1:n.cluster, function(i) {
    effect <- (CORE_rho.m.e * CORE_methyl.gene.level.mean + sqrt(1 - CORE_rho.m.e^2) * CORE_mean.expr) + DEG[, i] * delta.expr
    mvrnorm(n = n.sample.in.cluster[i], mu = effect, Sigma = cov.str)})

  biomark_DEG <- lapply(1:n.cluster, function(i) names(which(DEG[, i] == 1)))
  non_DEG <- names(rowSums(DEG))[which(rowSums(DEG) %in% 0)]
  
  if(noise.expr == FALSE) {
    d <- d
  } else {
    for (i in 1:n.cluster) {
      index <- match(non_DEG, colnames(d[[i]]))
      for (j in index) {
        d[[i]][, j] <- rnorm(nrow(d[[i]]), mean(d[[i]][, j]), noise.expr*sd(d[[i]][, j]))
      }
    }
  }
  
  sim.expr <- do.call(rbind, d)
  
  ################################# protein expression ##################################
  
  CORE_mean.protein <- CORE_intersim_data$CORE_mean.protein
  CORE_cov.protein <- CORE_intersim_data$CORE_cov.protein
  CORE_rho.expr.protein <- CORE_intersim_data$CORE_rho.expr.protein
  gene_to_protein <- CORE_intersim_data$gene_to_protein
  CORE_mean.expr.with.mapped.protein <- CORE_intersim_data$CORE_mean.expr.with.mapped.protein
  
  n.protein <- ncol(CORE_cov.protein)
  
  if (!is.null(sigma.protein)) {
    if (sigma.protein == "indep")
      cov.str <- diag(diag(CORE_cov.protein))
    else cov.str <- sigma.protein
  }else cov.str <- CORE_cov.protein
  
  
  if (!is.null(cor.expr.protein)){
    CORE_rho.e.p <- cor.expr.protein
  }else CORE_rho.e.p <- CORE_rho.expr.protein
  
  
  if (!is.null(p.DEP)) {
    DEP <- sapply(1:n.cluster, function(x) rbinom(n.protein, 1, prob = p.DEP))
    rownames(DEP) <- names(CORE_mean.protein)
  } else {
    DEP <- sapply(1:n.cluster, function(x) {
      gene.name <- rownames(subset(DEG, DEG[, x] == 1))
      protein.name <- rownames(gene_to_protein[gene_to_protein$gene %in% gene.name, ])
      as.numeric(names(CORE_mean.protein) %in% protein.name)
    })
    rownames(DEP) <- names(CORE_mean.protein)
  }
  
  d <- lapply(1:n.cluster, function(i) {
    effect <- (CORE_rho.e.p * CORE_mean.expr.with.mapped.protein + sqrt(1 - CORE_rho.e.p^2) * CORE_mean.protein) + DEP[, i] * delta.protein
    mvrnorm(n = n.sample.in.cluster[i], mu = effect, Sigma = cov.str)})

  biomark_DEP <- lapply(1:n.cluster, function(i) names(which(DEP[, i] == 1)))
  non_DEP <- names(rowSums(DEP))[which(rowSums(DEP) %in% 0)]
  
  if(noise.protein == FALSE) {
    d <- d 
  } else {
    for (i in 1:n.cluster) {
      index <- match(non_DEP, colnames(d[[i]]))
      for (j in index) {
        d[[i]][, j] <- rnorm(nrow(d[[i]]), mean(d[[i]][, j]), noise.protein*sd(d[[i]][, j]))
      }
    }
  }
  
  sim.protein <- do.call(rbind, d)
  
  indices <- sample(1:n.sample)
  cluster.id <- cluster.id[indices]
  sim.methyl <- sim.methyl[indices, ]
  sim.expr <- sim.expr[indices, ]
  sim.protein <- sim.protein[indices, ]
  rownames(sim.methyl) <- paste("subject", 1:nrow(sim.methyl), sep = "")
  rownames(sim.expr) <- paste("subject", 1:nrow(sim.expr), sep = "")
  rownames(sim.protein) <- paste("subject", 1:nrow(sim.protein), sep = "")
  d.cluster <- data.frame(rownames(sim.methyl), cluster.id)
  colnames(d.cluster)[1] <- "subjects"
  
  sim <- list(data = list(dat.methyl = sim.methyl, dat.expr = sim.expr, dat.protein = sim.protein), 
              biomark = list(biomark_DMP = biomark_DMP, biomark_DEG = biomark_DEG, biomark_DEP = biomark_DEP), 
              clustering.assignment = d.cluster)
  sim %>% saveRDS(file = file)
}
