
rm(list = ls())

library(InterSIM)
library(tidyverse)
library(impute)

################################## cpg to gene ##############################
CpG_to_gene <- CpG.gene.map.for.DEG

CORE_beta <- read_tsv("Data/HumanMethylation27.tsv")
CORE_beta <- tibble::column_to_rownames(CORE_beta, var = "sample")
index1 <- intersect(rownames(CORE_beta), CpG_to_gene$tmp.cg)

CORE_FPKM <- read_tsv("Data/AgilentG4502A_07_3.tsv")
CORE_FPKM <- tibble::column_to_rownames(CORE_FPKM, var = "sample")
index2 <- intersect(rownames(CORE_FPKM), CpG_to_gene$tmp.gene)

CpG_to_gene$tmp.cg <- as.character(CpG_to_gene$tmp.cg)
CpG_to_gene$tmp.gene <- as.character(CpG_to_gene$tmp.gene)

CpG_to_gene <- CpG_to_gene[which(CpG_to_gene$tmp.gene %in% index2), ]

####### missing value
cpg_na <- apply(CORE_beta[CpG_to_gene$tmp.cg,], 1, function(x) sum(is.na(x))/ncol(CORE_beta))
cpg_sel <- cpg_na[cpg_na < 0.2]

CpG_to_gene <- CpG_to_gene[which(CpG_to_gene$tmp.cg %in% names(cpg_sel)), ]

################################## gene to protein ##############################
gene_to_protein <- protein.gene.map.for.DEP
gene_to_protein$protein <- toupper(rownames(gene_to_protein))

gene_to_protein <- gene_to_protein[which(gene_to_protein$gene %in% CpG_to_gene$tmp.gene), ]

gene_to_protein$gene <- as.character(gene_to_protein$gene)
gene_to_protein$protein <- as.character(gene_to_protein$protein)

CORE_protein <- read.csv("CORE_proteomics.csv", row.names = 1, check.names = F)
CORE_protein <- CORE_protein[, -c(2, 3, 4)]
CORE_protein <- data.frame(t(CORE_protein), check.names = F)
colnames(CORE_protein) <- CORE_protein[1, ]
CORE_protein <- CORE_protein[-1, ]
rownames(CORE_protein) <- toupper(rownames(CORE_protein))

index3 <- intersect(rownames(CORE_protein), gene_to_protein$protein)
gene_to_protein <- gene_to_protein[which(gene_to_protein$protein %in% index3), ]

CpG_to_gene <- CpG_to_gene[which(CpG_to_gene$tmp.gene %in% gene_to_protein$gene), ]

########## filter cpg, gene and protein
CORE_beta_intersim <- CORE_beta[unique(CpG_to_gene$tmp.cg), ]
CORE_FPKM_intersim <- CORE_FPKM[unique(CpG_to_gene$tmp.gene), ]
CORE_protein_intersim <- CORE_protein[unique(gene_to_protein$protein), ]


### tumor sample
cpg_TCGAid <- colnames(CORE_beta)
cpg_tumor_id <- cpg_TCGAid[which(as.numeric(substr(cpg_TCGAid, 14, 15)) == 1)]

gene_TCGAid <- colnames(CORE_FPKM)
gene_tumor_id <- gene_TCGAid[which(as.numeric(substr(gene_TCGAid, 14, 15)) == 1)]

protein_TCGAid <- colnames(CORE_protein)
protein_tumor_id <- protein_TCGAid[which(as.numeric(substr(protein_TCGAid, 14, 15)) == 1)]

source("Script/Simulated study/TCGA_operation.R")

cpg_filter <- tcgaReplicateFilter(tsb = cpg_tumor_id, analyte_target = "DNA")
gene_filter <- tcgaReplicateFilter(tsb = gene_tumor_id, analyte_target = "RNA")
protein_filter <- tcgaReplicateFilter(tsb = protein_tumor_id, analyte_target = "RNA")

CORE_beta_intersim_filter <- CORE_beta_intersim[, cpg_filter]
colnames(CORE_beta_intersim_filter) <- substr(colnames(CORE_beta_intersim_filter), 1, 12)

CORE_FPKM_intersim_filter <- CORE_FPKM_intersim[, gene_filter]
colnames(CORE_FPKM_intersim_filter) <- substr(colnames(CORE_FPKM_intersim_filter), 1, 12)

CORE_protein_intersim_filter <- CORE_protein_intersim[, protein_filter]
colnames(CORE_protein_intersim_filter) <- substr(colnames(CORE_protein_intersim_filter), 1, 12)


intersect_id <- intersect(colnames(CORE_protein_intersim_filter), colnames(CORE_FPKM_intersim_filter)) %>% 
  intersect(colnames(CORE_protein_intersim_filter))

cpgdata_intersim <- CORE_beta_intersim_filter[, intersect_id]
genedata_intersim <- CORE_FPKM_intersim_filter[, intersect_id]
proteindata_intersim <- CORE_protein_intersim_filter[, intersect_id]

### KNN imputation
cpgdata_intersim_knn <- impute.knn(as.matrix(cpgdata_intersim))
cpgdata_intersim_knn <- cpgdata_intersim_knn$data
cpgdata_intersim_knn <- data.frame(cpgdata_intersim_knn, check.names = F)

cpgdata_intersim_t <- data.frame(t(cpgdata_intersim_knn))
proteindata_intersim_t <- data.frame(t(proteindata_intersim))
genedata_intersim_t <- data.frame(t(genedata_intersim))

######################## calculate correlation using in simulated study #############################

### mean value
for (i in 1:93) {
  proteindata_intersim_t[, i] <- as.numeric(proteindata_intersim_t[, i])
}

CORE_mean.protein <- apply(proteindata_intersim_t, 2, mean)
CORE_mean.expr <- apply(genedata_intersim_t, 2, mean)
CORE_mean.M <- apply(logit(cpgdata_intersim_t), 2, mean)

### variance-covariance matrix
CORE_cov.protein <- cov(proteindata_intersim_t)
CORE_cov.expr <- cov(genedata_intersim_t)
CORE_cov.M <- cov(logit(cpgdata_intersim_t))

### the corresponding methy.mean in gene
cpgdata_intersim_t_long <- tibble::rownames_to_column(cpgdata_intersim_t, var = "ID")
cpgdata_intersim_t_long <- reshape2::melt(cpgdata_intersim_t_long, id = "ID")

cpgdata_intersim_t_long_gene <- CpG_to_gene$tmp.gene[match(cpgdata_intersim_t_long$variable, CpG_to_gene$tmp.cg)]
cpgdata_intersim_t_long <- mutate(cpgdata_intersim_t_long, gene = cpgdata_intersim_t_long_gene)

methyl.gene.level <- cpgdata_intersim_t_long %>% group_by(ID, gene) %>% summarise(methyl.gene.level = median(value))

methyl.gene.level_wide <- reshape2::dcast(methyl.gene.level, ID ~ gene)
methyl.gene.level_wide <- tibble::column_to_rownames(methyl.gene.level_wide, var = "ID")
methyl.gene.level_wide <- logit(methyl.gene.level_wide)

CORE_methyl.gene.level.mean <- apply(methyl.gene.level_wide, 2, mean)

### the corresponding expr.mean in protein
index_gene <- match(gene_to_protein$gene, colnames(genedata_intersim_t))
CORE_mean.expr.with.mapped.protein <- apply(genedata_intersim_t, 2, mean)[index_gene]
names(CORE_mean.expr.with.mapped.protein) <- gene_to_protein$protein

### correlation between cpg and gene

cor.res <- data.frame()
for (i in 1:81) {
  cpg <- as.numeric(methyl.gene.level_wide[, i])
  gene <- as.numeric(genedata_intersim_t[, i])
  cor.fit <- cor.test(cpg, gene, method = "pearson", na.action = na.omit)
  
  cor.res <- rbind(cor.res, c(as.character(colnames(genedata_intersim_t)[i]), 
                              sprintf("%.9f", cor.fit$estimate), 
                              sprintf("%.3f", cor.fit$p.value)))
}

colnames(cor.res) <- c("gene_id", "cor", "p.val")

CORE_rho.methyl.expr <- as.numeric(cor.res$cor)
names(CORE_rho.methyl.expr) <- cor.res$gene_id

### correlation between gene and protein
cor.res2 <- data.frame()
for (i in 1:93) {
  gene <- as.numeric(genedata_intersim_t[, as.character(gene_to_protein$gene[i])])
  protein <- as.numeric(proteindata_intersim_t[, as.character(gene_to_protein$protein[i])])
  cor.fit <- cor.test(gene, protein, method = "pearson")
  
  cor.res2 <- rbind(cor.res2, c(as.character(gene_to_protein$gene[i]), 
                                as.character(gene_to_protein$protein[i]), 
                                sprintf("%.9f", cor.fit$estimate), 
                                sprintf("%.3f", cor.fit$p.value)))
}

colnames(cor.res2) <- c("gene_id", "protein_id", "cor", "p.val")

CORE_rho.expr.protein <- as.numeric(cor.res2$cor)
names(CORE_rho.expr.protein) <- cor.res2$protein_id

### product list
CORE_intersim_data<- list(CORE_mean.M = CORE_mean.M, 
                          CORE_mean.expr = CORE_mean.expr, 
                          CORE_mean.protein = CORE_mean.protein, 
                          CORE_cov.M = CORE_cov.M,
                          CORE_cov.expr = CORE_cov.expr, 
                          CORE_cov.protein = CORE_cov.protein, 
                          CpG_to_gene = CpG_to_gene,
                          gene_to_protein = gene_to_protein, 
                          CORE_methyl.gene.level.mean = CORE_methyl.gene.level.mean,
                          CORE_mean.expr.with.mapped.protein = CORE_mean.expr.with.mapped.protein, 
                          CORE_rho.methyl.expr = CORE_rho.methyl.expr,
                          CORE_rho.expr.protein = CORE_rho.expr.protein, 
                          cpgdata = cpgdata_intersim_t, 
                          genedata = genedata_intersim_t, 
                          proteindata = proteindata_intersim_t)


