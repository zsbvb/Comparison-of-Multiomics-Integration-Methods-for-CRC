
rm(list = ls())

library(tidyverse)
library(maftools)
library(bigstatsr)
library(ComplexHeatmap)
library(edgeR)

############################################## 1.Data processing ##################################

##### Clinical
Clinical <- read.csv("Output/Subtype_info_K=5.csv") %>% 
  rename(Tumor_Sample_Barcode = "Patient_ID") %>% 
  mutate(Age_group = ifelse(Age > 60, ">60", "<=60"), 
         Stage = gsub("STAGE", "Stage", Stage), 
         SNF_groups = factor(SNF_groups, levels = c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V"))) %>% 
  arrange(SNF_groups)

##### Mutation
COAD_maf <- read.maf("Data/COAD_SNV_20210405.maf", isTCGA = T, clinicalData = Clinical)
READ_maf <- read.maf("Data/READ_SNV_20210405.maf", isTCGA = T, clinicalData = Clinical)

CRC_maf <- merge_mafs(mafs = c(COAD_maf, READ_maf))
CRC_maf <- subsetMaf(CRC_maf, tsb = CRC_maf@clinical.data$Tumor_Sample_Barcode) ### 285 patients

Mut_mat <- CRC_maf@data
Mut_mat <- Mut_mat %>% 
  select(c("Hugo_Symbol", "Variant_Type", "Tumor_Sample_Barcode")) %>% 
  filter(Variant_Type == "SNP") %>% unique

Mut_mat_wide <- Mut_mat %>% pivot_wider(id_cols = "Hugo_Symbol",
                                        names_from = "Tumor_Sample_Barcode", 
                                        values_from = "Variant_Type")

SNP_num <- rowSums(!is.na(Mut_mat_wide))
names(SNP_num) <- Mut_mat_wide$Hugo_Symbol

##### TOP15
top15 <- sort(SNP_num, decreasing = T)[1:15] %>% names()

Mut_mat_wide_top15 <- Mut_mat_wide %>% filter(Hugo_Symbol %in% top15)
Mut_mat_wide_top15 <- Mut_mat_wide_top15 %>% 
  slice(match(top15, Mut_mat_wide_top15$Hugo_Symbol)) %>% 
  column_to_rownames(var = "Hugo_Symbol")


Clinical_mut <- Clinical %>% filter(Tumor_Sample_Barcode %in% names(Mut_mat_wide_top15))

Mut_mat_wide_top15 <- Mut_mat_wide_top15 %>% select(Clinical_mut$Tumor_Sample_Barcode)
Mut_mat_wide_top15 <- as.matrix(Mut_mat_wide_top15)
#write.csv(Mut_mat_wide_top15, "Output/Mutation.csv")


############################################## 2.Mutated frequency ##################################

snf_groups <- unique(Clinical$SNF_groups)

SNP_SNF <- do.call(rbind, lapply(snf_groups, function(ii) {
  ID <- Clinical_mut$Tumor_Sample_Barcode[which(Clinical_mut$SNF_groups == ii)]
  Mut_mat_sub <- Mut_mat_wide_top15[, ID]
  rowSums(!is.na(Mut_mat_sub))/ncol(Mut_mat_sub)
}))


############################################## 3.OncoPrint ##################################

Mut_mat_wide_top15[is.na(Mut_mat_wide_top15)] <- ""

col <- c("SNP" = "darkgreen")

alter_fun = list(background = alter_graphic("rect", fill = "white"),
                 SNP = alter_graphic("rect", fill = col["SNP"]))

heatmap_legend_param <- list(title = "Alternations",
                             at = c("SNP"),
                             labels = c("Mutation"))

Top = HeatmapAnnotation(`SNF Clusters` = Clinical_mut$SNF_groups,
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 10), border = T,
                                                       title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                       ncol = 1),
                        border = T,
                        col = list(`SNF Clusters` = c("SNF-I" = "#925E9FFF", 
                                                      "SNF-II" = "#0099B4FF", 
                                                      "SNF-III" = "#42B540FF", 
                                                      "SNF-IV" = "#E18727FF", 
                                                      "SNF-V" = "#DC0000CC")),
                        show_annotation_name = TRUE,
                        annotation_name_side = "right",
                        annotation_name_gp = gpar(fontsize = 10))

p1 <- oncoPrint(Mut_mat_wide_top15,
                alter_fun = alter_fun, 
                col = col,
                top_annotation = Top,
                right_annotation = NULL,
                column_order = 1:ncol(Mut_mat_wide_top15),
                column_split = c(rep(1, 70), rep(2, 67), rep(3, 61), rep(4, 39), rep(5, 48)),
                border_gp = gpar(col = "black", lwd = 1.5),
                heatmap_legend_param = heatmap_legend_param)
p1          
#graph2pdf(p1, file = "Output/Oncoprint.pdf", width = 8, height = 3)

########################################## 3.EdgeR #############################################

CRC_MAD <- readRDS("Output/CRC_MAD.rds")
CRC_MAD <- lapply(CRC_MAD, function(ii) {
  ii <- data.frame(ii, check.names = F)
})


load("Data/CORE_Counts.Rdata")

### top1000
Counts <- mergeRawCounts %>% 
  slice(which(rownames(mergeRawCounts) %in% names(CRC_MAD$FPKM)))

Counts <- Counts[, which(substr(names(Counts), 14, 15) %in% "01")]

source("Script/TCGA_operation.R")
id_filter <- tcgaReplicateFilter(tsb = names(Counts), full_barcode = T, analyte_target = "RNA")
Counts <- Counts %>% dplyr::select(all_of(id_filter))

names(Counts) <- substr(names(Counts), 1, 12)
Counts <- Counts %>% select(rownames(CRC_MAD$FPKM))

######## DA
SNF_group <- unique(Clinical$SNF_groups) %>% sort
comb <- combn(SNF_group, 2)

DA_res_counts <- do.call(rbind, lapply(1:ncol(comb), function(ii){
  
  comb_two <- comb[, ii]
  a <- Clinical$Tumor_Sample_Barcode[Clinical$SNF_groups == comb_two[1]]
  b <- Clinical$Tumor_Sample_Barcode[Clinical$SNF_groups == comb_two[2]]
  
  counts_snf <- Counts %>% select(c(a, b))
  group <- c(rep(1, length(a)), rep(2, length(b)))
  
  y <- DGEList(count = counts_snf, group = group)
  
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = F]
  y <- y %>% calcNormFactors() %>% estimateDisp()
  
  differ_res <- y %>% 
    exactTest() %>% 
    topTags(n = 1000) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Var") %>% 
    mutate(Group = paste0(comb_two[1], " VS ", comb_two[2]))
  
  return(differ_res)
}))

DA_res_counts <- DA_res_counts %>% rename(Pvalue = "PValue", Variable = "Var")
# write.csv(DA_res_counts, "Output/DA_res_counts.csv", row.names = F)

DA_res_counts_sig <- DA_res_counts %>% filter(FDR < 0.05 & abs(logFC) > 1.5)
diff_gene <- unique(DA_res_counts_sig$Variable)

library(UpSetR)

upset_df <- DA_res_counts_sig %>% select(c("Variable", "Group")) %>% mutate(Number = 1)
upset_df_wide <- upset_df %>% 
  pivot_wider(id_cols = "Variable", names_from = "Group", values_from = "Number") %>% 
  column_to_rownames(var = "Variable")

upset_df_wide[is.na(upset_df_wide)] <- 0

plot <- upset(upset_df_wide,
              sets = rev(colnames(upset_df_wide)),
              keep.order = T,
              nsets = 20, 
              nintersects = NA,
              set_size.show = T, 
              set_size.scale_max = 200,
              text.scale = 1.5,
              main.bar.color = "cadetblue4", 
              matrix.color = "grey60", 
              point.size = 2, 
              line.size = 0.8, 
              shade.color = 'grey', 
              shade.alpha = 0.2, 
              matrix.dot.alpha = 0.6, 
              mb.ratio = c(0.7, 0.3)) 

#graph2pdf(x = plot, file = "Output/Upset.pdf", width = 18, height = 9)

########################################### 4.Heatmap ##########################################

diff_gene_expr <- CRC_MAD$FPKM %>% 
  select(which(names(CRC_MAD$FPKM) %in% diff_gene)) %>% 
  slice(match(colnames(Mut_mat_wide_top15), rownames(CRC_MAD$FPKM))) %>% t

m <- diff_gene_expr
m[m > 2] <- 2
m[m < -2] <- -2

Cluster <- c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")
names(Cluster) <- levels(Clinical$SNF_groups)

Top <- HeatmapAnnotation(Cluster = Clinical_mut$SNF_groups,
                         annotation_legend_param = list(labels_gp = gpar(fontsize = 10), border = T,
                                                        title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                        ncol=1),
                         border = T,
                         col = list(Cluster = Cluster),
                         show_annotation_name = T,
                         annotation_name_side = "left",
                         annotation_name_gp = gpar(fontsize = 10))

p2 <- Heatmap(m, name = 'Z-score',
             top_annotation = Top,
             cluster_rows = T,
             col = colorRamp2(c(-2, 0, 2), c('paleturquoise', 'white', 'firebrick')),
             color_space = "RGB",
             cluster_columns = F, border = T,
             row_order = NULL,
             column_order = NULL,
             show_row_names = F,
             show_column_names = F,
             row_names_gp = gpar(fontsize = 9),
             column_split = c(rep(1, 70), rep(2, 67), rep(3, 61), rep(4, 39), rep(5, 48)),
             border_gp = gpar(col = "black", lwd = 1.5),
             column_title_gp = gpar(fontsize = 20),
             show_heatmap_legend = T,
             heatmap_legend_param = list(labels_gp = gpar(fontsize = 10), border = T,
                                         title_gp = gpar(fontsize = 10, fontface = "bold")))
p2

graph2pdf(p2, file = "Output/Heatmap_diffgene.pdf", width = 8, height = 4)

