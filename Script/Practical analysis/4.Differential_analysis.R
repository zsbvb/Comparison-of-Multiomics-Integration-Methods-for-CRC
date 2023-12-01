
rm(list = ls())

library(tidyverse)
library(GSVA)
library(ConsensusTME)
library(ggpubr)

Clinical <- read.csv("Output/Subtype_info_K=5.csv") %>% 
  rename(Tumor_Sample_Barcode = "Patient_ID") %>% 
  mutate(Age_group = ifelse(Age > 60, ">60", "<=60"), 
         Stage = gsub("STAGE", "Stage", Stage), 
         SNF_groups = factor(SNF_groups, levels = c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V"))) %>% 
  arrange(SNF_groups)

Clinical_phenotype <- read.table("Data/Clinical_sample.txt", 
                                 header = T, 
                                 sep = '\t',                 
                                 quote = "",                              
                                 comment.char = "#",         
                                 fill = T)

Clinical_phenotype <- Clinical_phenotype %>% 
  select(c(1, 14, 15, 17)) %>% 
  rename(Tumor_Sample_Barcode = "PATIENT_ID")

Clinical_phenotype <- merge(Clinical_phenotype, Clinical %>% select(c(1, 18)), by = "Tumor_Sample_Barcode")

########## TMB
Box_tmb <- ggplot(Clinical_phenotype, aes(SNF_groups, TMB_NONSYNONYMOUS, fill = SNF_groups)) + 
  geom_boxplot(outlier.alpha = 0) + 
  scale_fill_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
  scale_x_discrete(labels = c("SNF-I" = "SNF-I (n = 70)", 
                              "SNF-II" = "SNF-II (n = 71)", 
                              "SNF-III" = "SNF-III (n = 64)", 
                              "SNF-IV" = "SNF-IV (n = 50)", 
                              "SNF-V" = "SNF-V (n = 46)")) + 
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  theme_bw() +
  labs(x = "", y = "Nonsynonymous TMB") +
  theme(panel.border = element_rect(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15),
        legend.position = "none")
Box_tmb
ggsave("Output/Boxplot_TMB.pdf", Box_tmb, width = 3, height = 4)

########## MSI
Box_msi <- ggplot(Clinical_phenotype, aes(SNF_groups, MSI_SCORE_MANTIS, fill = SNF_groups)) + 
  geom_boxplot(outlier.alpha = 0.8) + 
  scale_fill_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
  scale_x_discrete(labels = c("SNF-I" = "SNF-I (n = 70)", 
                              "SNF-II" = "SNF-II (n = 69)", 
                              "SNF-III" = "SNF-III (n = 66)", 
                              "SNF-IV" = "SNF-IV (n = 57)", 
                              "SNF-V" = "SNF-V (n = 48)")) + 
  theme_bw() +
  labs(x = "", y = "MSI Score Mantis") +
  theme(panel.border = element_rect(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15),
        legend.position = "none") + 
  stat_compare_means(aes(group = SNF_groups, label = ..p.signif..), 
                     method = "kruskal.test", size = 5, label.x = 3, label.y = 1.3) 

Box_msi
ggsave("Output/Boxplot_MSI.pdf", Box_msi, width = 3, height = 4)


###### PD-L1 

PDL1 <- mRNA_FPKM %>% slice(which(rownames(mRNA_FPKM) %in% c("CD274", "TP53"))) %>% 
  t %>% data.frame(check.names = F) %>% rownames_to_column(var = "Tumor_Sample_Barcode")

PDL1 <- PDL1 %>% merge(Clinical %>% select(c(1, 18)), by = "Tumor_Sample_Barcode")

Box_PDL1 <- ggplot(PDL1, aes(SNF_groups, CD274, fill = SNF_groups)) + 
  geom_boxplot(outlier.alpha = 0.8) + 
  scale_fill_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
  scale_x_discrete(labels = c("SNF-I" = "SNF-I (n = 70)", 
                              "SNF-II" = "SNF-II (n = 72)", 
                              "SNF-III" = "SNF-III (n = 71)", 
                              "SNF-IV" = "SNF-IV (n = 66)", 
                              "SNF-V" = "SNF-V (n = 52)")) + 
  theme_bw() +
  labs(x = "", y = "PD-L1 Expression") +
  theme(panel.border = element_rect(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15),
        legend.position = "none") + 
  stat_compare_means(aes(group = SNF_groups, label = ..p.signif..), 
                     method = "kruskal.test", size = 5, label.x = 3, label.y = 4.5) 
Box_PDL1
ggsave("Output/CRC/Figure/Boxplot_PDL1.pdf", Box_PDL1, width = 3, height = 4)


###### TP53
Box_TP53 <- ggplot(PDL1, aes(SNF_groups, TP53, fill = SNF_groups)) + 
  geom_boxplot(outlier.alpha = 0.8) + 
  scale_fill_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
  scale_x_discrete(labels = c("SNF-I" = "SNF-I (n = 70)", 
                              "SNF-II" = "SNF-II (n = 72)", 
                              "SNF-III" = "SNF-III (n = 71)", 
                              "SNF-IV" = "SNF-IV (n = 66)", 
                              "SNF-V" = "SNF-V (n = 52)")) + 
  theme_bw() +
  labs(x = "", y = "TP53 Expression") +
  theme(panel.border = element_rect(linewidth = 1),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15),
        legend.position = "none") + 
  stat_compare_means(aes(group = SNF_groups, label = ..p.signif..), 
                     method = "kruskal.test", size = 5, label.x = 3, label.y = 6.3) 
Box_TP53
ggsave("Output/CRC/Figure/Boxplot_TP53.pdf", Box_TP53, width = 3, height = 4)

