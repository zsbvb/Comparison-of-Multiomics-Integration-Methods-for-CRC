
rm(list = ls())

library(tidyverse)
library(survival)
library(survminer)
library(CMSclassifier)
library(org.Hs.eg.db)
library(clusterProfiler)

load("Data/CRC/Mergedata.Rdata")

###################################### 1.CMS subtypes #####################################
fpkm <- Mergedata$FPKM

CMS_entrez <- CMSclassifier::listModelGenes("RF")
CMS.entrez_symbol <- bitr(CMS_entrez,
                          fromType = "ENTREZID",
                          toType = "SYMBOL", 
                          OrgDb = "org.Hs.eg.db")

CMS_genes <- CMS.entrez_symbol$SYMBOL %>% na.omit
CMS_data <- fpkm %>% slice(which(rownames(fpkm) %in% CMS_genes))

setdiff(CMS_genes, rownames(CMS_data)) ### Not matching MRC1

CMS_data[273, ] <- rep(NA, 331)
rownames(CMS_data)[273] <- "MRC1"

CMS_data <- CMS_data %>% slice(match(CMS.entrez_symbol$SYMBOL, rownames(CMS_data)))
rownames(CMS_data) <- CMS.entrez_symbol$ENTREZID

Rfcms <- classifyCMS(CMS_data, method = "RF")[[3]]
#saveRDS(Rfcms, "Output/CRC_CMS_classification.rds")


#################################### 2.Classification accuracy ###################################

Clinical <- Mergedata$Clinical
Clinical <- Clinical %>% 
  mutate(Stage = gsub("[A-D]$", "", Clinical$AJCC_stage)) %>%
  mutate(Stage_number = case_when(Stage == "STAGE I" ~ 1, 
                                  Stage == "STAGE II" ~ 2,
                                  Stage == "STAGE III" ~ 3,
                                  Stage == "STAGE IV" ~ 4))

Rfcms <- Rfcms %>% 
  mutate(CMS_groups = sapply(strsplit(RF.nearestCMS, ","), "[", 1))

table(Clinical$Patient_ID == rownames(Rfcms)) # TRUE 331

Clinical <- Clinical %>% 
  mutate(CMS_groups = Rfcms$CMS_groups, 
         CMS_num = case_when(CMS_groups == "CMS1" ~ 1, 
                             CMS_groups == "CMS2" ~ 2,
                             CMS_groups == "CMS3" ~ 3,
                             CMS_groups == "CMS4" ~ 4), 
         OS_status = ifelse(OS_status == "0:LIVING", 0, 1),
         PFS_status = ifelse(PFS_status == "0:CENSORED", 0, 1))
  
methods <- c("LRAcluster", "SNF", "CIMLR", "Mocluster", "intNMF", "MCIA", "iCluster", "PINSPlus")
path <- "Output/Clustering_res"

########## ARI,NMI,F1
CRC_res <- do.call(rbind, lapply(methods, function(mm){
  print(mm)
  ff <- file.path(path, sprintf("%s_res.rds", mm))
  r <- readRDS(ff)
  
  ARI <- r$clust %>% mclust::adjustedRandIndex(Clinical$CMS_num)
  NMI <- r$clust %>% aricode::NMI(Clinical$CMS_num)
  F1 <- as.factor(r$clust) %>% FlowSOM::FMeasure(Clinical$CMS_num)
  
  data.frame(method = mm, ARI = ARI, NMI = NMI, F1 = F1)
}))
names(CRC_res)[4] <- "F1-score"
#write.csv(CRC_res, "Output/CRC_res.csv", row.names = F)

CRC_res_long <- CRC_res %>% pivot_longer(2:4, names_to = "Index", values_to = "Value")

CRC_res_long <- CRC_res_long %>% 
  mutate(Value = as.numeric(sprintf("%.3f", Value)),
         method = gsub("iCluster", "iClusterPlus", method) %>% 
           factor(levels = c("LRAcluster", "SNF", "CIMLR", "Mocluster", "intNMF", "MCIA", "iClusterPlus", "PINSPlus")), 
         Index = factor(Index, levels = c("ARI", "NMI", "F1-score")))

g_CRC <- CRC_res_long %>% ggplot(aes(method, Value, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(label = Value, hjust = 0.5, vjust = -0.2), size = 3) +
  scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2", 8)) +
  theme_bw() + xlab("") + ylab("") + 
  facet_wrap(. ~ Index, ncol = 3, scales = "free_y", labeller = label_wrap_gen(multi_line = T)) + 
  theme(panel.border = element_rect(linewidth = 0.5),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12), 
        legend.position = "top", 
        legend.text = element_text(size = 12),
        legend.title = element_blank())
g_CRC


######################################### 3.SNF #######################################

CRC_SNF_res <- file.path(path, sprintf("SNF_res_k.rds")) %>% readRDS
CRC_SNF_res_k <- do.call(rbind, lapply(1:length(CRC_SNF_res), function(ii){

  r <- CRC_SNF_res[[ii]]
  
  ARI <- r$clust %>% mclust::adjustedRandIndex(Clinical$CMS_num)
  NMI <- r$clust %>% aricode::NMI(Clinical$CMS_num)
  F1 <- as.factor(r$clust) %>% FlowSOM::FMeasure(Clinical$CMS_num)
  
  data.frame(k_number = names(CRC_SNF_res)[ii], ARI = ARI, NMI = NMI, F1 = F1)
}))
names(CRC_SNF_res_k)[4] <- "F1-score"
#write.csv(CRC_SNF_res_k, "Output/CRC_SNF_res_k.csv", row.names = F)

CRC_SNF_res_k_long <- CRC_SNF_res_k %>% pivot_longer(2:4, names_to = "Index", values_to = "Value")

CRC_SNF_res_k_long <- CRC_SNF_res_k_long %>% 
  mutate(Value = as.numeric(sprintf("%.3f", Value)),
         Index = factor(Index, levels = c("ARI", "NMI", "F1-score")))

g_CRC_SNF <- CRC_SNF_res_k_long %>% ggplot(aes(k_number, Value, fill = k_number)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(label = Value, hjust = 0.5, vjust = -0.2), size = 4) +
  scale_fill_manual(values = MetBrewer::met.brewer("VanGogh2", 5)) +
  theme_bw() + xlab("") + ylab("") + 
  facet_wrap(. ~ Index, ncol = 3, scales = "free_y", labeller = label_wrap_gen(multi_line = T)) + 
  theme(panel.border = element_rect(linewidth = 0.5),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.position = "top", 
        legend.text = element_text(size = 12),
        legend.title = element_blank())
g_CRC_SNF


########################################### 4.Survival analysis ##########################################

######### K=5
Clinical <- Clinical %>% 
  mutate(SNF_groups = CRC_SNF_res$`k=5`$clust) %>% 
  mutate(SNF_groups = case_when(SNF_groups == 1 ~ "SNF-I", 
                                SNF_groups == 2 ~ "SNF-II",
                                SNF_groups == 3 ~ "SNF-III", 
                                SNF_groups == 4 ~ "SNF-IV", 
                                SNF_groups == 5 ~ "SNF-V"))

Clinical <- Clinical %>% 
  mutate(SNF_groups = factor(SNF_groups, levels = c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V")))

########## SNF Cox regression 
Clinical <- Clinical %>% 
  mutate(SNF_groups_relevel = relevel(SNF_groups, ref = "SNF-V"))

coxph(Surv(OS_followtime, OS_status) ~ SNF_groups_relevel + Age + Stage, data = Clinical)
exp(confint(coxph(Surv(OS_followtime, OS_status) ~ SNF_groups_relevel + Age + Stage, data = Clinical)))

######### SNF OS KM
fit_SNF_OS <- survfit(Surv(OS_followtime, OS_status) ~ SNF_groups, data = Clinical)
SNF_OS <- ggsurvplot(fit_SNF_OS, data = Clinical,
                     pval = T, pval.coord = c(0, 0.4), legend = "none",
                     ylab = "Overall Survival", xlab = "Time in years", 
                     conf.int = F,
                     palette = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC"),
                     legend.lab = c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V"),
                     risk.table = T, 
                     tables.height = 0.3,
                     tables.y.text.col = T,
                     tables.y.text = T,
                     risk.table.pos = "out",
                     font.x = 18,
                     font.y = 18,
                     font.tickslab = 15,
                     font.caption = 15,
                     censor.shape = 3,
                     censor.size = 2)
SNF_OS

########## CMS Cox regression 
Clinical <- Clinical %>% 
  mutate(CMS_groups = factor(CMS_groups, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))) %>% 
  mutate(CMS_groups_relevel = relevel(CMS_groups, ref = "CMS4"))
  
coxph(Surv(OS_followtime, OS_status) ~ CMS_groups_relevel + Age + Stage, data = Clinical)
exp(confint(coxph(Surv(OS_followtime, OS_status) ~ CMS_groups_relevel + Age + Stage, data = Clinical)))

######### CMS OS KM
fit_CMS_OS <- survfit(Surv(OS_followtime, OS_status) ~ CMS_groups, data = Clinical)
CMS_OS <- ggsurvplot(fit_CMS_OS, data = Clinical,
                     pval = T, pval.coord = c(0, 0.4), legend = "none",
                     ylab = "Overall Survival", xlab = "Time in years", 
                     conf.int = F, 
                     palette = c("#00A087CC", "#3C5488CC", "#E18727FF", "#DC0000CC"),
                     legend.lab = c("CMS1", "CMS2", "CMS3", "CMS4"),
                     risk.table = T, 
                     tables.height = 0.3,
                     tables.y.text.col = T,
                     tables.y.text = T,
                     risk.table.pos = "out",
                     font.x = 18,
                     font.y = 18,
                     font.tickslab = 15,
                     font.caption = 15,
                     censor.shape = 3,
                     censor.size = 2)
CMS_OS

######### K=2, 3, 4, 6

K <- c(2, 3, 4, 6)
palettes <- c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC", "#3C5488CC")
legend_labs <- c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V", "SNF-VI")

lapply(K, function(ii) {
  
  Clinical <- Clinical %>% 
    mutate(SNF_groups = CRC_SNF_res[[ii-1]]$clust) %>% 
    mutate(SNF_groups = case_when(SNF_groups == 1 ~ "SNF-I", 
                                  SNF_groups == 2 ~ "SNF-II",
                                  SNF_groups == 3 ~ "SNF-III", 
                                  SNF_groups == 4 ~ "SNF-IV", 
                                  SNF_groups == 5 ~ "SNF-V", 
                                  SNF_groups == 6 ~ "SNF-VI"))
  write.csv(Clinical, paste0("Output/Subtype_info_K=", ii, ".csv"), row.names = F)
  
  ######### SNF OS KM
  fit_SNF_OS <- survfit(Surv(OS_followtime, OS_status) ~ SNF_groups, data = Clinical)
  SNF_OS <- ggsurvplot(fit_SNF_OS, data = Clinical,
                       pval = T, pval.coord = c(0, 0.4), legend = "none",
                       ylab = "Overall Survival", xlab = "Time in years", 
                       conf.int = F,
                       palette = palettes[1:ii],
                       legend.lab = legend_labs[1:ii],
                       risk.table = T, 
                       tables.height = 0.3,
                       tables.y.text.col = T,
                       tables.y.text = T,
                       risk.table.pos = "out",
                       font.x = 18,
                       font.y = 18,
                       font.tickslab = 15,
                       font.caption = 15,
                       censor.shape = 3,
                       censor.size = 2)
  
  pdf(paste0("Output/CRC/Figure/KM_SNF_OS_K=", ii, ".pdf"), 
      height = 5.5, width = 6, onefile = FALSE, family = "ArialMT")
  print(SNF_OS)
  dev.off()
  
})

