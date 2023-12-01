
rm(list = ls())

library(tidyverse)
library(circlize)
library(factoextra)
library(FactoMineR)
library(Rtsne)
library(patchwork)
library(ggalluvial)
library(tableone)


########################################## 1.Circos ########################################

Clinical <- read.csv("Output/Subtype_info_K=5.csv")
Clinical <- Clinical %>% 
  mutate(Stage = gsub("STAGE", "Stage", Stage))

myVars <- names(Clinical)[c(3, 4, 5, 14, 16, 18)]
catVars <- names(Clinical)[c(3, 5, 14, 16, 18)]

tab <- CreateTableOne(vars = myVars, data = Clinical, factorVars = catVars)
tab1 <- print(tab, showAllLevels = T, quote = F, noSpaces = T, printToggle = F)
tab1 <- tab1 %>% data.frame() %>% tibble::rownames_to_column(var = "var")
#writexl::write_xlsx(tab1, "Output/Baseline.xlsx")

###### K=5
table_cms <- table(Clinical$CMS_groups, Clinical$SNF_groups) %>% as.data.frame()
colors <- c(`CMS1` = "#00A087CC", `CMS2` = "#3C5488CC",
            `CMS3` = "#E18727FF", `CMS4` = "#DC0000CC",
            `SNF-I` = "#925E9FFF", `SNF-II` = "#0099B4FF", `SNF-III` = "#42B540FF", 
            `SNF-IV` = "#E18727FF", `SNF-V` = "#DC0000CC")

pdf("Output/Circos_SNF_CMS_K=5.pdf", height = 4, width = 4)
chordDiagram(table_cms, order = c("CMS1", "CMS2", "CMS3", "CMS4", 
                                  "SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V"),
             grid.col = colors, annotationTrack =  c("name", "grid"))
dev.off()


########## K=2, 3, 4, 6
K <- c(2, 3, 4, 6)
SNF_colors <- c(`SNF-I` = "#925E9FFF", `SNF-II` = "#0099B4FF", `SNF-III` = "#42B540FF", 
                `SNF-IV` = "#E18727FF", `SNF-V` = "#DC0000CC", `SNF-VI` =  "#3C5488CC")
SNFs <- c("SNF-I", "SNF-II", "SNF-III", "SNF-IV", "SNF-V", "SNF-VI")

lapply(K, function(ii) {

  Clinical <- read.csv(paste0("Output/Subtype_info_K=", ii, ".csv"))
  
  table_cms <- table(Clinical$CMS_groups, Clinical$SNF_groups) %>% as.data.frame()
  colors <- c(`CMS1` = "#00A087CC", `CMS2` = "#3C5488CC",
              `CMS3` = "#E18727FF", `CMS4` = "#DC0000CC", SNF_colors[1:ii])
  
  pdf(paste0("Output/Circos_SNF_CMS_K=", ii, ".pdf"), height = 4, width = 4)
  chordDiagram(table_cms, order = c("CMS1", "CMS2", "CMS3", "CMS4", SNFs[1:6]),
               grid.col = colors, annotationTrack =  c("name", "grid"))
  dev.off()
  circos.clear()
})


########################################## 2.t-SNE ########################################
SNF_group <- Clinical %>% select(c(1, 18)) %>% column_to_rownames(var = "Patient_ID") 
CRC_MAD <- readRDS("Output/CRC_MAD.rds")

########## tSNE

title <- names(CRC_MAD)

tsne <- lapply(1:length(title), function (ii) {
  tSNEdata <- CRC_MAD[[ii]]
  
  set.seed(666)
  tsne_res <- Rtsne(tSNEdata, perplexity = 50)
  
  data <- data.frame(tsne_res$Y, Clinical$SNF_groups)
  ggplot(data, aes(X1, X2, fill = Clinical$SNF_groups,
                   color = Clinical$SNF_groups)) + 
    geom_point(size = 1.5, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.1, linetype = 1) +
    xlab("t-SNE Dim1") + ylab("t-SNE Dim2") + ggtitle(title[ii]) + 
    scale_fill_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
    scale_color_manual(values = c("#925E9FFF", "#0099B4FF", "#42B540FF", "#E18727FF", "#DC0000CC")) +
    theme_classic() + 
    theme(title = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 20, color = "black"), 
          axis.text = element_text(size = 18, color = "black", face = "bold"),
          panel.grid = element_blank(), 
          legend.title = element_blank(),
          legend.position = c(0.9, 0.95))
})

gg_tsne <- tsne[[1]] + tsne[[2]] + tsne[[3]]
ggsave("Output/tSNE_K=5.pdf", gg_tsne, width = 15, height = 5)

