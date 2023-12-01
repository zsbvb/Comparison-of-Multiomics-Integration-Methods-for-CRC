
rm(list = ls())

library(tidyverse)
library(bigstatsr)
library(VennDiagram)

##### clinical data

Clinical <- read.table("Data/Clinical_patient.txt", 
                       header = T, 
                       sep = '\t',                 
                       quote = "",                              
                       comment.char = "#",         
                       fill = T)

Clinical <- Clinical %>% dplyr::select(c(1:3, 5:7, 20:22, 30, 31, 36, 37))
names(Clinical) <- c("Patient_ID", "Subtype", "Cancer_type", "Age", "Sex", "AJCC_stage", "M_stage", "N_stage", 
                     "T_stage", "OS_status", "OS_followtime", "PFS_status", "PFS_followtime")

Clinical <- Clinical %>% filter(!AJCC_stage == "") %>% drop_na(OS_followtime) %>%
  mutate(OS_followtime = OS_followtime/12, 
         PFS_followtime = PFS_followtime/12) %>% 
  filter(OS_followtime > 0 & PFS_followtime > 0)


#### DNA methylation

Methy <- read.table("Data/data_methylation_hm27_hm450_merged.txt", 
                    header = T, 
                    sep = '\t',                 
                    quote = "",                              
                    comment.char = "#",         
                    fill = T)

Methy <- Methy %>% 
  select(c(1, which(substr(names(Methy), 14, 15) %in% "01"))) %>%
  column_to_rownames(var = "ENTITY_STABLE_ID")

names(Methy) <- gsub("[.]", "-", names(Methy))

source("TCGA_operation.R")
id_filter <- tcgaReplicateFilter(tsb = names(Methy), full_barcode = T, analyte_target = "DNA")

names(Methy) <- substr(names(Methy), 1, 12)


########################## FPKM

mRNA_FPKM <- data.table::fread("Data/CORE_FPKM.csv") %>% 
  tibble::column_to_rownames(var = "V1")

mRNA_FPKM <- mRNA_FPKM[, which(substr(names(mRNA_FPKM), 14, 15) %in% "01")]

id_filter <- tcgaReplicateFilter(tsb = names(mRNA_FPKM), full_barcode = T, analyte_target = "RNA")
mRNA_FPKM <- mRNA_FPKM %>% dplyr::select(all_of(id_filter))

### log
FPKM_log <- log2(mRNA_FPKM + 1)
names(FPKM_log) <- substr(names(FPKM_log), 1, 12)

########################### protein

RPPA <- read.table("Data/data_rppa_zscores.txt",
                   header = T, 
                   sep = '\t',                 
                   quote = "",                              
                   comment.char = "#",         
                   fill = T)

RPPA <- RPPA %>% 
  select(c(1, which(substr(names(RPPA), 14, 15) %in% "01"))) %>%
  column_to_rownames(var = "Composite.Element.REF")

names(RPPA) <- gsub("[.]", "-", names(RPPA)) %>% substr(1, 12)




sample_id <- Reduce(intersect, list(Clinical$Patient_ID, names(FPKM_log), names(Methy), names(RPPA)))

set.seed(1111)   
venndata <- list(`Clinic\n554` = Clinical$Patient_ID, 
                 `Transcriptomics\n472` = names(FPKM_log), 
                 `Epigenomics\n591` = names(Methy), 
                 `Proteomics\n464` = names(RPPA))

venn.plot <- venn.diagram(
  venndata,
  filename = NULL, 
  fill = MetBrewer::met.brewer("Egypt", 4), 
  disable.logging = T,
  alpha = 0.50,          
  lwd = 1,           
  cat.col = MetBrewer::met.brewer("Egypt", 4), 
  cat.cex = 1.7,      
  cat.fontface = "bold", 
  cat.pos = 0,          
  cex = 1.7,
  fontface = "bold",  
  margin = 0)

Mergedata <- list(Clinical = Clinical %>% filter(Patient_ID %in% all_of(sample_id)),
                  Methy = Methy %>% select(all_of(sample_id)),
                  FPKM = FPKM_log %>% select(all_of(sample_id)), 
                  RPPA = RPPA %>% select(all_of(sample_id)))

save(Mergedata, file = "Data/Mergedata.Rdata")
