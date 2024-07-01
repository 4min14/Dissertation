library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MVisAGe)
library(GenVisR)
library(convaq)
library(ggpubr)
##CNV Daten einlesen (von firebrowse)
CNV_data <- read_delim("~/Downloads/gdac.broadinstitute.org_HNSC.Merge_snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3.2016012800.0.0/HNSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

risk_high <- riskall %>%
  dplyr::filter(risk == "high") %>%
  dplyr::select(id, risk)
risk_low <- riskall %>%
  dplyr::filter(risk == "low") %>%
  dplyr::select(id, risk)
write.table(risk_low, file = "risk_low.txt", quote = F)
write.table(risk_high, file = "risk_high.txt", quote = F)

CNV_data_low <- CNV_data %>%
  #dplyr::filter(Chromosome != "23") %>%
  dplyr::rename(chr = Chromosome, start = Start, end = End, seg_mean = `Segment_Mean`) %>%
  mutate(name = str_sub(Sample, start = 1, end = 12)) %>%
  filter(name %in% risk_low$id) %>%
  dplyr::select(name, chr, start, end, seg_mean)

CNV_data_high <- CNV_data %>%
  #dplyr::filter(Chromosome != "23") %>%
  dplyr::rename(chr = Chromosome, start = Start, end = End, seg_mean = `Segment_Mean`) %>%
  mutate(name = str_sub(Sample, start = 1, end = 12)) %>%
  filter(name %in% risk_high$id) %>%
  dplyr::select(name, chr, start, end, seg_mean)

write.table(CNV_data_high, "CNV_data_high.seg.txt", sep="\t", row.names=FALSE, na="")
write.table(CNV_data_low, "CNV_data_low.seg.txt", sep="\t", row.names=FALSE, na="")

HNSC_convaq_high <- CNV_data_high %>%
  mutate(type = case_when(seg_mean > 0.2 ~ "Gain",
                          seg_mean < -0.2 ~ "Loss",
                          seg_mean <= 0.2 & seg_mean >= -0.2 ~ "")) %>%
  dplyr::filter(type!="") %>%
  dplyr::select(-c(seg_mean))
write.table(HNSC_convaq_high, "HNSC_convaq_high.txt", sep="\t", row.names=FALSE, na="")

HNSC_convaq_low <- CNV_data_low %>%
  mutate(type = case_when(seg_mean > 0.2 ~ "Gain",
                          seg_mean < -0.2 ~ "Loss",
                          seg_mean <= 0.2 & seg_mean >= -0.2 ~ "")) %>%
  dplyr::filter(type!="") %>%
  dplyr::select(-c(seg_mean))
write.table(HNSC_convaq_UCHL1, "HNSC_convaq_low.txt", sep="\t", row.names=FALSE, na="")

HNSC_convaq_high_Gain <- HNSC_convaq_high %>%
  dplyr::filter(type == "Gain")
write.table(HNSC_convaq_high_Gain, "HNSC_convaq_high_Gain.txt", sep="\t", row.names=FALSE, na="")
HNSC_convaq_high_Loss <- HNSC_convaq_high %>%
  dplyr::filter(type == "Loss")
write.table(HNSC_convaq_high_Loss, "HNSC_convaq_high_Loss.txt", sep="\t", row.names=FALSE, na="")
HNSC_convaq_low_Gain <- HNSC_convaq_low %>%
  dplyr::filter(type == "Gain")
write.table(HNSC_convaq_low_Gain, "HNSC_convaq_low_Gain.txt", sep="\t", row.names=FALSE, na="")
HNSC_convaq_low_Loss <- HNSC_convaq_low %>%
  dplyr::filter(type == "Loss")
write.table(HNSC_convaq_low_Loss, "HNSC_convaq_low_Loss.txt", sep="\t", row.names=FALSE, na="")
  
Sig00001_CN_Gain_Genes <- read_excel("~/Desktop/Doktorarbeit/Data/CN_gain_genes.xlsx", sheet = 1)
Sig00001_CN_Loss_Genes <- read_excel("~/Desktop/Doktorarbeit/Data/CN_loss_genes.xlsx", sheet = 1)
Sig00001_CN_Gain <- read_excel("~/Desktop/Doktorarbeit/Data/CN_gain_result.xlsx", sheet = 1)
Sig00001_CN_Loss <- read_excel("~/Desktop/Doktorarbeit/Data/CN_loss_result.xlsx", sheet = 1)

convaq_pval_001_gain<-Sig00001_CN_Gain %>%
  dplyr::filter(pvalue < 0.01) %>%
  add_column(., "P_val < 0.01") %>%
  add_column(., "5") %>%
  dplyr::rename(name = `"P_val < 0.01"`, seg_mean = `"5"`) %>%
  dplyr::select(name, chr, start, end, seg_mean) %>%
  mutate(seg_mean = as.double(seg_mean)) %>%
  mutate(chr = as.double(chr))
write.table(convaq_pval_001_gain, "convaq_pval_001_gain.seg.txt", sep="\t", row.names=FALSE, na="")

convaq_pval_001_loss<-Sig00001_CN_Loss %>%
  dplyr::filter(pvalue < 0.01) %>%
  add_column(., "P_val < 0.01") %>%
  add_column(., "-5") %>%
  dplyr::rename(name = `"P_val < 0.01"`, seg_mean = `"-5"`) %>%
  dplyr::select(name, chr, start, end, seg_mean) %>%
  mutate(seg_mean = as.double(seg_mean)) %>%
  mutate(chr = as.double(chr))
write.table(convaq_pval_001_loss, "convaq_pval_001_loss.seg.txt", sep="\t", row.names=FALSE, na="")

convaq_pval_00001_gain<-Sig00001_CN_Gain %>%
  dplyr::filter(pvalue < 0.001) %>%
  add_column(., "P_val < 0.0001") %>%
  add_column(., "5") %>%
  dplyr::rename(name = `"P_val < 0.0001"`, seg_mean = `"5"`) %>%
  dplyr::select(name, chr, start, end, seg_mean) %>%
  mutate(seg_mean = as.double(seg_mean)) %>%
  mutate(chr = as.double(chr))
write.table(convaq_pval_00001_gain, "convaq_pval_00001_gain.seg.txt", sep="\t", row.names=FALSE, na="")

convaq_pval_00001_loss<-Sig00001_CN_Loss %>%
  dplyr::filter(pvalue < 0.001) %>%
  add_column(., "P_val < 0.0001") %>%
  add_column(., "-5") %>%
  dplyr::rename(name = `"P_val < 0.0001"`, seg_mean = `"-5"`) %>%
  dplyr::select(name, chr, start, end, seg_mean) %>%
  mutate(seg_mean = as.double(seg_mean)) %>%
  mutate(chr = as.double(chr))
write.table(convaq_pval_00001_loss, "convaq_pval_00001_loss.seg.txt", sep="\t", row.names=FALSE, na="")


