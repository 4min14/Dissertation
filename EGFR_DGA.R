library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer) 

#####CommonList_EGFR#####

#Die Ergebnisse einlesen
EGFR_DESeq2_DGA_GeneSymbol <- read_excel("~/Desktop/Doktorarbeit/EGFR_DESeq2/EGFR_DESeq2_DGA_GeneSymbol.xlsx", sheet = 1)
EGFR_DESeq2_res <- EGFR_DESeq2_DGA_GeneSymbol %>%
 dplyr::rename(Entrez_ID = ID)

EGFR_edgeR_DGA_GeneSymbol <- read_excel("~/Desktop/Doktorarbeit/EGFR_EdgeR/EGFR_EdgeR_DGA_GeneSymbol.xlsx", sheet = 1)
EGFR_edgeR_res <- EGFR_edgeR_DGA_GeneSymbol %>%
  dplyr::rename(Entrez_ID = ID)

EGFR_voom_DGA_GeneSymbol <- read_excel("~/Desktop/Doktorarbeit/EGFR_VOOM/EGFR_VOOM/EGFR_VOOM_DGA_sorted.xlsx", sheet = 1)  
EGFR_voom_res <- EGFR_voom_DGA_GeneSymbol %>%
  dplyr::rename(Entrez_ID = ID)

#Die wichtigen Gene filtern (LogFC >= 2 und p.adj <=0,05)

EGFR_DESeq2_relevant <- filter(EGFR_DESeq2_res,padj < 0.05)
EGFR_DESeq2_relevant <- filter(EGFR_DESeq2_relevant, logFC >= 1 | logFC <= -1)
EGFR_DESeq2_relevant <- mutate(EGFR_DESeq2_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                       logFC < 0 ~ "down"))

EGFR_edgeR_relevant <- filter(EGFR_edgeR_res, FDR < 0.05)
EGFR_edgeR_relevant <- filter(EGFR_edgeR_relevant, logFC >= 1 | logFC <= -1)
EGFR_edgeR_relevant <- mutate(EGFR_edgeR_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                     logFC < 0 ~ "down"))

EGFR_voom_relevant <- filter(EGFR_voom_res, adj.P.Val < 0.05)
EGFR_voom_relevant <- filter(EGFR_voom_relevant, logFC >= 1 | logFC <= -1)
EGFR_voom_relevant <- mutate(EGFR_voom_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                   logFC < 0 ~ "down"))

#Common Genelist erstellen
EGFR_relevant_ids <- EGFR_voom_relevant %>%
  dplyr::select(Entrez_ID, Group)
EGFR_common_relevant <- EGFR_relevant_ids %>%
  filter(EGFR_relevant_ids$Entrez_ID %in% EGFR_DESeq2_relevant$Entrez_ID)
EGFR_common_relevant <- EGFR_common_relevant %>% 
  filter(EGFR_common_relevant$Entrez_ID %in% EGFR_edgeR_relevant$Entrez_ID)


######GSEA vorbereiten, Rankliste erstellen#####

#1. Tabelle mit verschienden logFCs

EGFR_DESeq2_join <- EGFR_DESeq2_relevant %>%
  dplyr::select(Entrez_ID, logFC) %>%
  dplyr::rename(logFC_DESeq2 = logFC)

EGFR_voom_join <- EGFR_voom_relevant %>%
  dplyr::select(Entrez_ID, logFC) %>%
  dplyr::rename(logFC_voom = logFC)
  
EGFR_edgeR_join <- EGFR_edgeR_relevant %>%
  dplyr::select(Entrez_ID, logFC) %>%
  dplyr::rename(logFC_edgeR = logFC)

EGFR_GSEA_ranklist <- left_join(EGFR_common_relevant, EGFR_edgeR_join, by = "Entrez_ID")
EGFR_GSEA_ranklist <- left_join(EGFR_GSEA_ranklist, EGFR_DESeq2_join, by = "Entrez_ID")
EGFR_GSEA_ranklist <- left_join(EGFR_GSEA_ranklist, EGFR_voom_join, by = "Entrez_ID")
EGFR_GSEA_ranklist$Group <- NULL

#Mittelwert berechnen

EGFR_GSEA_ranklist <- mutate(EGFR_GSEA_ranklist, logFC_mean = (EGFR_GSEA_ranklist$logFC_edgeR + EGFR_GSEA_ranklist$logFC_DESeq2 + EGFR_GSEA_ranklist$logFC_voom)/ 3)
EGFR_GSEA_ranklist <- EGFR_GSEA_ranklist %>%
    dplyr::select(-logFC_edgeR, -logFC_voom, -logFC_DESeq2)
#Sortieren
EGFR_GSEA_ranklist <- EGFR_GSEA_ranklist[order(-EGFR_GSEA_ranklist$logFC_mean),]

#Eine Tabelle mit Genesymbols
EGFR_GSEA_ranklist_symbol <- EGFR_GSEA_ranklist
EGFR_GSEA_ranklist_symbol <- as.data.frame(EGFR_GSEA_ranklist_symbol)
rownames(EGFR_GSEA_ranklist_symbol) <- EGFR_GSEA_ranklist_symbol$Entrez_ID
EGFR_GSEA_ranklist_symbol$ID <- c(rownames(EGFR_GSEA_ranklist_symbol))
EGFR_GSEA_ranklist_symbol$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = EGFR_GSEA_ranklist_symbol$ID, column= "SYMBOL", keytype =  "ENTREZID")
EGFR_GSEA_ranklist_symbol <- EGFR_GSEA_ranklist_symbol %>%
  dplyr::select(Gene_Symbol, logFC_mean)
rownames(EGFR_GSEA_ranklist_symbol) <- NULL

#Ausgabe

write.table(EGFR_GSEA_ranklist, "EGFR_GSEA_ranklist.txt", sep="\t", col.names = NA, quote = F)

write.table(EGFR_GSEA_ranklist_symbol, "EGFR_GSEA_ranklist_symbol.txt", sep="\t", col.names = NA, quote = F)






######CommonList_EGFR_all

#Die Ergebnisse einlesen
EGFR_DESeq2_DGA_GeneSymbol_all <- read_excel("~/Desktop/Doktorarbeit/Data/DEG_ALLALL.xlsx", sheet = 3)
EGFR_DESeq2_all_res <- EGFR_DESeq2_DGA_GeneSymbol_all %>%
  dplyr::rename(Entrez_ID = ID)

EGFR_edgeR_DGA_GeneSymbol_all <- read_excel("~/Desktop/Doktorarbeit/Data/DEG_ALLALL.xlsx", sheet = 2)
EGFR_edgeR_all_res <- EGFR_edgeR_DGA_GeneSymbol_all %>%
  dplyr::rename(Entrez_ID = ID)

EGFR_voom_DGA_GeneSymbol_all <- read_excel("~/Desktop/Doktorarbeit/Data/DEG_ALLALL.xlsx", sheet = 1)  
EGFR_voom_all_res <- EGFR_voom_DGA_GeneSymbol_all %>%
  dplyr::rename(Entrez_ID = ID)

#Die wichtigen Gene filtern (LogFC >= 2 und p.adj <=0,05)

EGFR_DESeq2_all_relevant <- filter(EGFR_DESeq2_all_res,padj < 0.05)
EGFR_DESeq2_all_relevant <- filter(EGFR_DESeq2_all_relevant, logFC >= 1 | logFC <= -1)
EGFR_DESeq2_all_relevant <- mutate(EGFR_DESeq2_all_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                       logFC < 0 ~ "down"))

EGFR_edgeR_all_relevant <- filter(EGFR_edgeR_all_res, FDR < 0.05)
EGFR_edgeR_all_relevant <- filter(EGFR_edgeR_all_relevant, logFC >= 1 | logFC <= -1)
EGFR_edgeR_all_relevant <- mutate(EGFR_edgeR_all_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                     logFC < 0 ~ "down"))

EGFR_voom_all_relevant <- filter(EGFR_voom_all_res, adj.P.Val < 0.05)
EGFR_voom_all_relevant <- filter(EGFR_voom_all_relevant, logFC >= 1 | logFC <= -1)
EGFR_voom_all_relevant <- mutate(EGFR_voom_all_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                   logFC < 0 ~ "down"))

#Common Genelist erstellen
EGFR_all_relevant_ids <- EGFR_voom_all_relevant %>%
  dplyr::select(Entrez_ID, Gene_Symbol, Group)
EGFR_common_all_relevant <- EGFR_all_relevant_ids %>%
  filter(EGFR_all_relevant_ids$Entrez_ID %in% EGFR_DESeq2_all_relevant$Entrez_ID)
EGFR_common_all_relevant <- EGFR_common_all_relevant %>% 
  filter(EGFR_common_all_relevant$Entrez_ID %in% EGFR_edgeR_all_relevant$Entrez_ID)

write.table(EGFR_common_all_relevant, "EGFR_common_all_relevant.txt", sep="\t", col.names = NA, quote = F)
