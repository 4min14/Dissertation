#Geneset generation and Clustering
library(tidyverse)
library(gridExtra)
library(futile.logger)
library(VennDiagram)
library(grid)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readr)
library(eulerr)
library(clusterProfiler)
library(enrichplot)


#Ihr braucht eine "gerankte" also eine geordnete liste für GSEA, mache das anhand der logFCs, hoffe ihr rafft das ganze, müsst natürlich paar sachen anpassen


#GSEA with Geneset==========
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


gsea_genelist <- EGFR_GSEA_ranklist$logFC_mean
names(gsea_genelist) <- (EGFR_GSEA_ranklist$Entrez_ID)

H_sig <- read.gmt("~/Desktop/Doktorarbeit/Data/GSEA/h.all.v7.4.entrez.gmt")
#c1_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c1.all.v7.2.entrez.gmt")
#c2_biocarta_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c2.cp.biocarta.v7.2.entrez.gmt")
#c2_kegg_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c2.cp.kegg.v7.2.entrez.gmt")
#c2_pid_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c2.cp.pid.v7.2.entrez.gmt")
#c2_reactome_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c2.cp.reactome.v7.2.entrez.gmt")
#c2_wikipathways_sig <- read.gmt(("~/Documents/Doktor 2/Data/MSigDb/c2.cp.wikipathways.v7.2.entrez.gmt"))
#c3tft_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c3.tft.v7.2.entrez.gmt")
#c5_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c5.go.v7.2.entrez.gmt")
#c5_bp_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c5.go.bp.v7.2.entrez.gmt")
#c5_cc_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c5.go.cc.v7.2.entrez.gmt")
#c5_mf_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c5.go.mf.v7.2.entrez.gmt")
c6_sig <- read.gmt("~/Desktop/Doktorarbeit/Data/GSEA/c6.all.v7.4.entrez.gmt")
#c7_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c7.all.v7.2.entrez.gmt") 
#c8_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c8.all.v7.2.entrez.gmt")

sig <- rbind(H_sig, c6_sig)#
#sig <- rbind(c2_kegg_sig,H_sig)
geneset_gsea  <- GSEA(geneList = gsea_genelist, exponent = 1, nPerm = 1000, minGSSize = 1,
                      maxGSSize = 1000, pvalueCutoff = 0.01, pAdjustMethod = "none",
                      TERM2GENE = sig, TERM2NAME = NA, verbose = TRUE, seed = FALSE,
                      by = "fgsea")

tiff(filename = "geneset_gsea.tiff", width = 25, height = 25, units = "cm", pointsize = 12, res = 75)
enrichplot::dotplot(geneset_gsea, x = "NES", showCategory=50,  orderBy = "x", color= "pvalue", font.size = 12)
dev.off()

#GO analysis=========
c5_bp_sig <- read.gmt("~/Documents/Doktor 2/Data/MSigDb/c5.bp.v7.1.entrez.gmt")

go_up <- Genes_regulation %>%
  dplyr::filter(reg == "Positive correlated Genes") %>%
  dplyr::select(Entrez_Gene_Id)

GO_DEG_up <- enricher(as.vector(go_up$Entrez_Gene_Id), pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = DEG_UCHL1_highvslow$ID,
                      minGSSize = 0, maxGSSize = 500, qvalueCutoff = 0.5, TERM2GENE = c5_bp_sig,
                      TERM2NAME = NA)
tiff(filename = "GO_DEG_up.tiff", width = 35, height = 10, units = "cm", pointsize = 12, res = 75)
barplot(GO_DEG_up, showCategory=5, font.size = 36)
dev.off()


go_down <- Genes_regulation %>%
  dplyr::filter(reg == "Negative correlated Genes") %>%
  dplyr::select(Entrez_Gene_Id)

GO_DEG_down <- enricher(as.vector(go_down$Entrez_Gene_Id), pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = DEG_UCHL1_highvslow$ID,
                        minGSSize = 0, maxGSSize = 500, qvalueCutoff = 0.5, TERM2GENE = c5_bp_sig,
                        TERM2NAME = NA)
tiff(filename = "GO_DEG_down.tiff", width = 35, height = 10, units = "cm", pointsize = 12, res = 75)
barplot(GO_DEG_down,x = 'Count' ,showCategory=5, font.size = 36)
dev.off()

temp <- DEG_UCHL1_highvslow %>%
  dplyr::filter(ID %in% Genes_regulation$Entrez_Gene_Id) %>%
  arrange(desc(logFC))
gsea_genelist <- temp$logFC
names(gsea_genelist) <- (temp$ID)

sig <- rbind(c2_biocarta_sig,c2_kegg_sig,H_sig,c2_pid_sig,c2_reactome_sig, c2_wiki_sig)#
#sig <- rbind(c2_kegg_sig,H_sig)
gs_gsea_H  <- GSEA(geneList = gsea_genelist, exponent = 1, nPerm = 1000, minGSSize = 1,
                   maxGSSize = 1000, pvalueCutoff = 0.01, pAdjustMethod = "none",
                   TERM2GENE =sig, TERM2NAME = NA, verbose = TRUE, seed = FALSE,
                   by = "fgsea")
tiff(filename = "gs_gsea_H.tiff", width = 25, height = 25, units = "cm", pointsize = 12, res = 75)
enrichplot::dotplot(gs_gsea_H, x = "NES", showCategory=50,  orderBy = "x", color= "pvalue", font.size = 12)
dev.off()


#Geneset for GSVA
gscoruchl_up <- Genes_regulation %>%
  dplyr::filter(., reg == "Positive correlation")

gscoruchl_down <- Genes_regulation %>%
  dplyr::filter(., reg == "Negative correlation")

gscoruchl <- list("Positive correlated Genes" = gscoruchl_up$Entrez_Gene_Id, "Negative correlated Genes" = gscoruchl_down$Entrez_Gene_Id)


#Validation Geneset

Hedgehog <- c2_kegg_sig %>%
  dplyr::filter(ont == "KEGG_HEDGEHOG_SIGNALING_PATHWAY")
IFNG <- H_sig %>%
  dplyr::filter(ont == "HALLMARK_INTERFERON_GAMMA_RESPONSE")
Hedgehog_H <-  H_sig %>%
  dplyr::filter(ont == "HALLMARK_HEDGEHOG_SIGNALING")

validation_gs <- list("Positive correlated Genes" = gscoruchl_up$Entrez_Gene_Id, "Negative correlated Genes" = gscoruchl_down$Entrez_Gene_Id)#, "KEGG_SHH" = Hedgehog$gene, "HALLMARK_IFNG" = IFNG$gene, "HALLMARK_SHH" = Hedgehog_H$gene )

