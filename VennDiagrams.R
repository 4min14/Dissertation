library(VennDiagram)

## DATEN EINLESEN
Sig00001_CN_Gain_Genes <- read_excel("~/Desktop/Doktorarbeit/Data/CN_gain_genes.xlsx", sheet = 1)
Sig00001_CN_Loss_Genes <- read_excel("~/Desktop/Doktorarbeit/Data/CN_loss_genes.xlsx", sheet = 1)
genelist <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)

venn_gene_deg_up <- genelist %>%
  dplyr::filter(Group == "up") %>%
  dplyr::select(GeneID_all)
venn_gene_deg_down <- genelist %>%
  dplyr::filter(Group == "down") %>%
  dplyr::select(GeneID_all)
venn_gene_cnv_loss <- Sig00001_CN_Loss_Genes %>%
  dplyr::select("Entrez ID")
venn_gene_cnv_gain <- Sig00001_CN_Gain_Genes %>%
  dplyr::select("Entrez ID")

colors <-  c("#6b7fff", "#c3db0f", "#ff4059", "#2cff21")

venn.diagram(x = list(venn_gene_cnv_gain$`Entrez ID`, venn_gene_cnv_loss$`Entrez ID`, venn_gene_deg_up$GeneID_all, venn_gene_deg_down$GeneID_all) ,
             category.names = c("CNV Gain Genes", "CNV Loss Genes", "EGFR DEG Up", "EGFR DEG Down"),
             filename = 'CNV_venn.tiff',
             output=TRUE,
             imagetype="tiff", 
             scaled = FALSE,
             col = "black",
             fill = colors,
             cat.col = colors,
             cat.cex = 0.9,
             margin = 0.15
)

# Venndiagramm neu (mit nur 440 genen statt 441)
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
EGFR_DESeq2_relevant <- EGFR_DESeq2_relevant %>%
  filter(Entrez_ID != "6706")

EGFR_edgeR_relevant <- filter(EGFR_edgeR_res, FDR < 0.05)
EGFR_edgeR_relevant <- filter(EGFR_edgeR_relevant, logFC >= 1 | logFC <= -1)
EGFR_edgeR_relevant <- mutate(EGFR_edgeR_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                     logFC < 0 ~ "down"))
EGFR_edgeR_relevant <- EGFR_edgeR_relevant %>%
  filter(Entrez_ID != "6706")

EGFR_voom_relevant <- filter(EGFR_voom_res, adj.P.Val < 0.05)
EGFR_voom_relevant <- filter(EGFR_voom_relevant, logFC >= 1 | logFC <= -1)
EGFR_voom_relevant <- mutate(EGFR_voom_relevant, Group = case_when(logFC > 0 ~ "up",
                                                                   logFC < 0 ~ "down"))
EGFR_voom_relevant <- EGFR_voom_relevant %>%
  filter(Entrez_ID != "6706")

EGFR_DESeq2_up <- EGFR_DESeq2_relevant %>%
  filter(Group == "up")
EGFR_DESeq2_down <- EGFR_DESeq2_relevant %>%
  filter(Group == "down")

EGFR_edgeR_up <- EGFR_edgeR_relevant %>%
  filter(Group == "up")
EGFR_edgeR_down <- EGFR_edgeR_relevant %>%
  filter(Group == "down")

EGFR_voom_up <- EGFR_voom_relevant %>%
  filter(Group == "up")
EGFR_voom_down <- EGFR_voom_relevant %>%
  filter(Group == "down")

colors <-  c("#6b7fff", "#c3db0f", "#ff4059")

venn.diagram(x = list(EGFR_DESeq2_relevant$Entrez_ID, EGFR_edgeR_relevant$Entrez_ID, EGFR_voom_relevant$Entrez_ID) ,
             category.names = c("DESeq2", "edgeR", "Limma Voom"),
             filename = 'DEG_all_venn.tiff',
             output=TRUE,
             imagetype="tiff", 
             scaled = FALSE,
             col = "black",
             fill = colors,
             #cat.col = colors,
             cat.cex = 0.8,
             margin = 0.20
)

venn.diagram(x = list(EGFR_DESeq2_up$Entrez_ID, EGFR_edgeR_up$Entrez_ID, EGFR_voom_up$Entrez_ID) ,
             category.names = c("DESeq2_up", "edgeR_up", "Limma Voom_up"),
             filename = 'DEG_up_venn.tiff',
             output=TRUE,
             imagetype="tiff", 
             scaled = FALSE,
             col = "black",
             fill = F,
             #cat.col = colors,
             cat.cex = 1,
             margin = 0.20
)

venn.diagram(x = list(EGFR_DESeq2_down$Entrez_ID, EGFR_edgeR_down$Entrez_ID, EGFR_voom_down$Entrez_ID) ,
             category.names = c("DESeq2_down", "edgeR_down", "Limma Voom_down"),
             filename = 'DEG_down_venn.tiff',
             output=TRUE,
             imagetype="tiff", 
             scaled = FALSE,
             col = "black",
             fill = F,
             #cat.col = colors,
             cat.cex = 1,
             margin = 0.20
)


write.table(EGFR_DESeq2_relevant, "EGFR_DESeq2_all.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_edgeR_relevant, "EGFR_edgeR_all.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_voom_relevant, "EGFR_voom_all.txt", sep="\t", col.names = NA, quote = F)

write.table(EGFR_DESeq2_up, "EGFR_DESeq2_up.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_edgeR_up, "EGFR_edgeR_up.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_voom_up, "EGFR_voom_up.txt", sep="\t", col.names = NA, quote = F)

write.table(EGFR_DESeq2_down, "EGFR_DESeq2_down.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_edgeR_down, "EGFR_edgeR_down.txt", sep="\t", col.names = NA, quote = F)
write.table(EGFR_voom_down, "EGFR_voom_down.txt", sep="\t", col.names = NA, quote = F)

