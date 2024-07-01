library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)

TCGA_HNSC_GDCHarmonized_HTSeqFPKM <- read_delim("~/Desktop/Doktorarbeit/Data/TCGA_HNSC_GDCHarmonized_HTSeqFPKM.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

rnaseq_data <- read_delim("~/Desktop/Doktorarbeit/Data/TCGA_HNSC_GDCHarmonized_HTSeqCounts.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

HER2 <- tibble(Patient_ID = str_replace_all(colnames(TCGA_HNSC_GDCHarmonized_HTSeqFPKM),"-01$", ""),
               "2064" = TCGA_HNSC_GDCHarmonized_HTSeqFPKM %>% 
                 filter(Entrez_Gene_ID == 2064) %>%
                 t(.),
               Sample_ID = colnames(TCGA_HNSC_GDCHarmonized_HTSeqFPKM)) 

HER2 <- mutate(HER2, Group = case_when(`2064` >= quantile(`2064`, probs = 0.75, na.rm = T)~ "high", 
                                       `2064` > quantile(`2064`, probs = 0.25, na.rm = T) & `2064` < quantile(`2064`, probs = 0.75, na.rm = T) ~"moderate",
                                       `2064` <= quantile(`2064`, probs = 0.25, na.rm = T)~ "low"))

eset <- rnaseq_data %>% 
  as.data.frame(.)
rownames(eset) <- eset$Entrez_Gene_ID
eset$Entrez_Gene_ID <- NULL

temp_DESeq2 <- HER2 %>%
  dplyr::select(Sample_ID, Group)
temp_DESeq2 <- as.data.frame(colnames(eset)) %>%
  left_join(.,temp_DESeq2, by = c("colnames(eset)" = "Sample_ID")) %>%
  na.omit(.)

eset <- eset %>%
  dplyr::select(temp_DESeq2$`colnames(eset)`)

#eine Matrix mit den Samples als Spalte mit NAmen
(colData <- data.frame(row.names=colnames(eset), temp_DESeq2 = temp_DESeq2 ))
#(colData <- colData[-1])  (Wenn Dopplung stören sollte)




dds <- DESeqDataSetFromMatrix(countData = eset,
                              colData = colData,
                              design = ~temp_DESeq2.Group)
  
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

plotDispEsts(dds)

res <- results(dds, contrast=c("temp_DESeq2.Group", "high","low"))
head(res)
head(results(dds, tidy=TRUE))

resOrdered <- res[order(res$pvalue),]
resOrdered

DEG <-  as.data.frame(resOrdered)
HER2_DESeq2 <-  na.omit(DEG)

write.table(HER2_DESeq2, file="HER2_DESeq2_DGA.txt", sep="\t", col.names = NA, quote = F)
#Gensymbole einfügen as usual
HER2_DESeq2_DGA_GeneSymbol <- HER2_DESeq2
HER2_DESeq2_DGA_GeneSymbol$ID <- c(rownames(HER2_DESeq2_DGA_GeneSymbol))
HER2_DESeq2_DGA_GeneSymbol$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = HER2_DESeq2_DGA_GeneSymbol$ID, column= "SYMBOL", keytype =  "ENTREZID")
HER2_DESeq2_DGA_GeneSymbol <- HER2_DESeq2_DGA_GeneSymbol[-7]

write.table(HER2_DESeq2_DGA_GeneSymbol, "HER2_DESeq2_DGA_GeneSymbol.txt", sep="\t", col.names = NA, quote = F)

#Liste nur mit relevanten Genen
#DESeq2_Relevant <- EGFR_DESeq2_DGA_GeneSymbol %>% 
#filter(EGFR_DESeq2_DGA_GeneSymbol$Entrez_Gene_Id %in% genelist$GeneID_all)

volcanoplot_DESeq2 <- HER2_DESeq2 %>%
  mutate(class = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up",
                           log2FoldChange <= -1 & padj < 0.05 ~ "Down",
                           (log2FoldChange < 1 & log2FoldChange > -1) | padj >= 0.05 ~ "Not Sig"))

DEG_Volcanoplot_DESeq2 <- ggplot(data = volcanoplot_DESeq2, mapping = aes(x=log2FoldChange, y=-log(padj), colour = class)) +
  geom_point() +
  scale_colour_manual(values = c(Up = "red", Down = "blue", `Not Sig` = "grey")) +
  geom_vline(xintercept = 1, color ="pink") +
  geom_vline(xintercept = -1, color ="pink") +
  geom_hline(yintercept = -log(0.05), color = "pink") + 
  coord_cartesian(ylim = c(0,150),xlim = c(-10, 10)) + 
  xlab("log2FC") + 
  ylab("-log ad.P")

DEG_HER2_DESeq2_vp <- DEG_Volcanoplot_DESeq2 + theme_bw()

DEG_HER2_DESeq2_vp
ggsave("HER2_DESeq2_DGA.tiff", plot = DEG_HER2_DESeq2_vp)



