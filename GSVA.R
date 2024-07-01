# gene set variation analysis codes in R
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(readr) #import txt
library(readxl) #import excel
library(tidyverse) # data manipulation
library(limma) # DEG
library(edgeR) # DEG
library(AnnotationDbi) # Converting Entrez Gene ID to Symbol
library(org.Hs.eg.db) # Converting Entrez Gene ID to Symbol
library(ComplexHeatmap) #Making the Heatmap
library(circlize) #colours
library(RColorBrewer) # colours
library(dplyr)
library(msigdbr)
library(clusterProfiler)


RiskScore<- read_excel("~/Desktop/PROMOTION/DEG/MALE/Lasso cox:spearman/OS:violin:spearman/all/new/RISK MALE OS UNBE.xlsx")

Pos_HR_Cell_Lines_GSVA <- sample_info %>%
  dplyr::filter(disease %in% c("Head and Neck Cancer","Cervical Cancer","Sarcoma","Bladder Cancer","Pancreatic Cancer","Endometrial/Uterine Cancer","Liver Cancer","Lymphoma","Breast Cancer", "Lung Cancer","Brain Cancer", "Bile Duct Cancer","Ovarian Cancer", "Thyroid Cancer" )) %>%
  dplyr::filter(!str_detect(CCLE_Name, "SALIVARY_GLAND")) %>%
  dplyr::filter(disease_subtype != "Small Cell Lung Cancer (SCLC)") %>%
  dplyr::filter(disease_subtype != "Non-Small Cell Lung Cancer (NSCLC), unspecified") %>%
  dplyr::filter(disease_subtype != "Astrocytoma") %>%
  dplyr::filter(disease_subtype != "Astrocytoma, anaplastic") %>%
  dplyr::filter(disease_subtype != "Medulloblastoma") %>%
  dplyr::filter(disease_subtype != "Meningioma") %>%
  dplyr::filter(disease_subtype != "Primitive Neuroectodermal Tumor (PNET)") %>%
  dplyr::filter(disease_subtype != "Hepatoblastoma") %>%
  dplyr::select(DepMap_ID) %>%
  dplyr::filter(DepMap_ID %in% sanger_dose_response$ARXSPAN_ID | DepMap_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  left_join(., CCLE_expression, by = c("DepMap_ID" = "Cell_ID")) %>%
  na.omit(.) %>% 
  transpose_tibble(., col_names = DepMap_ID, id_col = "Gene") %>%
  separate(.,Gene,into = c("Gene", "Entrez_Gene_ID"), sep = " \\(")
Pos_HR_Cell_Lines_GSVA$Entrez_Gene_ID <-   str_remove_all(Pos_HR_Cell_Lines_GSVA$Entrez_Gene_ID, pattern = "\\)")

Pos_HR_Cell_Lines_GSVA <- Pos_HR_Cell_Lines_GSVA %>%
  dplyr::select(-Gene) %>%
  transpose_tibble(., col_names = Entrez_Gene_ID, id_col = "Cell_ID")
rownames(Pos_HR_Cell_Lines_GSVA) <- Pos_HR_Cell_Lines_GSVA$Cell_ID
Pos_HR_Cell_Lines_GSVA <- as.data.frame(Pos_HR_Cell_Lines_GSVA)
Pos_HR_Cell_Lines_GSVA$Cell_ID <- NULL
Pos_HR_Cell_Lines_GSVA <- t(Pos_HR_Cell_Lines_GSVA)




expr <- Pos_HR_Cell_Lines_GSVA
expr <- as.matrix.data.frame(expr)




geneset_up <- genelist %>%
  dplyr::filter(Group == "up") %>%
  dplyr::select(GeneID_all)
#geneset_up <- as.list(geneset_up)
geneset_down <- genelist %>%
  dplyr::filter(Group == "down") %>%
  dplyr::select(GeneID_all)
#geneset_down <- as.list(geneset_down)


gsva_up <- gsva(expr, geneset_up, method="gsva", kcdf="Gaussian",min.sz=1,max.sz=Inf,tau=1,verbose=TRUE, parallel.sz=0)
gsva_down <- gsva(expr, geneset_down, method="gsva", kcdf="Gaussian",min.sz=1,max.sz=Inf,tau=1,verbose=TRUE, parallel.sz=0)


HNSC_Cell_Lines_GSVA <- sample_info %>%
  dplyr::filter(disease == "Head and Neck Cancer") %>%
  dplyr::filter(!str_detect(CCLE_Name, "SALIVARY_GLAND")) %>%
  dplyr::select(DepMap_ID) %>%
  dplyr::filter(DepMap_ID %in% sanger_dose_response$ARXSPAN_ID | DepMap_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  left_join(., CCLE_expression, by = c("DepMap_ID" = "Cell_ID")) %>%
  na.omit(.) %>% 
  transpose_tibble(., col_names = DepMap_ID, id_col = "Gene") %>%
  separate(.,Gene,into = c("Gene", "Entrez_Gene_ID"), sep = " \\(")
HNSC_Cell_Lines_GSVA$Entrez_Gene_ID <-   str_remove_all(HNSC_Cell_Lines_GSVA$Entrez_Gene_ID, pattern = "\\)")

HNSC_Cell_Lines_GSVA <- HNSC_Cell_Lines_GSVA %>%
  dplyr::select(-Gene) %>%
  transpose_tibble(., col_names = Entrez_Gene_ID, id_col = "Cell_ID")
rownames(HNSC_Cell_Lines_GSVA) <- HNSC_Cell_Lines_GSVA$Cell_ID
HNSC_Cell_Lines_GSVA <- as.data.frame(HNSC_Cell_Lines_GSVA)
HNSC_Cell_Lines_GSVA$Cell_ID <- NULL
HNSC_Cell_Lines_GSVA <- t(HNSC_Cell_Lines_GSVA)

expr_HNSC <- HNSC_Cell_Lines_GSVA

gsva_up_HNSC <- gsva(expr_HNSC, geneset_up, method="gsva", kcdf="Gaussian",min.sz=1,max.sz=Inf,tau=1,verbose=TRUE, parallel.sz=0)
gsva_down_HNSC <- gsva(expr_HNSC, geneset_down, method="gsva", kcdf="Gaussian",min.sz=1,max.sz=Inf,tau=1,verbose=TRUE, parallel.sz=0)

##spearman
gsva_down <- t(gsva_down)
gsva_down <- as.data.frame(gsva_down)
gsva_down$Cell_ID <- rownames(gsva_down)
rownames(gsva_down) <- NULL
gsva_down <- gsva_down %>%
  dplyr::rename(`gsva_down` = GeneID_all)
gsva_up <- t(gsva_up)
gsva_up <- as.data.frame(gsva_up)
gsva_up$Cell_ID <- rownames(gsva_up)
rownames(gsva_up) <- NULL
gsva_up <- gsva_up %>%
  dplyr::rename(`gsva_up` = GeneID_all)
gsva_spearman_all <- left_join(Pos_HR_Cell_Lines_DepMap_Risk_model, gsva_down, by = "Cell_ID")
gsva_spearman_all <- left_join(gsva_spearman_all, gsva_up, by = "Cell_ID")

cor_matrix <- gsva_spearman_all %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  dplyr::select( `Risk Score`, gsva_up, gsva_down)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

library(RColorBrewer)
mypal2 <- brewer.pal(9, "Set1")

k <- ggcorrplot(cor_matrix$r, 
                method = "circle", 
                colors = c("#6D9EC1", "white", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                hc.order = F, 
                legend.title = "Spearman's ρ",
                lab = F, 
                ggtheme =   theme_classic() +
                  theme(axis.title.y = element_text(size = 18),
                        axis.title.x = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16),
                        legend.text.align = 0,
                        legend.title.align = 0.5,
                        legend.box.just = "left",
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14))) + 
  
  guides(fill = guide_legend(label.position = "right", label.hjust = 1, title.position = "top"))

k
cor_matrix$r
ggsave("pos_hr_PRISM_plot_GSVA_RISK_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


gsva_down_HNSC <- t(gsva_down_HNSC)
gsva_down_HNSC <- as.data.frame(gsva_down_HNSC)
gsva_down_HNSC$Cell_ID <- rownames(gsva_down_HNSC)
rownames(gsva_down_HNSC) <- NULL
gsva_down_HNSC <- gsva_down_HNSC %>%
  dplyr::rename(`gsva_down` = GeneID_all)
gsva_up_HNSC <- t(gsva_up_HNSC)
gsva_up_HNSC <- as.data.frame(gsva_up_HNSC)
gsva_up_HNSC$Cell_ID <- rownames(gsva_up_HNSC)
rownames(gsva_up_HNSC) <- NULL
gsva_up_HNSC <- gsva_up_HNSC %>%
  dplyr::rename(`gsva_up` = GeneID_all)
gsva_spearman_HNSC <- left_join(Pos_HR_Cell_Lines_DepMap_Risk_model, gsva_down_HNSC, by = "Cell_ID")
gsva_spearman_HNSC <- left_join(gsva_spearman_HNSC, gsva_up_HNSC, by = "Cell_ID")

cor_matrix <- gsva_spearman_HNSC %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  dplyr::select( `Risk Score`, gsva_up, gsva_down)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

library(RColorBrewer)
mypal2 <- brewer.pal(9, "Set1")

k <- ggcorrplot(cor_matrix$r, 
                method = "circle", 
                colors = c("#6D9EC1", "white", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                hc.order = F, 
                legend.title = "Spearman's ρ",
                lab = F, 
                ggtheme =   theme_classic() +
                  theme(axis.title.y = element_text(size = 18),
                        axis.title.x = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16),
                        legend.text.align = 0,
                        legend.title.align = 0.5,
                        legend.box.just = "left",
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14))) + 
  
  guides(fill = guide_legend(label.position = "right", label.hjust = 1, title.position = "top"))

k
cor_matrix$r
ggsave("HNSC_PRISM_plot_GSVA_RISK_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)

###GDSC
cor_matrix <- gsva_spearman_HNSC %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID) %>%
  dplyr::select( `Risk Score`, gsva_up, gsva_down)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

library(RColorBrewer)
mypal2 <- brewer.pal(9, "Set1")

k <- ggcorrplot(cor_matrix$r, 
                method = "circle", 
                colors = c("#6D9EC1", "white", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                hc.order = F, 
                legend.title = "Spearman's ρ",
                lab = F, 
                ggtheme =   theme_classic() +
                  theme(axis.title.y = element_text(size = 18),
                        axis.title.x = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16),
                        legend.text.align = 0,
                        legend.title.align = 0.5,
                        legend.box.just = "left",
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14))) + 
  
  guides(fill = guide_legend(label.position = "right", label.hjust = 1, title.position = "top"))

k
cor_matrix$r
ggsave("HNSC_GDSC_plot_GSVA_RISK_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


cor_matrix <- gsva_spearman_all %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID) %>%
  dplyr::select( `Risk Score`, gsva_up, gsva_down)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

library(RColorBrewer)
mypal2 <- brewer.pal(9, "Set1")

k <- ggcorrplot(cor_matrix$r, 
                method = "circle", 
                colors = c("#6D9EC1", "white", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                hc.order = F, 
                legend.title = "Spearman's ρ",
                lab = F, 
                ggtheme =   theme_classic() +
                  theme(axis.title.y = element_text(size = 18),
                        axis.title.x = element_blank(),
                        legend.text=element_text(size=14),
                        legend.title=element_text(size=16),
                        legend.text.align = 0,
                        legend.title.align = 0.5,
                        legend.box.just = "left",
                        axis.text.x = element_text(size = 14),
                        axis.text.y = element_text(size = 14))) + 
  
  guides(fill = guide_legend(label.position = "right", label.hjust = 1, title.position = "top"))

k
cor_matrix$r
ggsave("pos_HR_GDSC_plot_GSVA_RISK_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)

