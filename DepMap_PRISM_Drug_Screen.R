#Prism Drug screen analysis=========
library(readxl)
library(readxl)
library(tidyverse)
library(tblhelpr)
library(broom)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(ggcorrplot)


#Calculating Risk Model 
# DATEN EINLESEN
sample_info <- read_csv("~/Desktop/Doktorarbeit/Data/sample_info.csv")
CCLE_expression <- read_csv("~/Desktop/Doktorarbeit/Data/CCLE_expression.csv") 
CCLE_expression <- CCLE_expression %>%
  dplyr::rename(Cell_ID = "...1")
sanger_dose_response <- read_csv("~/Desktop/Doktorarbeit/Data/sanger-dose-response.csv")
primary_screen_replicate_collapsed_logfold_change <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-logfold-change.csv")
primary_screen_replicate_collapsed_logfold_change <- primary_screen_replicate_collapsed_logfold_change %>%
  dplyr::rename(Cell_ID = "...1")

HNSC_Cell_Lines_DepMap <- sample_info %>%
  dplyr::filter(disease == "Head and Neck Cancer") %>%
  dplyr::filter(!str_detect(CCLE_Name, "SALIVARY_GLAND")) %>%
  dplyr::select(DepMap_ID) %>%
  dplyr::filter(DepMap_ID %in% sanger_dose_response$ARXSPAN_ID | DepMap_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  left_join(., CCLE_expression, by = c("DepMap_ID" = "Cell_ID")) %>%
  na.omit(.) %>% 
  transpose_tibble(., col_names = DepMap_ID, id_col = "Gene") %>%
  separate(.,Gene,into = c("Gene", "Entrez_Gene_ID"), sep = " \\(")
HNSC_Cell_Lines_DepMap$Entrez_Gene_ID <-   str_remove_all(HNSC_Cell_Lines_DepMap$Entrez_Gene_ID, pattern = "\\)")

HNSC_Cell_Lines_DepMap <- HNSC_Cell_Lines_DepMap %>%
  dplyr::select(-Entrez_Gene_ID) %>%
  transpose_tibble(., col_names = Gene, id_col = "Cell_ID")

write.table(HNSC_Cell_Lines_DepMap, "HNSC_Cell_Lines_DepMap_TPM_expression.txt", sep="\t", row.names=FALSE, na="")

HNSC_Cell_Lines_DepMap_Risk_model <- HNSC_Cell_Lines_DepMap %>%
  dplyr::mutate(., Risk_score = ZAP70 * (-0.0582152196382938) +  EFNB2 * (0.00470780696956776) + CCBE1 * (0.11219886435049) + ARRDC5 * (-0.267800991188706) + CALML5 * (-7.50310759764005e-05) + CCDC177 * (-0.0844511106749959))

write.table(HNSC_Cell_Lines_DepMap_Risk_model, "HNSC_Cell_Lines_DepMap_Risk_model.txt", sep="\t", row.names=FALSE, na="")


#Drug Screen


primary_screen_replicate_collapsed_logfold_change <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-logfold-change.csv") %>%
  dplyr::rename(Cell_ID = "...1")
primary_screen_replicate_collapsed_treatment_info <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-treatment-info.csv")

drug_screen_df <- HNSC_Cell_Lines_DepMap_Risk_model %>%
  dplyr::select(Cell_ID, Risk_score) %>%
  left_join(., primary_screen_replicate_collapsed_logfold_change, by = "Cell_ID") %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID)

#reparing a data frame for further analysis

Drug_sensitivity_scores <- drug_screen_df %>%
  pivot_longer(., cols = c(-Cell_ID, -Risk_score), names_to = "Drug", values_to = "DSS") %>%
  nest(-Drug) %>% 
  mutate(cor=map(data,~cor.test(.x$Risk_score, .x$DSS, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  dplyr::select(Drug, estimate, p.value) %>%
  dplyr::rename(spearman = estimate) %>%
  left_join(., primary_screen_replicate_collapsed_treatment_info, by=c("Drug" = "column_name"))

write.table(Drug_sensitivity_scores, "HNSC_DepMap_PRISM_Drug_sensitivity_risk_score_spearman.txt", sep="\t", row.names=FALSE, na="")



plot_Drug_sensitivity_scores <- Drug_sensitivity_scores %>%
  mutate(class = case_when(p.value <= 0.005 ~ "Sig",
                           p.value > 0.005 ~ "Not Sig"))
plot_Drug_sensitivity_scores$name <- str_to_sentence(plot_Drug_sensitivity_scores$name)



volcanoplot_DSS_spearman <- ggplot(data = plot_Drug_sensitivity_scores, mapping = aes(x= spearman, y=-log(p.value), color = class)) +
  geom_point() +
  #geom_text_repel(mapping = aes(label = ifelse(p.value <= 0.05, name, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
  #geom_text_repel(mapping = aes(label = ifelse(name %in% gse108061_drug_candidates$Drug, name, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
  scale_colour_manual(values = c(Sig = "red", `Not Sig` = "grey"), name = "Regulation") +
  #geom_vline(xintercept = 6, color ="gray") +
  #geom_vline(xintercept = -3, color ="gray") +
  geom_hline(yintercept = -log(0.005), color = "gray") + 
  #coord_cartesian(xlim = c(-50,100),ylim = c(0, 30)) + 
  xlab("Spearman's ρ") + 
  ylab("-log10 adj.P") +  
  #scale_x_log10() + 
  theme_classic() +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) + 
  #theme(axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
  theme(legend.position = "none") 
#ggtitle("Differential Expressed Genes") + 

volcanoplot_DSS_spearman
ggsave("DepMap_PRISM_volcanoplot_DSS_spearman_005.tiff", plot = volcanoplot_DSS_spearman, width = 26, height= 26, units = "cm", dpi = 125)


cor_matrix <- HNSC_Cell_Lines_DepMap_Risk_model %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  dplyr::select(EGFR, `Risk Score`)

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
ggsave("DepMap_PRISM_plot_EGFR_correlation_HNSC.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


f <- as.data.frame(cor_matrix$r)
f$Other <-rownames(f)
write.table(f, "Spearman_rho_PRISM.txt", sep="\t", row.names=FALSE, na="")

f <- as.data.frame(cor_matrix$P)
f$Other <-rownames(f)
write.table(f, "Spearman_pvalue_PRISM.txt", sep="\t", row.names=FALSE, na="")

### Drugscreen mit allten positiven HR models
#Calculating Risk Model 
# DATEN EINLESEN
sample_info <- read_csv("~/Desktop/Doktorarbeit/Data/sample_info.csv")
CCLE_expression <- read_csv("~/Desktop/Doktorarbeit/Data/CCLE_expression.csv") 
CCLE_expression <- CCLE_expression %>%
  dplyr::rename(Cell_ID = "...1")
sanger_dose_response <- read_csv("~/Desktop/Doktorarbeit/Data/sanger-dose-response.csv")
primary_screen_replicate_collapsed_logfold_change <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-logfold-change.csv")
primary_screen_replicate_collapsed_logfold_change <- primary_screen_replicate_collapsed_logfold_change %>%
  dplyr::rename(Cell_ID = "...1")

Pos_HR_Cell_Lines_DepMap <- sample_info %>%
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
Pos_HR_Cell_Lines_DepMap$Entrez_Gene_ID <-   str_remove_all(Pos_HR_Cell_Lines_DepMap$Entrez_Gene_ID, pattern = "\\)")

Pos_HR_Cell_Lines_DepMap <- Pos_HR_Cell_Lines_DepMap %>%
  dplyr::select(-Entrez_Gene_ID) %>%
  transpose_tibble(., col_names = Gene, id_col = "Cell_ID")

write.table(Pos_HR_Cell_Lines_DepMap, "Pos_HR_Cell_Lines_DepMap_TPM_expression.txt", sep="\t", row.names=FALSE, na="")

Pos_HR_Cell_Lines_DepMap_Risk_model <- Pos_HR_Cell_Lines_DepMap %>%
  dplyr::mutate(., Risk_score = ZAP70 * (-0.0582152196382938) +  EFNB2 * (0.00470780696956776) + CCBE1 * (0.11219886435049) + ARRDC5 * (-0.267800991188706) + CALML5 * (-7.50310759764005e-05) + CCDC177 * (-0.0844511106749959))

write.table(Pos_HR_Cell_Lines_DepMap_Risk_model, "Pos_HR_Cell_Lines_DepMap_Risk_model.txt", sep="\t", row.names=FALSE, na="")


#Drug Screen


primary_screen_replicate_collapsed_logfold_change <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-logfold-change.csv") %>%
  dplyr::rename(Cell_ID = "...1")
primary_screen_replicate_collapsed_treatment_info <- read_csv("~/Desktop/Doktorarbeit/Data/primary-screen-replicate-collapsed-treatment-info.csv")

pos_hr_drug_screen_df <- Pos_HR_Cell_Lines_DepMap_Risk_model %>%
  dplyr::select(Cell_ID, Risk_score) %>%
  left_join(., primary_screen_replicate_collapsed_logfold_change, by = "Cell_ID") %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID)

#reparing a data frame for further analysis

pos_hr_Drug_sensitivity_scores <- pos_hr_drug_screen_df %>%
  pivot_longer(., cols = c(-Cell_ID, -Risk_score), names_to = "Drug", values_to = "DSS") %>%
  nest(-Drug) %>% 
  mutate(cor=map(data,~cor.test(.x$Risk_score, .x$DSS, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  dplyr::select(Drug, estimate, p.value) %>%
  dplyr::rename(spearman = estimate) %>%
  left_join(., primary_screen_replicate_collapsed_treatment_info, by=c("Drug" = "column_name"))

write.table(pos_hr_Drug_sensitivity_scores, "pos_hr_DepMap_PRISM_Drug_sensitivity_risk_score_spearman.txt", sep="\t", row.names=FALSE, na="")



pos_hr_plot_Drug_sensitivity_scores <- pos_hr_Drug_sensitivity_scores %>%
  mutate(class = case_when(p.value <= 0.05 ~ "Sig",
                           p.value > 0.05 ~ "Not Sig"))
pos_hr_plot_Drug_sensitivity_scores$name <- str_to_sentence(pos_hr_plot_Drug_sensitivity_scores$name)



volcanoplot_DSS_spearman <- ggplot(data = pos_hr_plot_Drug_sensitivity_scores, mapping = aes(x= spearman, y=-log(p.value), color = class)) +
  geom_point() +
  #geom_text_repel(mapping = aes(label = ifelse(p.value <= 0.05, name, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
  #geom_text_repel(mapping = aes(label = ifelse(name %in% gse108061_drug_candidates$Drug, name, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
  scale_colour_manual(values = c(Sig = "red", `Not Sig` = "grey"), name = "Regulation") +
  #geom_vline(xintercept = 6, color ="gray") +
  #geom_vline(xintercept = -3, color ="gray") +
  geom_hline(yintercept = -log(0.05), color = "gray") + 
  #coord_cartesian(xlim = c(-50,100),ylim = c(0, 30)) + 
  xlab("Spearman's ρ") + 
  ylab("-log10 adj.P") +  
  #scale_x_log10() + 
  theme_classic() +
  theme(axis.title.y = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) + 
  #theme(axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
  theme(legend.position = "none") 
#ggtitle("Differential Expressed Genes") + 

volcanoplot_DSS_spearman
ggsave("pos_hr_DepMap_PRISM_volcanoplot_DSS_spearman_05.tiff", plot = volcanoplot_DSS_spearman, width = 26, height= 26, units = "cm", dpi = 125)


cor_matrix <- Pos_HR_Cell_Lines_DepMap_Risk_model %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% primary_screen_replicate_collapsed_logfold_change$Cell_ID) %>%
  dplyr::select(EGFR, `Risk Score`)

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
ggsave("pos_hr_DepMap_PRISM_plot_EGFR_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


f <- as.data.frame(cor_matrix$r)
f$Other <-rownames(f)
write.table(f, "Spearman_rho_PRISM.txt", sep="\t", row.names=FALSE, na="")

f <- as.data.frame(cor_matrix$P)
f$Other <-rownames(f)
write.table(f, "Spearman_pvalue_PRISM.txt", sep="\t", row.names=FALSE, na="")



#### ding mit den meds als spearman corrplot
ds_res_prism <- read_excel("~/Desktop/Doktorarbeit/Data/Drugscreen_posHR.xlsx", sheet = 3)
ds_res_prism <- as.data.frame(ds_res_prism)
ds_res_prism <- ds_res_prism %>%
  mutate(broad_id = str_sub(BROAD_ID, start = 1, end = 13))
ds_res_prism$BROAD_ID <- NULL
ds_res_gdsc <- read_excel("~/Desktop/Doktorarbeit/Data/Drugscreen_posHR.xlsx", sheet = 4)
ds_res_gdsc <- ds_res_gdsc %>%
  mutate(broad_id = str_sub(BROAD_ID, start = 1, end = 12)) 
ds_res_all1 <- left_join(ds_res_gdsc, ds_res_prism, by = "broad_id") %>%
  dplyr::filter(ds_res_gdsc$DRUG_ID %in% ds_res_prism$DRUG_ID) 
ds_res_all1 <- ds_res_all1 %>%  
  dplyr::filter(ds_res_all1$p.value.x < 0.05)
ds_res_all1 <- ds_res_all1 %>%
  dplyr::filter(ds_res_all1$p.value.y < 0.05)
ds_res_all1 <- ds_res_all1 %>%
  dplyr::select(DRUG_NAME.x, spearman_prism, spearman_gdsc)
rownames(ds_res_all1) <- ds_res_all1$DRUG_NAME.x
ds_res_all1 <- as.data.frame(ds_res_all1)
ds_res_all1$DRUG_NAME.x <- NULL
ds_res_all1 <- ds_res_all1 %>%
  rename(`PRISM` = spearman_prism) %>%
  rename(`GDSC` = spearman_gdsc)
attach(ds_res_all1)
ds_res_all1 <- ds_res_all1[order(PRISM),]
detach(ds_res_all1)
ds_res_all1 <- t(ds_res_all1)

p <- ggcorrplot(ds_res_all1,  
                method = "circle", 
                #colors = c("blue", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                #hc.order = F, 
                legend.title = "Spearman's ρ",
                #lab = F, 
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

p
ggsave("pos_hr_DS_res.tiff", plot = p, width = 14, height= 14, units = "cm", dpi = 125)


ds_res_prism <- read_excel("~/Desktop/Doktorarbeit/Data/Drugscreen_HNSC.xlsx", sheet = 1)
ds_res_gdsc <- read_excel("~/Desktop/Doktorarbeit/Data/Drugscreen_HNSC.xlsx", sheet = 2)
ds_res_all1 <- left_join(ds_res_gdsc, ds_res_prism, by = "DRUG_ID") %>%
  dplyr::filter(ds_res_gdsc$DRUG_ID %in% ds_res_prism$DRUG_ID) 
ds_res_all1 <- ds_res_all1 %>%  
  dplyr::filter(ds_res_all1$p.value.x < 0.05)
ds_res_all1 <- ds_res_all1 %>%
  dplyr::filter(ds_res_all1$p.value.y < 0.05)
ds_res_all1 <- ds_res_all1 %>%
  dplyr::select(DRUG_NAME.x, spearman_prism, spearman_gdsc)
rownames(ds_res_all1) <- ds_res_all1$DRUG_NAME.x
ds_res_all1 <- as.data.frame(ds_res_all1)
ds_res_all1$DRUG_NAME.x <- NULL
ds_res_all1 <- ds_res_all1 %>%
  rename(`PRISM` = spearman_prism) %>%
  rename(`GDSC` = spearman_gdsc)
attach(ds_res_all1)
ds_res_all1 <- ds_res_all1[order(PRISM),]
detach(ds_res_all1)
ds_res_all1 <- t(ds_res_all1)

p <- ggcorrplot(ds_res_all1,  
                method = "circle", 
                #colors = c("blue", mypal2[1]), #c(mypal2[2], "white", mypal2[1]),
                #hc.order = F, 
                legend.title = "Spearman's ρ",
                #lab = F, 
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

p
ggsave("HNSC_DS_res.tiff", plot = p, width = 14, height= 14, units = "cm", dpi = 125)


