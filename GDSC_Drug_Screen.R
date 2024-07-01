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


sanger_dose_response <- read_csv("~/Documents/Doktor 2/Data/DepMap/sanger-dose-response.csv")
HNSC_Cell_Lines_DepMap_Risk_model <- read_delim("~/Nextcloud/Doktor/Sharing/Christoph Drug Screen/HNSC_Cell_Lines_DepMap_Risk_model.txt", "\t", escape_double = FALSE, trim_ws = TRUE)



pos_hr_drug_screen_df <- Pos_HR_Cell_Lines_DepMap_Risk_model %>%
  dplyr::select(Cell_ID, Risk_score) %>%
  left_join(., sanger_dose_response, by = c("Cell_ID" = "ARXSPAN_ID")) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID)


pos_hr_onlyone <- pos_hr_drug_screen_df %>%  
  group_by(`DRUG_ID`) %>%
  dplyr::summarize(n=n()) %>%
  dplyr::filter(n <=5 )

pos_hr_drug_screen_df <- pos_hr_drug_screen_df %>%  
  dplyr::filter(! (`DRUG_ID` %in% pos_hr_onlyone$`DRUG_ID`)) 

#reparing a data frame for further analysis

pos_hr_Drug_sensitivity_scores <- pos_hr_drug_screen_df %>%
  #pivot_longer(., cols = c(-Cell_ID, -Risk_score), names_to = "Drug", values_to = "DSS") %>%
  nest(-DRUG_ID) %>% 
  mutate(cor=map(data,~cor.test(.x$Risk_score, .x$IC50_PUBLISHED, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  unnest(data, .drop = T) %>%
  dplyr::select(DRUG_ID, BROAD_ID, Cell_ID, Risk_score, IC50_PUBLISHED, estimate, p.value, DRUG_NAME) %>%
  dplyr::rename(spearman = estimate)

write.table(pos_hr_Drug_sensitivity_scores, "pos_hr_DepMap_GDSC1_2_Drug_sensitivity_risk_score_spearman.txt", sep="\t", row.names=FALSE, na="")



pos_hr_plot_Drug_sensitivity_scores <- pos_hr_Drug_sensitivity_scores %>%
  mutate(class = case_when(p.value <= 0.005 ~ "Sig",
                           p.value > 0.005 ~ "Not Sig"))




volcanoplot_DSS_spearman <- ggplot(data = pos_hr_plot_Drug_sensitivity_scores, mapping = aes(x= spearman, y=-log(p.value), color = class)) +
  geom_point() +
  #geom_text_repel(mapping = aes(label = ifelse(p.value <= 0.05, DRUG_NAME, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
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
ggsave("pos_hr_DepMap_GDSC1_2_volcanoplot_DSS_spearman_005.tiff", plot = volcanoplot_DSS_spearman, width = 26, height= 26, units = "cm", dpi = 125)

cor_matrix <- Pos_HR_Cell_Lines_DepMap_Risk_model %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID) %>%
  dplyr::select(EGFR, `Risk Score`) %>%
  na.omit(.)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

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
ggsave("pos_hr_DepMap_GDSC_plot_EGFR_correlation.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


f <- as.data.frame(cor_matrix$r)
f$Other <-rownames(f)
write.table(f, "Spearman_rho_GDSC.txt", sep="\t", row.names=FALSE, na="")

f <- as.data.frame(cor_matrix$P)
f$Other <-rownames(f)
write.table(f, "Spearman_pvalue_GDSC.txt", sep="\t", row.names=FALSE, na="")

## Drugscreen mit allen 19 pos HR und rel pval
drug_screen_df <- HNSC_Cell_Lines_DepMap_Risk_model %>%
  dplyr::select(Cell_ID, Risk_score) %>%
  left_join(., sanger_dose_response, by = c("Cell_ID" = "ARXSPAN_ID")) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID)


onlyone <- drug_screen_df %>%  
  group_by(`DRUG_ID`) %>%
  dplyr::summarize(n=n()) %>%
  dplyr::filter(n <=5 )

drug_screen_df <- drug_screen_df %>%  
  dplyr::filter(! (`DRUG_ID` %in% onlyone$`DRUG_ID`)) 

#reparing a data frame for further analysis

Drug_sensitivity_scores <- drug_screen_df %>%
  #pivot_longer(., cols = c(-Cell_ID, -Risk_score), names_to = "Drug", values_to = "DSS") %>%
  nest(-DRUG_ID) %>% 
  mutate(cor=map(data,~cor.test(.x$Risk_score, .x$IC50_PUBLISHED, method = "sp"))) %>%
  mutate(tidied = map(cor, tidy)) %>% 
  unnest(tidied, .drop = T) %>%
  unnest(data, .drop = T) %>%
  dplyr::select(DRUG_ID, Cell_ID, Risk_score, IC50_PUBLISHED, estimate, p.value, DRUG_NAME) %>%
  dplyr::rename(spearman = estimate)

write.table(Drug_sensitivity_scores, "DepMap_GDSC1_2_Drug_sensitivity_risk_score_spearman.txt", sep="\t", row.names=FALSE, na="")



plot_Drug_sensitivity_scores <- Drug_sensitivity_scores %>%
  mutate(class = case_when(p.value <= 0.005 ~ "Sig",
                           p.value > 0.005 ~ "Not Sig"))




volcanoplot_DSS_spearman <- ggplot(data = plot_Drug_sensitivity_scores, mapping = aes(x= spearman, y=-log(p.value), color = class)) +
  geom_point() +
  #geom_text_repel(mapping = aes(label = ifelse(p.value <= 0.05, DRUG_NAME, ""),  hjust = 0.5, vjust = -1), color = "#999999", size = 4,angle = 0) +
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
ggsave("DepMap_GDSC1_2_volcanoplot_DSS_spearman_005.tiff", plot = volcanoplot_DSS_spearman, width = 26, height= 26, units = "cm", dpi = 125)

cor_matrix <- HNSC_Cell_Lines_DepMap_Risk_model %>%
  dplyr::rename(`Risk Score` = Risk_score) %>%
  dplyr::filter(Cell_ID %in% sanger_dose_response$ARXSPAN_ID) %>%
  dplyr::select(EGFR, `Risk Score`) %>%
  na.omit(.)

cor_matrix <- as.data.frame(cor_matrix) 
#rownames(cor_matrix) <- cor_matrix$Sample_ID
#cor_matrix$Sample_ID <- NULL
cor_matrix <- rcorr(as.matrix(cor_matrix), type = "spearman")

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
ggsave("DepMap_GDSC_plot_EGFR_correlation_HNSC.tiff", plot = k, width = 14, height= 14, units = "cm", dpi = 125)


f <- as.data.frame(cor_matrix$r)
f$Other <-rownames(f)
write.table(f, "Spearman_rho_GDSC.txt", sep="\t", row.names=FALSE, na="")

f <- as.data.frame(cor_matrix$P)
f$Other <-rownames(f)
write.table(f, "Spearman_pvalue_GDSC.txt", sep="\t", row.names=FALSE, na="")
