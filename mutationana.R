library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap) #Heatmap
library(circlize) #colours
library(RColorBrewer) # colours
library(countToFPKM)
library(glmnet)
library(survival)
library(maxstat)
library(survminer)

#####schauen ob die EGFR mutataion in den gruppen unterschiedlich ist####
mutationdata <- read_delim("~/Desktop/Doktorarbeit/Data/sample_matrix.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(mutationdata) <- mutationdata$`studyID:sampleId`
mutationdata <- as.data.frame(mutationdata)
mutationdata$`studyID:sampleId` <- NULL
mutationdata <- t(mutationdata)
colnames(mutationdata) <- gsub("-01","",colnames(mutationdata), fixed=TRUE)
colnames(mutationdata) <- gsub("hnsc_tcga_pan_can_atlas_2018:","",colnames(mutationdata), fixed=TRUE)
mutationdata <- t(mutationdata)
mutationdata <- as.data.frame(mutationdata)
mutationdata$Patient_ID <- rownames(mutationdata)
mutationdata <- left_join(mutationdata, risks, by = "Patient_ID")
chisq.test(mutationdata$EGFR, mutationdata$risk)

###### somatic mutation analyses#####
mat  <-  read_excel("~/Desktop/Doktorarbeit/Data/somatic_mutation_analysis.xlsx", sheet = 2)
mat <- as.data.frame(mat)
rownames(mat) <- mat$Gene
mat$Gene <- NULL
mat <- t(mat)
mat$Patient_ID <- rownames(mat)
mat <- mat %>%
  dplyr::select(Patient_ID, TP53, FAT1, CDKN2A, NOTCH1, PIK3CA, NSD1, CASP8, EP300, NFE2L2, FBXW7, HRAS)
risks <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1)
risks$Patient_ID <- risks$id
risks <- risks %>%
  dplyr::select(Patient_ID, risk)
risks_high <- risks
risks_high <- risks_high %>%
  filter(risks_high$risk == "high")
risks_low <- risks
risks_low <- risks_low %>%
  filter(risks_low$risk == "low")
mat_high <- mat %>%
  filter(mat$Patient_ID %in% risks_high$Patient_ID)
mat_high$Patient_ID <- NULL
mat_high <- t(mat_high)
mat_low <- mat %>%
  filter(mat$Patient_ID %in% risks_low$Patient_ID)
mat_low$Patient_ID <- NULL
mat_low <- t(mat_low)


col = c("MISS" = "blue", "TRU" = "red", "INF" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  MISS = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["MISS"], col = NA))
  },
  # bug red
  TRU = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["TRU"], col = NA))
  },
  # small green
  INF = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["INF"], col = NA))
  }
)
column_title_high = "OncoPrint for high risk"
heatmap_legend_param = list(title = "Mutation", at = c("MISS", "TRU", "INF"), 
                            labels = c("Missense Mutation", "Truncating mutation","Inframe Mutation"))

oncoPrint(mat_high,
          alter_fun = alter_fun, col = col, 
          column_title = column_title_high, heatmap_legend_param = heatmap_legend_param,row_order = 1:nrow(mat_high), remove_empty_columns = TRUE, remove_empty_rows = TRUE,) 

column_title_low = "OncoPrint for low risk"
oncoPrint(mat_low,
          alter_fun = alter_fun, col = col, 
          column_title = column_title_low, heatmap_legend_param = heatmap_legend_param,row_order = 1:nrow(mat_low), remove_empty_columns = TRUE, remove_empty_rows = TRUE,) 


###mutation altered data analyis
mutation_altered_data <- read_excel("~/Desktop/Doktorarbeit/Data/somatic_mutation_analysis.xlsx", sheet = 3)
rownames(mutation_altered_data) <- mutation_altered_data$`studyID:sampleId`
mutation_altered_data <- as.data.frame(mutation_altered_data)
mutation_altered_data$`studyID:sampleId` <- NULL
mutation_altered_data <- t(mutation_altered_data)
mutation_altered_data <- as.data.frame(mutation_altered_data)
colnames(mutation_altered_data) <- gsub("-01","",colnames(mutation_altered_data), fixed=TRUE)
colnames(mutation_altered_data) <- gsub("hnsc_tcga:","",colnames(mutation_altered_data), fixed=TRUE)
mutation_altered_data <- t(mutation_altered_data)
mutation_altered_data <- as.data.frame(mutation_altered_data)
mutation_altered_data$Patient_ID <- rownames(mutation_altered_data)

######HNSC riskscores aus dem cancerproof skript#####
riskscoresHNSC <- HNSC_risk %>%
  dplyr::select(bcr_patient_barcode, risk)
riskscoresHNSC <- as.data.frame(riskscoresHNSC)
riskscoresHNSC$Patient_ID <- riskscoresHNSC$bcr_patient_barcode
riskscoresHNSC$bcr_patient_barcode <- NULL

mutation_altered_data <- left_join(mutation_altered_data, riskscoresHNSC, by = "Patient_ID")
mutation_altered_data <- mutation_altered_data %>%
  dplyr::filter(!is.na(mutation_altered_data$risk))
rownames(mutation_altered_data) <- mutation_altered_data$Patient_ID
mutation_altered_data$Patient_ID <- NULL

##### Chisqure test isv crosstab machen mit jedem gene woooop#####
sn.tab <- table(mutation_altered_data$CASP8, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NOTCH1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$CDKN2A, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$FBXW7, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$GATA2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$HRAS, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$KRT17, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NEDD8, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NFE2L2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NRF1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NSD1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$PIK3CA, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$TP53, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ULK2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ZNF750, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$FAT1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR6C65, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$`HLA-B`, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$TGFBR2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$`HLA-A`, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$RAC1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$EP300, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR2M5, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$POM121L12, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$MAPK1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$EPHA2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$PTEN, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$DPPA2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$RHOA, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$AGTR1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$KIR3DL2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$B2M, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NPFFR2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$RASA1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$EYA1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$EPDR1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$DOK6, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR5D13, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$CD1E, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$CUL3, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$TRIM43, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ZNF99, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$C6, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ITGA8, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$KEAP1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ZNF676, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR5L1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$POTEG, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ZNF716, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$PDE10A, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$SMAD4, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$REG1A, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$LIN28B, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$CBWD1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$MEF2C, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR2M2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR8J3, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$PRB2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$PSG8, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$SPATA16, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NXPH1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$RARG, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$ALKAL1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$LCP1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$REG1B, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$AK5, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$FCRL4, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$H2BC8, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$H3C3, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$P2RY11, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$HBG2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$KHDRBS2, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$RSRC1, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$OR4K5, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)

sn.tab <- table(mutation_altered_data$NEK5, mutation_altered_data$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
chisq.test(sn.tab)
##### Barplot mit den mutationen
plotdata  <-  read_excel("~/Desktop/Doktorarbeit/Data/somatic_mutation_analysis.xlsx", sheet = 5)

plotdata %>% ggplot(aes(x = gene, fill = color)) +
  geom_bar(color = "black", position = "dodge") +
  scale_fill_manual(values = colors)
