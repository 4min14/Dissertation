library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap) #Heatmap
library(circlize) #colours
library(RColorBrewer) # colours

#Daten erstellen
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)
crosstab_data <- TCGA_clinicaldata
risks <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1)
risks <- risks %>%
  dplyr::select(id, risk) %>%
  dplyr::rename(`Patient_ID` = `id`)
crosstab_data <- left_join(crosstab_data, risks, by = "Patient_ID")
crosstab_data <- left_join(crosstab_data, cluster, by = "Patient_ID")
crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, Age_61 = case_when(Age > 61 ~ "yes",
                                   Age <= 61 ~ "no"))
crosstab_data <- crosstab_data %>%  
  mutate(crosstab_data, cN = case_when(cN == "N0" ~ "cN0",
                                       cN != "N0" ~ "cN+"))

crosstab_data <- crosstab_data %>%  
  mutate(crosstab_data, pN = case_when(pN == "N0" ~ "pN0",
                                       pN != "N0" ~ "pN+"))

crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, cT = case_when(cT == "T1" ~ "cT1/2",
                                       cT == "T2" ~ "cT1/2",
                                       cT == "T3" ~ "cT3/4",
                                       cT == "TX" ~ "cT3/4",
                                       cT == "T4" ~ "cT3/4",
                                       cT == "T4a" ~ "cT3/4",
                                       cT == "T4b" ~ "cT3/4"))
crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, pT = case_when(pT == "T0" ~ "pT1/2",
                                              pT == "T1" ~ "pT1/2",
                                              pT == "T2" ~ "pT1/2",
                                              pT == "T3" ~ "pT3/4",
                                              pT == "TX" ~ "pT3/4",
                                              pT == "T4" ~ "pT3/4",
                                              pT == "T4a" ~ "pT3/4",
                                              pT == "T4b" ~ "pT3/4"))
crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, pathologic_stage = case_when(pathologic_stage == "Stage I" ~ "G1/2",
                                                     pathologic_stage == "Stage II" ~ "G1/2",
                                                     pathologic_stage == "Stage III" ~ "G3/4",
                                                     pathologic_stage == "Stage IVA" ~ "G3/4",
                                                     pathologic_stage == "Stage IVB" ~ "G3/4",
                                                     pathologic_stage == "Stage IVC" ~ "G3/4"))
crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, margin_status = case_when(margin_status == "Negative" ~ "R0",
                                                         margin_status == "Positive" ~ "R+",
                                                         margin_status == "Close" ~ "R+"))

crosstab_data <- crosstab_data %>%
  mutate(crosstab_data, Subsite = case_when(Subsite == "Lip" ~ "Oral cavity",
                                            Subsite == "Hypopharynx" ~ "Hypopharynx",
                                            Subsite == "Larynx" ~ "Larynx",
                                            Subsite == "Oral_cavity" ~ "Oral cavity",
                                            Subsite == "Oropharynx" ~ "Oropharynx"))

sn <- crosstab_data %>%
  dplyr::select(Cluster, Patient_ID, Gender, Age_61, Smoking, Alcohol, HPV16, Subsite, Radiation, PNI, ALI, pN, pT, cN, cT, margin_status, pathologic_stage, risk)




sn.tab <- table(sn$HPV16,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)
fisher.test(sn$HPV16, sn$Cluster)


chisq.test(sn.tab)

sn.tab <- table(sn$Gender,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Smoking,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Alcohol,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Radiation,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$EGFR,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Age_61,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Subsite,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)
fisher.test(sn$Subsite, sn$Cluster)

sn.tab <- table(sn$PNI,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$ALI,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$cN,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$pN,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$cT,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$pT,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$margin_status,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)
fisher.test(sn$margin_status, sn$Cluster)

sn.tab <- table(sn$pathologic_stage,sn$Cluster)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)
fisher.test(sn$pathologic_stage, sn$Cluster)

### Crosstab with riskmodel

sn.tab <- table(sn$Smoking,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Alcohol,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$HPV16,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Age_61,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Gender,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$Radiation,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$PNI,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$ALI,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)
fisher.test(sn$ALI, sn$risk)

sn.tab <- table(sn$cN,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$cT,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$pN,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$pT,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$margin_status,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)

sn.tab <- table(sn$pathologic_stage,sn$risk)
sn.tab
round(prop.table(sn.tab), 2)
round(prop.table(sn.tab,1),2)
round(prop.table(sn.tab,2),2)
#names(sn)
chisq.test(sn.tab)



