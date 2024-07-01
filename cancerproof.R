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


### Daten einlesen
clinical_seinvater <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-CDR-SupplementalTableS1.xlsx", sheet = 1)
clinicaldata <- clinical_seinvater %>%
  dplyr::select(bcr_patient_barcode, OS, OS.time)



#####ACC##### 
## 1. Countdaten einlesen

ACC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/ACC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(ACC_counts) <- ACC_counts$Entrez_Gene_ID
ACC_counts <- as.data.frame(ACC_counts)
ACC_counts$Entrez_Gene_ID <- NULL
colnames(ACC_counts) <- gsub("-01","",colnames(ACC_counts), fixed=TRUE)
ACC_counts <- t(ACC_counts)
ACC_counts <- as.data.frame(ACC_counts)
ACC_counts$riskscore <- ACC_counts$'7535' * (-0.0582152196382938) + ACC_counts$'494514' * 0.013184642798603 + ACC_counts$'1948' * 0.00470780696956776 + ACC_counts$'147372' * 0.11219886435049 + ACC_counts$'645432' * -0.267800991188706 + ACC_counts$'51806' * -7.50310759764005e-05 + ACC_counts$'56936' * -0.0844511106749959
ACC_counts$bcr_patient_barcode <- rownames(ACC_counts)
ACC_counts <- ACC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

ACC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% ACC_counts$bcr_patient_barcode)

ACC_clinicdata <- ACC_clinicdata[match(rownames(ACC_counts),ACC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

ACC_risk <- left_join(ACC_clinicdata, ACC_counts, by = "bcr_patient_barcode")
ACC_risk$OS <- as.numeric(ACC_risk$OS)
ACC_risk$OS.time <- as.numeric(ACC_risk$OS.time)

stat <- maxstat.test(Surv(ACC_risk$OS.time , ACC_risk$OS) ~ ACC_risk$riskscore, data = ACC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(ACC_risk$riskscore > cutoff,"high","low"))
ACC_risk$risk <- risk
## Kaplanmeier
ACC_risk$OS.time <- ACC_risk$OS.time / 365 * 12
surv_obj <- Surv(ACC_risk$OS.time , ACC_risk$OS)
fit <- survfit(surv_obj~risk, data = ACC_risk) 
summary(fit)
ggsurvplot(fit, 
     title = "ACC",
     xlab = "Survival in months",
     xlim = c(0, 60),
     break.x.by = 12,
     ylab = "Survival Rate",
     mark.time = T,
     risk.table = T,
     risk.table.title = "",
     risk.table.height = .25,
     palette =c("red1","blue"),
     pval = T,
     legend.title = "Risk",
     legend.labs = c("high", "low"))
## univariate cox
risk_cancerproof <- ACC_risk 
risk_cancerproof <- risk_cancerproof %>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                 risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)

#####BLCA####
BLCA_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/BLCA_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(BLCA_counts) <- BLCA_counts$Entrez_Gene_ID
BLCA_counts <- as.data.frame(BLCA_counts)
BLCA_counts$Entrez_Gene_ID <- NULL
colnames(BLCA_counts) <- gsub("-01","",colnames(BLCA_counts), fixed=TRUE)
BLCA_counts <- t(BLCA_counts)
BLCA_counts <- as.data.frame(BLCA_counts)
BLCA_counts$riskscore <- BLCA_counts$'7535' * (-0.0582152196382938) + BLCA_counts$'494514' * 0.013184642798603 + BLCA_counts$'1948' * 0.00470780696956776 + BLCA_counts$'147372' * 0.11219886435049 + BLCA_counts$'645432' * -0.267800991188706 + BLCA_counts$'51806' * -7.50310759764005e-05 + BLCA_counts$'56936' * -0.0844511106749959
BLCA_counts$bcr_patient_barcode <- rownames(BLCA_counts)
BLCA_counts <- BLCA_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

BLCA_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% BLCA_counts$bcr_patient_barcode)

BLCA_clinicdata <- BLCA_clinicdata[match(rownames(BLCA_counts),BLCA_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

BLCA_risk <- left_join(BLCA_clinicdata, BLCA_counts, by = "bcr_patient_barcode")
BLCA_risk$OS <- as.numeric(BLCA_risk$OS)
BLCA_risk$OS.time <- as.numeric(BLCA_risk$OS.time)

stat <- maxstat.test(Surv(BLCA_risk$OS.time , BLCA_risk$OS) ~ BLCA_risk$riskscore, data = BLCA_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(BLCA_risk$riskscore > cutoff,"high","low"))
BLCA_risk$risk <- risk
## Kaplanmeier
BLCA_risk$OS.time <- BLCA_risk$OS.time / 365 * 12
surv_obj <- Surv(BLCA_risk$OS.time , BLCA_risk$OS)
fit <- survfit(surv_obj~risk, data = BLCA_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "BLCA",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           conf.int = T,
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univariate cox
risk_cancerproof <- BLCA_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####BRCA####
BRCA_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/BRCA_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(BRCA_counts) <- BRCA_counts$Entrez_Gene_ID
BRCA_counts <- as.data.frame(BRCA_counts)
BRCA_counts$Entrez_Gene_ID <- NULL
colnames(BRCA_counts) <- gsub("-01","",colnames(BRCA_counts), fixed=TRUE)
BRCA_counts <- t(BRCA_counts)
BRCA_counts <- as.data.frame(BRCA_counts)
BRCA_counts$riskscore <- BRCA_counts$'7535' * (-0.0582152196382938) + BRCA_counts$'494514' * 0.013184642798603 + BRCA_counts$'1948' * 0.00470780696956776 + BRCA_counts$'147372' * 0.11219886435049 + BRCA_counts$'645432' * -0.267800991188706 + BRCA_counts$'51806' * -7.50310759764005e-05 + BRCA_counts$'56936' * -0.0844511106749959
BRCA_counts$bcr_patient_barcode <- rownames(BRCA_counts)
BRCA_counts <- BRCA_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

BRCA_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% BRCA_counts$bcr_patient_barcode)

BRCA_clinicdata <- BRCA_clinicdata[match(rownames(BRCA_counts),BRCA_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

BRCA_risk <- left_join(BRCA_clinicdata, BRCA_counts, by = "bcr_patient_barcode")
BRCA_risk$OS <- as.numeric(BRCA_risk$OS)
BRCA_risk$OS.time <- as.numeric(BRCA_risk$OS.time)

stat <- maxstat.test(Surv(BRCA_risk$OS.time , BRCA_risk$OS) ~ BRCA_risk$riskscore, data = BRCA_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(BRCA_risk$riskscore > cutoff,"high","low"))
BRCA_risk$risk <- risk
## Kaplanmeier
BRCA_risk$OS.time <- BRCA_risk$OS.time / 365 * 12
surv_obj <- Surv(BRCA_risk$OS.time , BRCA_risk$OS)
fit <- survfit(surv_obj~risk, data = BRCA_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "BRCA",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- BRCA_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####CESC####
CESC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/CESC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(CESC_counts) <- CESC_counts$Entrez_Gene_ID
CESC_counts <- as.data.frame(CESC_counts)
CESC_counts$Entrez_Gene_ID <- NULL
colnames(CESC_counts) <- gsub("-01","",colnames(CESC_counts), fixed=TRUE)
CESC_counts <- t(CESC_counts)
CESC_counts <- as.data.frame(CESC_counts)
CESC_counts$riskscore <- CESC_counts$'7535' * (-0.0582152196382938) + CESC_counts$'494514' * 0.013184642798603 + CESC_counts$'1948' * 0.00470780696956776 + CESC_counts$'147372' * 0.11219886435049 + CESC_counts$'645432' * -0.267800991188706 + CESC_counts$'51806' * -7.50310759764005e-05 + CESC_counts$'56936' * -0.0844511106749959
CESC_counts$bcr_patient_barcode <- rownames(CESC_counts)
CESC_counts <- CESC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

CESC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% CESC_counts$bcr_patient_barcode)

CESC_clinicdata <- CESC_clinicdata[match(rownames(CESC_counts),CESC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

CESC_risk <- left_join(CESC_clinicdata, CESC_counts, by = "bcr_patient_barcode")
CESC_risk$OS <- as.numeric(CESC_risk$OS)
CESC_risk$OS.time <- as.numeric(CESC_risk$OS.time)

stat <- maxstat.test(Surv(CESC_risk$OS.time , CESC_risk$OS) ~ CESC_risk$riskscore, data = CESC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(CESC_risk$riskscore > cutoff,"high","low"))
CESC_risk$risk <- risk
## Kaplanmeier
CESC_risk$OS.time <- CESC_risk$OS.time / 365 * 12
surv_obj <- Surv(CESC_risk$OS.time , CESC_risk$OS)
fit <- survfit(surv_obj~risk, data = CESC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "CESC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           conf.int = T,
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- CESC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####CHOL####
CHOL_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/CHOL_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(CHOL_counts) <- CHOL_counts$Entrez_Gene_ID
CHOL_counts <- as.data.frame(CHOL_counts)
CHOL_counts$Entrez_Gene_ID <- NULL
colnames(CHOL_counts) <- gsub("-01","",colnames(CHOL_counts), fixed=TRUE)
CHOL_counts <- t(CHOL_counts)
CHOL_counts <- as.data.frame(CHOL_counts)
CHOL_counts$riskscore <- CHOL_counts$'7535' * (-0.0582152196382938) + CHOL_counts$'494514' * 0.013184642798603 + CHOL_counts$'1948' * 0.00470780696956776 + CHOL_counts$'147372' * 0.11219886435049 + CHOL_counts$'645432' * -0.267800991188706 + CHOL_counts$'51806' * -7.50310759764005e-05 + CHOL_counts$'56936' * -0.0844511106749959
CHOL_counts$bcr_patient_barcode <- rownames(CHOL_counts)
CHOL_counts <- CHOL_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

CHOL_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% CHOL_counts$bcr_patient_barcode)

CHOL_clinicdata <- CHOL_clinicdata[match(rownames(CHOL_counts),CHOL_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

CHOL_risk <- left_join(CHOL_clinicdata, CHOL_counts, by = "bcr_patient_barcode")
CHOL_risk$OS <- as.numeric(CHOL_risk$OS)
CHOL_risk$OS.time <- as.numeric(CHOL_risk$OS.time)

stat <- maxstat.test(Surv(CHOL_risk$OS.time , CHOL_risk$OS) ~ CHOL_risk$riskscore, data = CHOL_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(CHOL_risk$riskscore > cutoff,"high","low"))
CHOL_risk$risk <- risk
## Kaplanmeier
CHOL_risk$OS.time <- CHOL_risk$OS.time / 365 * 12
surv_obj <- Surv(CHOL_risk$OS.time , CHOL_risk$OS)
fit <- survfit(surv_obj~risk, data = CHOL_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "CHOL",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- CHOL_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####COAD####
COAD_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/COAD_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(COAD_counts) <- COAD_counts$Entrez_Gene_ID
COAD_counts <- as.data.frame(COAD_counts)
COAD_counts$Entrez_Gene_ID <- NULL
colnames(COAD_counts) <- gsub("-01","",colnames(COAD_counts), fixed=TRUE)
COAD_counts <- t(COAD_counts)
COAD_counts <- as.data.frame(COAD_counts)
COAD_counts$riskscore <- COAD_counts$'7535' * (-0.0582152196382938) + COAD_counts$'494514' * 0.013184642798603 + COAD_counts$'1948' * 0.00470780696956776 + COAD_counts$'147372' * 0.11219886435049 + COAD_counts$'645432' * -0.267800991188706 + COAD_counts$'51806' * -7.50310759764005e-05 + COAD_counts$'56936' * -0.0844511106749959
COAD_counts$bcr_patient_barcode <- rownames(COAD_counts)
COAD_counts <- COAD_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

COAD_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% COAD_counts$bcr_patient_barcode)

COAD_clinicdata <- COAD_clinicdata[match(rownames(COAD_counts),COAD_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

COAD_risk <- left_join(COAD_clinicdata, COAD_counts, by = "bcr_patient_barcode")
COAD_risk$OS <- as.numeric(COAD_risk$OS)
COAD_risk$OS.time <- as.numeric(COAD_risk$OS.time)

stat <- maxstat.test(Surv(COAD_risk$OS.time , COAD_risk$OS) ~ COAD_risk$riskscore, data = COAD_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(COAD_risk$riskscore > cutoff,"high","low"))
COAD_risk$risk <- risk
## Kaplanmeier
COAD_risk$OS.time <- COAD_risk$OS.time / 365 * 12
surv_obj <- Surv(COAD_risk$OS.time , COAD_risk$OS)
fit <- survfit(surv_obj~risk, data = COAD_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "COAD",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- COAD_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####DLBC####
DLBC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/DLBC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(DLBC_counts) <- DLBC_counts$Entrez_Gene_ID
DLBC_counts <- as.data.frame(DLBC_counts)
DLBC_counts$Entrez_Gene_ID <- NULL
colnames(DLBC_counts) <- gsub("-01","",colnames(DLBC_counts), fixed=TRUE)
DLBC_counts <- t(DLBC_counts)
DLBC_counts <- as.data.frame(DLBC_counts)
DLBC_counts$riskscore <- DLBC_counts$'7535' * (-0.0582152196382938) + DLBC_counts$'494514' * 0.013184642798603 + DLBC_counts$'1948' * 0.00470780696956776 + DLBC_counts$'147372' * 0.11219886435049 + DLBC_counts$'645432' * -0.267800991188706 + DLBC_counts$'51806' * -7.50310759764005e-05 + DLBC_counts$'56936' * -0.0844511106749959
DLBC_counts$bcr_patient_barcode <- rownames(DLBC_counts)
DLBC_counts <- DLBC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

DLBC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% DLBC_counts$bcr_patient_barcode)

DLBC_clinicdata <- DLBC_clinicdata[match(rownames(DLBC_counts),DLBC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

DLBC_risk <- left_join(DLBC_clinicdata, DLBC_counts, by = "bcr_patient_barcode")
DLBC_risk$OS <- as.numeric(DLBC_risk$OS)
DLBC_risk$OS.time <- as.numeric(DLBC_risk$OS.time)

stat <- maxstat.test(Surv(DLBC_risk$OS.time , DLBC_risk$OS) ~ DLBC_risk$riskscore, data = DLBC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(DLBC_risk$riskscore > cutoff,"high","low"))
DLBC_risk$risk <- risk
## Kaplanmeier
DLBC_risk$OS.time <- DLBC_risk$OS.time / 365 * 12
surv_obj <- Surv(DLBC_risk$OS.time , DLBC_risk$OS)
fit <- survfit(surv_obj~risk, data = DLBC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "DLBC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- DLBC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####ESCA####
ESCA_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/ESCA_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(ESCA_counts) <- ESCA_counts$Entrez_Gene_ID
ESCA_counts <- as.data.frame(ESCA_counts)
ESCA_counts$Entrez_Gene_ID <- NULL
colnames(ESCA_counts) <- gsub("-01","",colnames(ESCA_counts), fixed=TRUE)
ESCA_counts <- t(ESCA_counts)
ESCA_counts <- as.data.frame(ESCA_counts)
ESCA_counts$riskscore <- ESCA_counts$'7535' * (-0.0582152196382938) + ESCA_counts$'494514' * 0.013184642798603 + ESCA_counts$'1948' * 0.00470780696956776 + ESCA_counts$'147372' * 0.11219886435049 + ESCA_counts$'645432' * -0.267800991188706 + ESCA_counts$'51806' * -7.50310759764005e-05 + ESCA_counts$'56936' * -0.0844511106749959
ESCA_counts$bcr_patient_barcode <- rownames(ESCA_counts)
ESCA_counts <- ESCA_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

ESCA_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% ESCA_counts$bcr_patient_barcode)

ESCA_clinicdata <- ESCA_clinicdata[match(rownames(ESCA_counts),ESCA_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

ESCA_risk <- left_join(ESCA_clinicdata, ESCA_counts, by = "bcr_patient_barcode")
ESCA_risk$OS <- as.numeric(ESCA_risk$OS)
ESCA_risk$OS.time <- as.numeric(ESCA_risk$OS.time)

stat <- maxstat.test(Surv(ESCA_risk$OS.time , ESCA_risk$OS) ~ ESCA_risk$riskscore, data = ESCA_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(ESCA_risk$riskscore > cutoff,"high","low"))
ESCA_risk$risk <- risk
## Kaplanmeier
ESCA_risk$OS.time <- ESCA_risk$OS.time / 365 * 12
surv_obj <- Surv(ESCA_risk$OS.time , ESCA_risk$OS)
fit <- survfit(surv_obj~risk, data = ESCA_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "ESCA",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- ESCA_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####GBM####
GBM_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/GBM_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(GBM_counts) <- GBM_counts$Entrez_Gene_ID
GBM_counts <- as.data.frame(GBM_counts)
GBM_counts$Entrez_Gene_ID <- NULL
colnames(GBM_counts) <- gsub("-01","",colnames(GBM_counts), fixed=TRUE)
GBM_counts <- t(GBM_counts)
GBM_counts <- as.data.frame(GBM_counts)
GBM_counts$riskscore <- GBM_counts$'7535' * (-0.0582152196382938) + GBM_counts$'494514' * 0.013184642798603 + GBM_counts$'1948' * 0.00470780696956776 + GBM_counts$'147372' * 0.11219886435049 + GBM_counts$'645432' * -0.267800991188706 + GBM_counts$'51806' * -7.50310759764005e-05 + GBM_counts$'56936' * -0.0844511106749959
GBM_counts$bcr_patient_barcode <- rownames(GBM_counts)
GBM_counts <- GBM_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

GBM_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% GBM_counts$bcr_patient_barcode)

GBM_clinicdata <- GBM_clinicdata[match(rownames(GBM_counts),GBM_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

GBM_risk <- left_join(GBM_clinicdata, GBM_counts, by = "bcr_patient_barcode")
GBM_risk$OS <- as.numeric(GBM_risk$OS)
GBM_risk$OS.time <- as.numeric(GBM_risk$OS.time)

stat <- maxstat.test(Surv(GBM_risk$OS.time , GBM_risk$OS) ~ GBM_risk$riskscore, data = GBM_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(GBM_risk$riskscore > cutoff,"high","low"))
GBM_risk$risk <- risk
## Kaplanmeier
GBM_risk$OS.time <- GBM_risk$OS.time / 365 * 12
surv_obj <- Surv(GBM_risk$OS.time , GBM_risk$OS)
fit <- survfit(surv_obj~risk, data = GBM_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "GBM",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = F,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- GBM_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####HNSC####
HNSC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/HNSC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(HNSC_counts) <- HNSC_counts$Entrez_Gene_ID
HNSC_counts <- as.data.frame(HNSC_counts)
HNSC_counts$Entrez_Gene_ID <- NULL
colnames(HNSC_counts) <- gsub("-01","",colnames(HNSC_counts), fixed=TRUE)
HNSC_counts <- t(HNSC_counts)
HNSC_counts <- as.data.frame(HNSC_counts)
HNSC_counts$riskscore <- HNSC_counts$'7535' * (-0.0582152196382938) + HNSC_counts$'494514' * 0.013184642798603 + HNSC_counts$'1948' * 0.00470780696956776 + HNSC_counts$'147372' * 0.11219886435049 + HNSC_counts$'645432' * -0.267800991188706 + HNSC_counts$'51806' * -7.50310759764005e-05 + HNSC_counts$'56936' * -0.0844511106749959
HNSC_counts$bcr_patient_barcode <- rownames(HNSC_counts)
HNSC_counts <- HNSC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

HNSC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% HNSC_counts$bcr_patient_barcode)

HNSC_clinicdata <- HNSC_clinicdata[match(rownames(HNSC_counts),HNSC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

HNSC_risk <- left_join(HNSC_clinicdata, HNSC_counts, by = "bcr_patient_barcode")
HNSC_risk$OS <- as.numeric(HNSC_risk$OS)
HNSC_risk$OS.time <- as.numeric(HNSC_risk$OS.time)

stat <- maxstat.test(Surv(HNSC_risk$OS.time , HNSC_risk$OS) ~ HNSC_risk$riskscore, data = HNSC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(HNSC_risk$riskscore > cutoff,"high","low"))
HNSC_risk$risk <- risk


## mit alten riskeinteilung weil HS
HNSC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% riskall$id)
HNSC_risk_fpkm <- riskall %>%
  dplyr::select(id, risk) %>%
  dplyr::rename(bcr_patient_barcode = id)
HNSC_risk <- left_join(HNSC_clinicdata, HNSC_risk_fpkm, by = "bcr_patient_barcode")
HNSC_risk$OS <- as.numeric(HNSC_risk$OS)
HNSC_risk$OS.time <- as.numeric(HNSC_risk$OS.time)
survdiff(surv_obj~ risk, data = HNSC_risk)
## Kaplanmeier
HNSC_risk$OS.time <- HNSC_risk$OS.time / 365 * 12
surv_obj <- Surv(HNSC_risk$OS.time , HNSC_risk$OS)
fit <- survfit(surv_obj~risk, data = HNSC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "HNSC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           conf.int = T,
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- HNSC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####KICH####
KICH_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/KICH_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(KICH_counts) <- KICH_counts$Entrez_Gene_ID
KICH_counts <- as.data.frame(KICH_counts)
KICH_counts$Entrez_Gene_ID <- NULL
colnames(KICH_counts) <- gsub("-01","",colnames(KICH_counts), fixed=TRUE)
KICH_counts <- t(KICH_counts)
KICH_counts <- as.data.frame(KICH_counts)
KICH_counts$riskscore <- KICH_counts$'7535' * (-0.0582152196382938) + KICH_counts$'494514' * 0.013184642798603 + KICH_counts$'1948' * 0.00470780696956776 + KICH_counts$'147372' * 0.11219886435049 + KICH_counts$'645432' * -0.267800991188706 + KICH_counts$'51806' * -7.50310759764005e-05 + KICH_counts$'56936' * -0.0844511106749959
KICH_counts$bcr_patient_barcode <- rownames(KICH_counts)
KICH_counts <- KICH_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

KICH_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% KICH_counts$bcr_patient_barcode)

KICH_clinicdata <- KICH_clinicdata[match(rownames(KICH_counts),KICH_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

KICH_risk <- left_join(KICH_clinicdata, KICH_counts, by = "bcr_patient_barcode")
KICH_risk$OS <- as.numeric(KICH_risk$OS)
KICH_risk$OS.time <- as.numeric(KICH_risk$OS.time)

stat <- maxstat.test(Surv(KICH_risk$OS.time , KICH_risk$OS) ~ KICH_risk$riskscore, data = KICH_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(KICH_risk$riskscore > cutoff,"high","low"))
KICH_risk$risk <- risk
## Kaplanmeier
KICH_risk$OS.time <- KICH_risk$OS.time / 365 * 12
surv_obj <- Surv(KICH_risk$OS.time , KICH_risk$OS)
fit <- survfit(surv_obj~risk, data = KICH_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "KICH",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- KICH_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####KIRC####
KIRC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/KIRC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(KIRC_counts) <- KIRC_counts$Entrez_Gene_ID
KIRC_counts <- as.data.frame(KIRC_counts)
KIRC_counts$Entrez_Gene_ID <- NULL
colnames(KIRC_counts) <- gsub("-01","",colnames(KIRC_counts), fixed=TRUE)
KIRC_counts <- t(KIRC_counts)
KIRC_counts <- as.data.frame(KIRC_counts)
KIRC_counts$riskscore <- KIRC_counts$'7535' * (-0.0582152196382938) + KIRC_counts$'494514' * 0.013184642798603 + KIRC_counts$'1948' * 0.00470780696956776 + KIRC_counts$'147372' * 0.11219886435049 + KIRC_counts$'645432' * -0.267800991188706 + KIRC_counts$'51806' * -7.50310759764005e-05 + KIRC_counts$'56936' * -0.0844511106749959
KIRC_counts$bcr_patient_barcode <- rownames(KIRC_counts)
KIRC_counts <- KIRC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

KIRC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% KIRC_counts$bcr_patient_barcode)

KIRC_clinicdata <- KIRC_clinicdata[match(rownames(KIRC_counts),KIRC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

KIRC_risk <- left_join(KIRC_clinicdata, KIRC_counts, by = "bcr_patient_barcode")
KIRC_risk$OS <- as.numeric(KIRC_risk$OS)
KIRC_risk$OS.time <- as.numeric(KIRC_risk$OS.time)

stat <- maxstat.test(Surv(KIRC_risk$OS.time , KIRC_risk$OS) ~ KIRC_risk$riskscore, data = KIRC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(KIRC_risk$riskscore > cutoff,"high","low"))
KIRC_risk$risk <- risk
## Kaplanmeier
KIRC_risk$OS.time <- KIRC_risk$OS.time / 365 * 12
surv_obj <- Surv(KIRC_risk$OS.time , KIRC_risk$OS)
fit <- survfit(surv_obj~risk, data = KIRC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "KIRC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           conf.int = T,
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- KIRC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####KIRP####
KIRP_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/KIRP_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(KIRP_counts) <- KIRP_counts$Entrez_Gene_ID
KIRP_counts <- as.data.frame(KIRP_counts)
KIRP_counts$Entrez_Gene_ID <- NULL
colnames(KIRP_counts) <- gsub("-01","",colnames(KIRP_counts), fixed=TRUE)
KIRP_counts <- t(KIRP_counts)
KIRP_counts <- as.data.frame(KIRP_counts)
KIRP_counts$riskscore <- KIRP_counts$'7535' * (-0.0582152196382938) + KIRP_counts$'494514' * 0.013184642798603 + KIRP_counts$'1948' * 0.00470780696956776 + KIRP_counts$'147372' * 0.11219886435049 + KIRP_counts$'645432' * -0.267800991188706 + KIRP_counts$'51806' * -7.50310759764005e-05 + KIRP_counts$'56936' * -0.0844511106749959
KIRP_counts$bcr_patient_barcode <- rownames(KIRP_counts)
KIRP_counts <- KIRP_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

KIRP_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% KIRP_counts$bcr_patient_barcode)

KIRP_clinicdata <- KIRP_clinicdata[match(rownames(KIRP_counts),KIRP_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

KIRP_risk <- left_join(KIRP_clinicdata, KIRP_counts, by = "bcr_patient_barcode")
KIRP_risk$OS <- as.numeric(KIRP_risk$OS)
KIRP_risk$OS.time <- as.numeric(KIRP_risk$OS.time)

stat <- maxstat.test(Surv(KIRP_risk$OS.time , KIRP_risk$OS) ~ KIRP_risk$riskscore, data = KIRP_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(KIRP_risk$riskscore > cutoff,"high","low"))
KIRP_risk$risk <- risk
## Kaplanmeier
KIRP_risk$OS.time <- KIRP_risk$OS.time / 365 * 12
surv_obj <- Surv(KIRP_risk$OS.time , KIRP_risk$OS)
fit <- survfit(surv_obj~risk, data = KIRP_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "KIRP",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- KIRP_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####LGG####
LGG_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/LGG_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(LGG_counts) <- LGG_counts$Entrez_Gene_ID
LGG_counts <- as.data.frame(LGG_counts)
LGG_counts$Entrez_Gene_ID <- NULL
colnames(LGG_counts) <- gsub("-01","",colnames(LGG_counts), fixed=TRUE)
LGG_counts <- t(LGG_counts)
LGG_counts <- as.data.frame(LGG_counts)
LGG_counts$riskscore <- LGG_counts$'7535' * (-0.0582152196382938) + LGG_counts$'494514' * 0.013184642798603 + LGG_counts$'1948' * 0.00470780696956776 + LGG_counts$'147372' * 0.11219886435049 + LGG_counts$'645432' * -0.267800991188706 + LGG_counts$'51806' * -7.50310759764005e-05 + LGG_counts$'56936' * -0.0844511106749959
LGG_counts$bcr_patient_barcode <- rownames(LGG_counts)
LGG_counts <- LGG_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

LGG_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% LGG_counts$bcr_patient_barcode)

LGG_clinicdata <- LGG_clinicdata[match(rownames(LGG_counts),LGG_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

LGG_risk <- left_join(LGG_clinicdata, LGG_counts, by = "bcr_patient_barcode")
LGG_risk$OS <- as.numeric(LGG_risk$OS)
LGG_risk$OS.time <- as.numeric(LGG_risk$OS.time)

stat <- maxstat.test(Surv(LGG_risk$OS.time , LGG_risk$OS) ~ LGG_risk$riskscore, data = LGG_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(LGG_risk$riskscore > cutoff,"high","low"))
LGG_risk$risk <- risk
## Kaplanmeier
LGG_risk$OS.time <- LGG_risk$OS.time / 365 * 12
surv_obj <- Surv(LGG_risk$OS.time , LGG_risk$OS)
fit <- survfit(surv_obj~risk, data = LGG_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "LGG",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- LGG_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####LIHC####
LIHC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/LIHC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(LIHC_counts) <- LIHC_counts$Entrez_Gene_ID
LIHC_counts <- as.data.frame(LIHC_counts)
LIHC_counts$Entrez_Gene_ID <- NULL
colnames(LIHC_counts) <- gsub("-01","",colnames(LIHC_counts), fixed=TRUE)
LIHC_counts <- t(LIHC_counts)
LIHC_counts <- as.data.frame(LIHC_counts)
LIHC_counts$riskscore <- LIHC_counts$'7535' * (-0.0582152196382938) + LIHC_counts$'494514' * 0.013184642798603 + LIHC_counts$'1948' * 0.00470780696956776 + LIHC_counts$'147372' * 0.11219886435049 + LIHC_counts$'645432' * -0.267800991188706 + LIHC_counts$'51806' * -7.50310759764005e-05 + LIHC_counts$'56936' * -0.0844511106749959
LIHC_counts$bcr_patient_barcode <- rownames(LIHC_counts)
LIHC_counts <- LIHC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

LIHC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% LIHC_counts$bcr_patient_barcode)

LIHC_clinicdata <- LIHC_clinicdata[match(rownames(LIHC_counts),LIHC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

LIHC_risk <- left_join(LIHC_clinicdata, LIHC_counts, by = "bcr_patient_barcode")
LIHC_risk$OS <- as.numeric(LIHC_risk$OS)
LIHC_risk$OS.time <- as.numeric(LIHC_risk$OS.time)

stat <- maxstat.test(Surv(LIHC_risk$OS.time , LIHC_risk$OS) ~ LIHC_risk$riskscore, data = LIHC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(LIHC_risk$riskscore > cutoff,"high","low"))
LIHC_risk$risk <- risk
## Kaplanmeier
LIHC_risk$OS.time <- LIHC_risk$OS.time / 365 * 12
surv_obj <- Surv(LIHC_risk$OS.time , LIHC_risk$OS)
fit <- survfit(surv_obj~risk, data = LIHC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "LIHC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- LIHC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####LUAD####
LUAD_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/LUAD_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(LUAD_counts) <- LUAD_counts$Entrez_Gene_ID
LUAD_counts <- as.data.frame(LUAD_counts)
LUAD_counts$Entrez_Gene_ID <- NULL
colnames(LUAD_counts) <- gsub("-01","",colnames(LUAD_counts), fixed=TRUE)
LUAD_counts <- t(LUAD_counts)
LUAD_counts <- as.data.frame(LUAD_counts)
LUAD_counts$riskscore <- LUAD_counts$'7535' * (-0.0582152196382938) + LUAD_counts$'494514' * 0.013184642798603 + LUAD_counts$'1948' * 0.00470780696956776 + LUAD_counts$'147372' * 0.11219886435049 + LUAD_counts$'645432' * -0.267800991188706 + LUAD_counts$'51806' * -7.50310759764005e-05 + LUAD_counts$'56936' * -0.0844511106749959
LUAD_counts$bcr_patient_barcode <- rownames(LUAD_counts)
LUAD_counts <- LUAD_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

LUAD_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% LUAD_counts$bcr_patient_barcode)

LUAD_clinicdata <- LUAD_clinicdata[match(rownames(LUAD_counts),LUAD_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

LUAD_risk <- left_join(LUAD_clinicdata, LUAD_counts, by = "bcr_patient_barcode")
LUAD_risk$OS <- as.numeric(LUAD_risk$OS)
LUAD_risk$OS.time <- as.numeric(LUAD_risk$OS.time)

stat <- maxstat.test(Surv(LUAD_risk$OS.time , LUAD_risk$OS) ~ LUAD_risk$riskscore, data = LUAD_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(LUAD_risk$riskscore > cutoff,"high","low"))
LUAD_risk$risk <- risk
## Kaplanmeier
LUAD_risk$OS.time <- LUAD_risk$OS.time / 365 * 12
surv_obj <- Surv(LUAD_risk$OS.time , LUAD_risk$OS)
fit <- survfit(surv_obj~risk, data = LUAD_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "LUAD",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- LUAD_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####LUSC####
LUSC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/LUSC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(LUSC_counts) <- LUSC_counts$Entrez_Gene_ID
LUSC_counts <- as.data.frame(LUSC_counts)
LUSC_counts$Entrez_Gene_ID <- NULL
colnames(LUSC_counts) <- gsub("-01","",colnames(LUSC_counts), fixed=TRUE)
LUSC_counts <- t(LUSC_counts)
LUSC_counts <- as.data.frame(LUSC_counts)
LUSC_counts$riskscore <- LUSC_counts$'7535' * (-0.0582152196382938) + LUSC_counts$'494514' * 0.013184642798603 + LUSC_counts$'1948' * 0.00470780696956776 + LUSC_counts$'147372' * 0.11219886435049 + LUSC_counts$'645432' * -0.267800991188706 + LUSC_counts$'51806' * -7.50310759764005e-05 + LUSC_counts$'56936' * -0.0844511106749959
LUSC_counts$bcr_patient_barcode <- rownames(LUSC_counts)
LUSC_counts <- LUSC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

LUSC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% LUSC_counts$bcr_patient_barcode)

LUSC_clinicdata <- LUSC_clinicdata[match(rownames(LUSC_counts),LUSC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

LUSC_risk <- left_join(LUSC_clinicdata, LUSC_counts, by = "bcr_patient_barcode")
LUSC_risk$OS <- as.numeric(LUSC_risk$OS)
LUSC_risk$OS.time <- as.numeric(LUSC_risk$OS.time)

stat <- maxstat.test(Surv(LUSC_risk$OS.time , LUSC_risk$OS) ~ LUSC_risk$riskscore, data = LUSC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(LUSC_risk$riskscore > cutoff,"high","low"))
LUSC_risk$risk <- risk
## Kaplanmeier
LUSC_risk$OS.time <- LUSC_risk$OS.time / 365 * 12
surv_obj <- Surv(LUSC_risk$OS.time , LUSC_risk$OS)
fit <- survfit(surv_obj~risk, data = LUSC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "LUSC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- LUSC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####MESO####
MESO_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/MESO_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(MESO_counts) <- MESO_counts$Entrez_Gene_ID
MESO_counts <- as.data.frame(MESO_counts)
MESO_counts$Entrez_Gene_ID <- NULL
colnames(MESO_counts) <- gsub("-01","",colnames(MESO_counts), fixed=TRUE)
MESO_counts <- t(MESO_counts)
MESO_counts <- as.data.frame(MESO_counts)
MESO_counts$riskscore <- MESO_counts$'7535' * (-0.0582152196382938) + MESO_counts$'494514' * 0.013184642798603 + MESO_counts$'1948' * 0.00470780696956776 + MESO_counts$'147372' * 0.11219886435049 + MESO_counts$'645432' * -0.267800991188706 + MESO_counts$'51806' * -7.50310759764005e-05 + MESO_counts$'56936' * -0.0844511106749959
MESO_counts$bcr_patient_barcode <- rownames(MESO_counts)
MESO_counts <- MESO_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

MESO_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% MESO_counts$bcr_patient_barcode)

MESO_clinicdata <- MESO_clinicdata[match(rownames(MESO_counts),MESO_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

MESO_risk <- left_join(MESO_clinicdata, MESO_counts, by = "bcr_patient_barcode")
MESO_risk$OS <- as.numeric(MESO_risk$OS)
MESO_risk$OS.time <- as.numeric(MESO_risk$OS.time)

stat <- maxstat.test(Surv(MESO_risk$OS.time , MESO_risk$OS) ~ MESO_risk$riskscore, data = MESO_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(MESO_risk$riskscore > cutoff,"high","low"))
MESO_risk$risk <- risk
## Kaplanmeier
MESO_risk$OS.time <- MESO_risk$OS.time / 365 * 12
surv_obj <- Surv(MESO_risk$OS.time , MESO_risk$OS)
fit <- survfit(surv_obj~risk, data = MESO_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "MESO",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- MESO_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####OV####
OV_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/OV_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(OV_counts) <- OV_counts$Entrez_Gene_ID
OV_counts <- as.data.frame(OV_counts)
OV_counts$Entrez_Gene_ID <- NULL
colnames(OV_counts) <- gsub("-01","",colnames(OV_counts), fixed=TRUE)
OV_counts <- t(OV_counts)
OV_counts <- as.data.frame(OV_counts)
OV_counts$riskscore <- OV_counts$'7535' * (-0.0582152196382938) + OV_counts$'494514' * 0.013184642798603 + OV_counts$'1948' * 0.00470780696956776 + OV_counts$'147372' * 0.11219886435049 + OV_counts$'645432' * -0.267800991188706 + OV_counts$'51806' * -7.50310759764005e-05 + OV_counts$'56936' * -0.0844511106749959
OV_counts$bcr_patient_barcode <- rownames(OV_counts)
OV_counts <- OV_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

OV_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% OV_counts$bcr_patient_barcode)

OV_clinicdata <- OV_clinicdata[match(rownames(OV_counts),OV_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

OV_risk <- left_join(OV_clinicdata, OV_counts, by = "bcr_patient_barcode")
OV_risk$OS <- as.numeric(OV_risk$OS)
OV_risk$OS.time <- as.numeric(OV_risk$OS.time)

stat <- maxstat.test(Surv(OV_risk$OS.time , OV_risk$OS) ~ OV_risk$riskscore, data = OV_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(OV_risk$riskscore > cutoff,"high","low"))
OV_risk$risk <- risk
## Kaplanmeier
OV_risk$OS.time <- OV_risk$OS.time / 365 * 12
surv_obj <- Surv(OV_risk$OS.time , OV_risk$OS)
fit <- survfit(surv_obj~risk, data = OV_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "OV",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- OV_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####PAAD####
PAAD_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/PAAD_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(PAAD_counts) <- PAAD_counts$Entrez_Gene_ID
PAAD_counts <- as.data.frame(PAAD_counts)
PAAD_counts$Entrez_Gene_ID <- NULL
colnames(PAAD_counts) <- gsub("-01","",colnames(PAAD_counts), fixed=TRUE)
PAAD_counts <- t(PAAD_counts)
PAAD_counts <- as.data.frame(PAAD_counts)
PAAD_counts$riskscore <- PAAD_counts$'7535' * (-0.0582152196382938) + PAAD_counts$'494514' * 0.013184642798603 + PAAD_counts$'1948' * 0.00470780696956776 + PAAD_counts$'147372' * 0.11219886435049 + PAAD_counts$'645432' * -0.267800991188706 + PAAD_counts$'51806' * -7.50310759764005e-05 + PAAD_counts$'56936' * -0.0844511106749959
PAAD_counts$bcr_patient_barcode <- rownames(PAAD_counts)
PAAD_counts <- PAAD_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

PAAD_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% PAAD_counts$bcr_patient_barcode)

PAAD_clinicdata <- PAAD_clinicdata[match(rownames(PAAD_counts),PAAD_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

PAAD_risk <- left_join(PAAD_clinicdata, PAAD_counts, by = "bcr_patient_barcode")
PAAD_risk$OS <- as.numeric(PAAD_risk$OS)
PAAD_risk$OS.time <- as.numeric(PAAD_risk$OS.time)

stat <- maxstat.test(Surv(PAAD_risk$OS.time , PAAD_risk$OS) ~ PAAD_risk$riskscore, data = PAAD_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(PAAD_risk$riskscore > cutoff,"high","low"))
PAAD_risk$risk <- risk
## Kaplanmeier
PAAD_risk$OS.time <- PAAD_risk$OS.time / 365 * 12
surv_obj <- Surv(PAAD_risk$OS.time , PAAD_risk$OS)
fit <- survfit(surv_obj~risk, data = PAAD_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "PAAD",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- PAAD_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####PCPG####
PCPG_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/PCPG_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(PCPG_counts) <- PCPG_counts$Entrez_Gene_ID
PCPG_counts <- as.data.frame(PCPG_counts)
PCPG_counts$Entrez_Gene_ID <- NULL
colnames(PCPG_counts) <- gsub("-01","",colnames(PCPG_counts), fixed=TRUE)
PCPG_counts <- t(PCPG_counts)
PCPG_counts <- as.data.frame(PCPG_counts)
PCPG_counts$riskscore <- PCPG_counts$'7535' * (-0.0582152196382938) + PCPG_counts$'494514' * 0.013184642798603 + PCPG_counts$'1948' * 0.00470780696956776 + PCPG_counts$'147372' * 0.11219886435049 + PCPG_counts$'645432' * -0.267800991188706 + PCPG_counts$'51806' * -7.50310759764005e-05 + PCPG_counts$'56936' * -0.0844511106749959
PCPG_counts$bcr_patient_barcode <- rownames(PCPG_counts)
PCPG_counts <- PCPG_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

PCPG_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% PCPG_counts$bcr_patient_barcode)

PCPG_clinicdata <- PCPG_clinicdata[match(rownames(PCPG_counts),PCPG_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

PCPG_risk <- left_join(PCPG_clinicdata, PCPG_counts, by = "bcr_patient_barcode")
PCPG_risk$OS <- as.numeric(PCPG_risk$OS)
PCPG_risk$OS.time <- as.numeric(PCPG_risk$OS.time)

stat <- maxstat.test(Surv(PCPG_risk$OS.time , PCPG_risk$OS) ~ PCPG_risk$riskscore, data = PCPG_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(PCPG_risk$riskscore > cutoff,"high","low"))
PCPG_risk$risk <- risk
## Kaplanmeier
PCPG_risk$OS.time <- PCPG_risk$OS.time / 365 * 12
surv_obj <- Surv(PCPG_risk$OS.time , PCPG_risk$OS)
fit <- survfit(surv_obj~risk, data = PCPG_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "PCPG",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- PCPG_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####PRAD####
PRAD_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/PRAD_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(PRAD_counts) <- PRAD_counts$Entrez_Gene_ID
PRAD_counts <- as.data.frame(PRAD_counts)
PRAD_counts$Entrez_Gene_ID <- NULL
colnames(PRAD_counts) <- gsub("-01","",colnames(PRAD_counts), fixed=TRUE)
PRAD_counts <- t(PRAD_counts)
PRAD_counts <- as.data.frame(PRAD_counts)
PRAD_counts$riskscore <- PRAD_counts$'7535' * (-0.0582152196382938) + PRAD_counts$'494514' * 0.013184642798603 + PRAD_counts$'1948' * 0.00470780696956776 + PRAD_counts$'147372' * 0.11219886435049 + PRAD_counts$'645432' * -0.267800991188706 + PRAD_counts$'51806' * -7.50310759764005e-05 + PRAD_counts$'56936' * -0.0844511106749959
PRAD_counts$bcr_patient_barcode <- rownames(PRAD_counts)
PRAD_counts <- PRAD_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

PRAD_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% PRAD_counts$bcr_patient_barcode)

PRAD_clinicdata <- PRAD_clinicdata[match(rownames(PRAD_counts),PRAD_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

PRAD_risk <- left_join(PRAD_clinicdata, PRAD_counts, by = "bcr_patient_barcode")
PRAD_risk$OS <- as.numeric(PRAD_risk$OS)
PRAD_risk$OS.time <- as.numeric(PRAD_risk$OS.time)

stat <- maxstat.test(Surv(PRAD_risk$OS.time , PRAD_risk$OS) ~ PRAD_risk$riskscore, data = PRAD_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(PRAD_risk$riskscore > cutoff,"high","low"))
PRAD_risk$risk <- risk
## Kaplanmeier
PRAD_risk$OS.time <- PRAD_risk$OS.time / 365 * 12
surv_obj <- Surv(PRAD_risk$OS.time , PRAD_risk$OS)
fit <- survfit(surv_obj~risk, data = PRAD_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "PRAD",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- PRAD_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####READ####
READ_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/READ_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(READ_counts) <- READ_counts$Entrez_Gene_ID
READ_counts <- as.data.frame(READ_counts)
READ_counts$Entrez_Gene_ID <- NULL
colnames(READ_counts) <- gsub("-01","",colnames(READ_counts), fixed=TRUE)
READ_counts <- t(READ_counts)
READ_counts <- as.data.frame(READ_counts)
READ_counts$riskscore <- READ_counts$'7535' * (-0.0582152196382938) + READ_counts$'494514' * 0.013184642798603 + READ_counts$'1948' * 0.00470780696956776 + READ_counts$'147372' * 0.11219886435049 + READ_counts$'645432' * -0.267800991188706 + READ_counts$'51806' * -7.50310759764005e-05 + READ_counts$'56936' * -0.0844511106749959
READ_counts$bcr_patient_barcode <- rownames(READ_counts)
READ_counts <- READ_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

READ_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% READ_counts$bcr_patient_barcode)

READ_clinicdata <- READ_clinicdata[match(rownames(READ_counts),READ_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

READ_risk <- left_join(READ_clinicdata, READ_counts, by = "bcr_patient_barcode")
READ_risk$OS <- as.numeric(READ_risk$OS)
READ_risk$OS.time <- as.numeric(READ_risk$OS.time)

stat <- maxstat.test(Surv(READ_risk$OS.time , READ_risk$OS) ~ READ_risk$riskscore, data = READ_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(READ_risk$riskscore > cutoff,"high","low"))
READ_risk$risk <- risk
## Kaplanmeier
READ_risk$OS.time <- READ_risk$OS.time / 365 * 12
surv_obj <- Surv(READ_risk$OS.time , READ_risk$OS)
fit <- survfit(surv_obj~risk, data = READ_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "READ",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- READ_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####SARC####
SARC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/SARC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(SARC_counts) <- SARC_counts$Entrez_Gene_ID
SARC_counts <- as.data.frame(SARC_counts)
SARC_counts$Entrez_Gene_ID <- NULL
colnames(SARC_counts) <- gsub("-01","",colnames(SARC_counts), fixed=TRUE)
SARC_counts <- t(SARC_counts)
SARC_counts <- as.data.frame(SARC_counts)
SARC_counts$riskscore <- SARC_counts$'7535' * (-0.0582152196382938) + SARC_counts$'494514' * 0.013184642798603 + SARC_counts$'1948' * 0.00470780696956776 + SARC_counts$'147372' * 0.11219886435049 + SARC_counts$'645432' * -0.267800991188706 + SARC_counts$'51806' * -7.50310759764005e-05 + SARC_counts$'56936' * -0.0844511106749959
SARC_counts$bcr_patient_barcode <- rownames(SARC_counts)
SARC_counts <- SARC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

SARC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% SARC_counts$bcr_patient_barcode)

SARC_clinicdata <- SARC_clinicdata[match(rownames(SARC_counts),SARC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

SARC_risk <- left_join(SARC_clinicdata, SARC_counts, by = "bcr_patient_barcode")
SARC_risk$OS <- as.numeric(SARC_risk$OS)
SARC_risk$OS.time <- as.numeric(SARC_risk$OS.time)

stat <- maxstat.test(Surv(SARC_risk$OS.time , SARC_risk$OS) ~ SARC_risk$riskscore, data = SARC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(SARC_risk$riskscore > cutoff,"high","low"))
SARC_risk$risk <- risk
## Kaplanmeier
SARC_risk$OS.time <- SARC_risk$OS.time / 365 * 12
surv_obj <- Surv(SARC_risk$OS.time , SARC_risk$OS)
fit <- survfit(surv_obj~risk, data = SARC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "SARC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- SARC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####SKCM####
SKCM_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/SKCM_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(SKCM_counts) <- SKCM_counts$Entrez_Gene_ID
SKCM_counts <- as.data.frame(SKCM_counts)
SKCM_counts$Entrez_Gene_ID <- NULL
colnames(SKCM_counts) <- gsub("-01","",colnames(SKCM_counts), fixed=TRUE)
SKCM_counts <- t(SKCM_counts)
SKCM_counts <- as.data.frame(SKCM_counts)
SKCM_counts$riskscore <- SKCM_counts$'7535' * (-0.0582152196382938) + SKCM_counts$'494514' * 0.013184642798603 + SKCM_counts$'1948' * 0.00470780696956776 + SKCM_counts$'147372' * 0.11219886435049 + SKCM_counts$'645432' * -0.267800991188706 + SKCM_counts$'51806' * -7.50310759764005e-05 + SKCM_counts$'56936' * -0.0844511106749959
SKCM_counts$bcr_patient_barcode <- rownames(SKCM_counts)
SKCM_counts <- SKCM_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

SKCM_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% SKCM_counts$bcr_patient_barcode)

SKCM_clinicdata <- SKCM_clinicdata[match(rownames(SKCM_counts),SKCM_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

SKCM_risk <- left_join(SKCM_clinicdata, SKCM_counts, by = "bcr_patient_barcode")
SKCM_risk$OS <- as.numeric(SKCM_risk$OS)
SKCM_risk$OS.time <- as.numeric(SKCM_risk$OS.time)

stat <- maxstat.test(Surv(SKCM_risk$OS.time , SKCM_risk$OS) ~ SKCM_risk$riskscore, data = SKCM_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(SKCM_risk$riskscore > cutoff,"high","low"))
SKCM_risk$risk <- risk
## Kaplanmeier
SKCM_risk$OS.time <- SKCM_risk$OS.time / 365 * 12
surv_obj <- Surv(SKCM_risk$OS.time , SKCM_risk$OS)
fit <- survfit(surv_obj~risk, data = SKCM_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "SKCM",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- SKCM_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####STAD####
STAD_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/STAD_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(STAD_counts) <- STAD_counts$Entrez_Gene_ID
STAD_counts <- as.data.frame(STAD_counts)
STAD_counts$Entrez_Gene_ID <- NULL
colnames(STAD_counts) <- gsub("-01","",colnames(STAD_counts), fixed=TRUE)
STAD_counts <- t(STAD_counts)
STAD_counts <- as.data.frame(STAD_counts)
STAD_counts$riskscore <- STAD_counts$'7535' * (-0.0582152196382938) + STAD_counts$'494514' * 0.013184642798603 + STAD_counts$'1948' * 0.00470780696956776 + STAD_counts$'147372' * 0.11219886435049 + STAD_counts$'645432' * -0.267800991188706 + STAD_counts$'51806' * -7.50310759764005e-05 + STAD_counts$'56936' * -0.0844511106749959
STAD_counts$bcr_patient_barcode <- rownames(STAD_counts)
STAD_counts <- STAD_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

STAD_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% STAD_counts$bcr_patient_barcode)

STAD_clinicdata <- STAD_clinicdata[match(rownames(STAD_counts),STAD_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

STAD_risk <- left_join(STAD_clinicdata, STAD_counts, by = "bcr_patient_barcode")
STAD_risk$OS <- as.numeric(STAD_risk$OS)
STAD_risk$OS.time <- as.numeric(STAD_risk$OS.time)

stat <- maxstat.test(Surv(STAD_risk$OS.time , STAD_risk$OS) ~ STAD_risk$riskscore, data = STAD_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(STAD_risk$riskscore > cutoff,"high","low"))
STAD_risk$risk <- risk
## Kaplanmeier
STAD_risk$OS.time <- STAD_risk$OS.time / 365 * 12
surv_obj <- Surv(STAD_risk$OS.time , STAD_risk$OS)
fit <- survfit(surv_obj~risk, data = STAD_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "STAD",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- STAD_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####TGCT####
TGCT_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/TGCT_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(TGCT_counts) <- TGCT_counts$Entrez_Gene_ID
TGCT_counts <- as.data.frame(TGCT_counts)
TGCT_counts$Entrez_Gene_ID <- NULL
colnames(TGCT_counts) <- gsub("-01","",colnames(TGCT_counts), fixed=TRUE)
TGCT_counts <- t(TGCT_counts)
TGCT_counts <- as.data.frame(TGCT_counts)
TGCT_counts$riskscore <- TGCT_counts$'7535' * (-0.0582152196382938) + TGCT_counts$'494514' * 0.013184642798603 + TGCT_counts$'1948' * 0.00470780696956776 + TGCT_counts$'147372' * 0.11219886435049 + TGCT_counts$'645432' * -0.267800991188706 + TGCT_counts$'51806' * -7.50310759764005e-05 + TGCT_counts$'56936' * -0.0844511106749959
TGCT_counts$bcr_patient_barcode <- rownames(TGCT_counts)
TGCT_counts <- TGCT_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

TGCT_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% TGCT_counts$bcr_patient_barcode)

TGCT_clinicdata <- TGCT_clinicdata[match(rownames(TGCT_counts),TGCT_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

TGCT_risk <- left_join(TGCT_clinicdata, TGCT_counts, by = "bcr_patient_barcode")
TGCT_risk$OS <- as.numeric(TGCT_risk$OS)
TGCT_risk$OS.time <- as.numeric(TGCT_risk$OS.time)

stat <- maxstat.test(Surv(TGCT_risk$OS.time , TGCT_risk$OS) ~ TGCT_risk$riskscore, data = TGCT_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(TGCT_risk$riskscore > cutoff,"high","low"))
TGCT_risk$risk <- risk
## Kaplanmeier
TGCT_risk$OS.time <- TGCT_risk$OS.time / 365 * 12
surv_obj <- Surv(TGCT_risk$OS.time , TGCT_risk$OS)
fit <- survfit(surv_obj~risk, data = TGCT_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "TGCT",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- TGCT_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####THCA####
THCA_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/THCA_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(THCA_counts) <- THCA_counts$Entrez_Gene_ID
THCA_counts <- as.data.frame(THCA_counts)
THCA_counts$Entrez_Gene_ID <- NULL
colnames(THCA_counts) <- gsub("-01","",colnames(THCA_counts), fixed=TRUE)
THCA_counts <- t(THCA_counts)
THCA_counts <- as.data.frame(THCA_counts)
THCA_counts$riskscore <- THCA_counts$'7535' * (-0.0582152196382938) + THCA_counts$'494514' * 0.013184642798603 + THCA_counts$'1948' * 0.00470780696956776 + THCA_counts$'147372' * 0.11219886435049 + THCA_counts$'645432' * -0.267800991188706 + THCA_counts$'51806' * -7.50310759764005e-05 + THCA_counts$'56936' * -0.0844511106749959
THCA_counts$bcr_patient_barcode <- rownames(THCA_counts)
THCA_counts <- THCA_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

THCA_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% THCA_counts$bcr_patient_barcode)

THCA_clinicdata <- THCA_clinicdata[match(rownames(THCA_counts),THCA_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

THCA_risk <- left_join(THCA_clinicdata, THCA_counts, by = "bcr_patient_barcode")
THCA_risk$OS <- as.numeric(THCA_risk$OS)
THCA_risk$OS.time <- as.numeric(THCA_risk$OS.time)

stat <- maxstat.test(Surv(THCA_risk$OS.time , THCA_risk$OS) ~ THCA_risk$riskscore, data = THCA_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(THCA_risk$riskscore > cutoff,"high","low"))
THCA_risk$risk <- risk
## Kaplanmeier
THCA_risk$OS.time <- THCA_risk$OS.time / 365 * 12
surv_obj <- Surv(THCA_risk$OS.time , THCA_risk$OS)
fit <- survfit(surv_obj~risk, data = THCA_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "THCA",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- THCA_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####THYM####
THYM_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/THYM_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(THYM_counts) <- THYM_counts$Entrez_Gene_ID
THYM_counts <- as.data.frame(THYM_counts)
THYM_counts$Entrez_Gene_ID <- NULL
colnames(THYM_counts) <- gsub("-01","",colnames(THYM_counts), fixed=TRUE)
THYM_counts <- t(THYM_counts)
THYM_counts <- as.data.frame(THYM_counts)
THYM_counts$riskscore <- THYM_counts$'7535' * (-0.0582152196382938) + THYM_counts$'494514' * 0.013184642798603 + THYM_counts$'1948' * 0.00470780696956776 + THYM_counts$'147372' * 0.11219886435049 + THYM_counts$'645432' * -0.267800991188706 + THYM_counts$'51806' * -7.50310759764005e-05 + THYM_counts$'56936' * -0.0844511106749959
THYM_counts$bcr_patient_barcode <- rownames(THYM_counts)
THYM_counts <- THYM_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

THYM_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% THYM_counts$bcr_patient_barcode)

THYM_clinicdata <- THYM_clinicdata[match(rownames(THYM_counts),THYM_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

THYM_risk <- left_join(THYM_clinicdata, THYM_counts, by = "bcr_patient_barcode")
THYM_risk$OS <- as.numeric(THYM_risk$OS)
THYM_risk$OS.time <- as.numeric(THYM_risk$OS.time)

stat <- maxstat.test(Surv(THYM_risk$OS.time , THYM_risk$OS) ~ THYM_risk$riskscore, data = THYM_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(THYM_risk$riskscore > cutoff,"high","low"))
THYM_risk$risk <- risk
## Kaplanmeier
THYM_risk$OS.time <- THYM_risk$OS.time / 365 * 12
surv_obj <- Surv(THYM_risk$OS.time , THYM_risk$OS)
fit <- survfit(surv_obj~risk, data = THYM_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "THYM",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- THYM_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####UCEC####
UCEC_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/UCEC_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(UCEC_counts) <- UCEC_counts$Entrez_Gene_ID
UCEC_counts <- as.data.frame(UCEC_counts)
UCEC_counts$Entrez_Gene_ID <- NULL
colnames(UCEC_counts) <- gsub("-01","",colnames(UCEC_counts), fixed=TRUE)
UCEC_counts <- t(UCEC_counts)
UCEC_counts <- as.data.frame(UCEC_counts)
UCEC_counts$riskscore <- UCEC_counts$'7535' * (-0.0582152196382938) + UCEC_counts$'494514' * 0.013184642798603 + UCEC_counts$'1948' * 0.00470780696956776 + UCEC_counts$'147372' * 0.11219886435049 + UCEC_counts$'645432' * -0.267800991188706 + UCEC_counts$'51806' * -7.50310759764005e-05 + UCEC_counts$'56936' * -0.0844511106749959
UCEC_counts$bcr_patient_barcode <- rownames(UCEC_counts)
UCEC_counts <- UCEC_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

UCEC_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% UCEC_counts$bcr_patient_barcode)

UCEC_clinicdata <- UCEC_clinicdata[match(rownames(UCEC_counts),UCEC_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

UCEC_risk <- left_join(UCEC_clinicdata, UCEC_counts, by = "bcr_patient_barcode")
UCEC_risk$OS <- as.numeric(UCEC_risk$OS)
UCEC_risk$OS.time <- as.numeric(UCEC_risk$OS.time)

stat <- maxstat.test(Surv(UCEC_risk$OS.time , UCEC_risk$OS) ~ UCEC_risk$riskscore, data = UCEC_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(UCEC_risk$riskscore > cutoff,"high","low"))
UCEC_risk$risk <- risk
## Kaplanmeier
UCEC_risk$OS.time <- UCEC_risk$OS.time / 365 * 12
surv_obj <- Surv(UCEC_risk$OS.time , UCEC_risk$OS)
fit <- survfit(surv_obj~risk, data = UCEC_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "UCEC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- UCEC_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####UCS####
UCS_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/UCS_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(UCS_counts) <- UCS_counts$Entrez_Gene_ID
UCS_counts <- as.data.frame(UCS_counts)
UCS_counts$Entrez_Gene_ID <- NULL
colnames(UCS_counts) <- gsub("-01","",colnames(UCS_counts), fixed=TRUE)
UCS_counts <- t(UCS_counts)
UCS_counts <- as.data.frame(UCS_counts)
UCS_counts$riskscore <- UCS_counts$'7535' * (-0.0582152196382938) + UCS_counts$'494514' * 0.013184642798603 + UCS_counts$'1948' * 0.00470780696956776 + UCS_counts$'147372' * 0.11219886435049 + UCS_counts$'645432' * -0.267800991188706 + UCS_counts$'51806' * -7.50310759764005e-05 + UCS_counts$'56936' * -0.0844511106749959
UCS_counts$bcr_patient_barcode <- rownames(UCS_counts)
UCS_counts <- UCS_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

UCS_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% UCS_counts$bcr_patient_barcode)

UCS_clinicdata <- UCS_clinicdata[match(rownames(UCS_counts),UCS_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

UCS_risk <- left_join(UCS_clinicdata, UCS_counts, by = "bcr_patient_barcode")
UCS_risk$OS <- as.numeric(UCS_risk$OS)
UCS_risk$OS.time <- as.numeric(UCS_risk$OS.time)

stat <- maxstat.test(Surv(UCS_risk$OS.time , UCS_risk$OS) ~ UCS_risk$riskscore, data = UCS_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(UCS_risk$riskscore > cutoff,"high","low"))
UCS_risk$risk <- risk
## Kaplanmeier
UCS_risk$OS.time <- UCS_risk$OS.time / 365 * 12
surv_obj <- Surv(UCS_risk$OS.time , UCS_risk$OS)
fit <- survfit(surv_obj~risk, data = UCS_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "UCS",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- UCS_risk 
risk_cancerproof <- risk_cancerproof%>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)
#####UVM####
UVM_counts <- read_delim("~/Desktop/Doktorarbeit/Data/Counts/UVM_counts_TCGA.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(UVM_counts) <- UVM_counts$Entrez_Gene_ID
UVM_counts <- as.data.frame(UVM_counts)
UVM_counts$Entrez_Gene_ID <- NULL
colnames(UVM_counts) <- gsub("-01","",colnames(UVM_counts), fixed=TRUE)
UVM_counts <- t(UVM_counts)
UVM_counts <- as.data.frame(UVM_counts)
UVM_counts$riskscore <- UVM_counts$'7535' * (-0.0582152196382938) + UVM_counts$'494514' * 0.013184642798603 + UVM_counts$'1948' * 0.00470780696956776 + UVM_counts$'147372' * 0.11219886435049 + UVM_counts$'645432' * -0.267800991188706 + UVM_counts$'51806' * -7.50310759764005e-05 + UVM_counts$'56936' * -0.0844511106749959
UVM_counts$bcr_patient_barcode <- rownames(UVM_counts)
UVM_counts <- UVM_counts %>%
  dplyr::select(bcr_patient_barcode, riskscore)

UVM_clinicdata <- clinicaldata %>%
  dplyr::filter(clinicaldata$bcr_patient_barcode %in% UVM_counts$bcr_patient_barcode)

UVM_clinicdata <- UVM_clinicdata[match(rownames(UVM_counts),UVM_clinicdata$bcr_patient_barcode),]  ## kp ob man das brauchr aber sicher ist sicher

UVM_risk <- left_join(UVM_clinicdata, UVM_counts, by = "bcr_patient_barcode")
UVM_risk$OS <- as.numeric(UVM_risk$OS)
UVM_risk$OS.time <- as.numeric(UVM_risk$OS.time)

stat <- maxstat.test(Surv(UVM_risk$OS.time , UVM_risk$OS) ~ UVM_risk$riskscore, data = UVM_risk, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(UVM_risk$riskscore > cutoff,"high","low"))
UVM_risk$risk <- risk
## Kaplanmeier
UVM_risk$OS.time <- UVM_risk$OS.time / 365 * 12
surv_obj <- Surv(UVM_risk$OS.time , UVM_risk$OS)
fit <- survfit(surv_obj~risk, data = UVM_risk) 
summary(fit)
ggsurvplot(fit, 
           title = "UVM",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,
           risk.table = T,
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))
## univriate
risk_cancerproof <- UVM_risk 
risk_cancerproof <- risk_cancerproof %>%
  mutate(risk_cancerproof, risk = case_when(risk == "low" ~ "1",
                                            risk == "high" ~ "2"))
res.cox <- coxph(Surv(OS.time, OS) ~ risk, data = risk_cancerproof)
res.cox
summary(res.cox)





#####Forestplots####
## forestplot mit den KM ergebnissen
library(forestplot)

forestplot_data <- structure(list(
  mean = c(NA, 0.3776,
           1.922,
           1.995,
           4.308,
           5.274,
           1.395,
           8.188,
           0.5924,
           1.744,
           2.235,
           3.333,
           0.546,
           0.3268,
           1.658,
           2.065,
           1.504,
           1.489,
           2.584,
           1.42,
           2.344,
           0.3766,
           3.118,
           0.6443,
           2.288,
           1.822,
           1.442,
           0.1449,
           4.903,
           14.35,
           2.262,
           2.064,
           0.2744),
  lower = c(NA, 0.1699,
            1.429,
            1.218,
            2.201,
            1.2,
            0.931,
            1.65,
            0.2918,
            1.163,
            1.645,
            0.8937,
            0.4047,
            0.1506,
            1.121,
            1.343,
            1.103,
            1.136,
            1.506,
            1.083,
            1.508,
            0.06286,
            0.5558,
            0.1918,
            1.495,
            0.8662,
            0.9827,
            0.0131,
            1.113,
            2.884,
            1.46,
            0.9572,
            0.112),
  upper = c(NA, 0.8391,
            2.586,
            3.265,
            8.433,
            23.17,
            2.089,
            40.65,
            1.203,
            2.614,
            3.038,
            12.43,
            0.7365,
            0.709,
            2.453,
            3.176,
            2.051,
            1.952,
            4.432,
            1.863,
            3.643,
            2.257,
            17.5,
            2.164,
            3.504,
            3.834,
            2.117,
            1.602,
            21.61,
            71.43,
            3.503,
            4.452,
            0.6728)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -33L),
  class = "data.frame"
)

tabletext <- cbind(
  c("Tumortype", "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS","UVM"),
  c("p Value", "2.56E-02",
    "2.35E-05",
    "2.78E-03",
"6.58E-07",
"7.68E-03",
"1.13E-01",
"1.72E-02",
"1.24E-01",
"6.64E-03",
"5.59E-08",
"7.51E-02",
"8.90E-05",
"1.19E-02",
"1.48E-02",
"4.21E-04",
"8.47E-03",
"4.25E-03",
"2.92E-04",
"1.01E-02",
"3.02E-04",
"3.09E-01",
"1.60E-01",
"5.01E-01",
"3.14E-04",
"1.23E-01",
"5.40E-02",
"1.03E-01",
"1.19E-02",
"2.39E-04",
"1.55E-04",
"5.26E-02",
"2.90E-03")
)

forestplot(tabletext,
           boxsize = 0.1,
           hrzl_lines = gpar(col = "#444444"),
           forestplot_data,new_page = TRUE,
           is.summary = c(TRUE,rep(FALSE,32)),
           #clip = c(0.1,1.0), 
           xlog = TRUE, 
           xlab = "Hazard ratio",
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           vertices = TRUE)

## forestplot mit den univariate ergebnissen sortiert nach HR
library(forestplot)

forestplot_data <- structure(list(
  mean = c(NA, 14.35,
           8.188,
           5.274,
           4.903,
           4.308,
           3.333,
           3.118,
           2.584,
           2.344,
           2.288,
           2.262,
           2.235,
           2.065,
           2.064,
           1.995,
           1.922,
           1.822,
           1.744,
           1.658,
           1.504,
           1.489,
           1.442,
           1.42,
           1.395,
           0.6443,
           0.5924,
           0.546,
           0.3776,
           0.3766,
           0.3268,
           0.2744,
           0.1449),
  lower = c(NA, 2.884,
            1.65,
            1.2,
            1.113,
            2.201,
            0.8937,
            0.5558,
            1.506,
            1.508,
            1.495,
            1.46,
            1.645,
            1.343,
            0.9572,
            1.218,
            1.429,
            0.8662,
            1.163,
            1.121,
            1.103,
            1.136,
            0.9827,
            1.083,
            0.931,
            0.1918,
            0.2918,
            0.4047,
            0.1699,
            0.06286,
            0.1506,
            0.112,
            0.0131),
  upper = c(NA, 71.43,
            40.65,
            23.17,
            21.61,
            8.433,
            12.43,
            17.5,
            4.432,
            3.643,
            3.504,
            3.503,
            3.038,
            3.176,
            4.452,
            3.265,
            2.586,
            3.834,
            2.614,
            2.453,
            2.051,
            1.952,
            2.117,
            1.863,
            2.089,
            2.164,
            1.203,
            0.7365,
            0.8391,
            2.257,
            0.709,
            0.6728,
            1.602)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -33L),
  class = "data.frame"
)

tabletext <- cbind(
  c("Tumortype", "THYM",
    "DLBC",
    "CHOL",
    "THCA",
    "CESC",
    "KICH",
    "PRAD",
    "MESO",
    "PAAD",
    "SARC",
    "UCEC",
    "HNSC",
    "LIHC",
    "UCS",
    "BRCA",
    "BLCA",
    "SKCM",
    "GBM",
    "LGG",
    "LUAD",
    "LUSC",
    "STAD",
    "OV",
    "COAD",
    "READ",
    "ESCA",
    "KIRC",
    "ACC",
    "PCPG",
    "KIRP",
    "UVM",
    "TGCT"),
  c("p Value", "3.00E-05",
    "2.00E-03",
    "1.00E-02",
    "2.00E-02",
    "4.00E-06",
    "6.00E-02",
    "2.00E-01",
    "4.00E-04",
    "1.00E-04",
    "9.00E-05",
    "2.00E-04",
    "1.00E-07",
    "7.00E-04",
    "6.00E-02",
    "5.00E-03",
    "1.00E-05",
    "1.00E-01",
    "6.00E-03",
    "1.00E-02",
    "9.00E-03",
    "4.00E-03",
    "6.00E-02",
    "1.00E-02",
    "1.00E-01",
    "5.00E-01",
    "1.00E-01",
    "6.00E-05",
    "1.00E-02",
    "3.00E-01",
    "3.00E-03",
    "3.00E-03",
    "7.00E-02")
)

forestplot(tabletext,
           boxsize = 0.1,
           hrzl_lines = gpar(col = "#444444"),
           forestplot_data,new_page = TRUE,
           is.summary = c(TRUE,rep(FALSE,32)),
           #clip = c(0.1,1.0), 
           xlog = TRUE, 
           xlab = "Hazard ratio",
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           vertices = TRUE)

## forestplot mit den KM ergebnissen sortiert nach pvalue
library(forestplot)

forestplot_data <- structure(list(
  mean = c(NA, 2.235,
           4.308,
           1.922,
           0.546,
           2.262,
           14.35,
           2.584,
           2.344,
           2.288,
           2.065,
           1.995,
           0.2744,
           1.489,
           1.744,
           5.274,
           1.504,
           1.42,
           4.903,
           0.3268,
           1.658,
           8.188,
           0.3776,
           2.064,
           1.442,
           3.333,
           0.1449,
           1.395,
           1.822,
           0.5924,
           3.118,
           0.3766,
           0.6443),
  lower = c(NA, 1.645,
            2.201,
            1.429,
            0.4047,
            1.46,
            2.884,
            1.506,
            1.508,
            1.495,
            1.343,
            1.218,
            0.112,
            1.136,
            1.163,
            1.2,
            1.103,
            1.083,
            1.113,
            0.1506,
            1.121,
            1.65,
            0.1699,
            0.9572,
            0.9827,
            0.8937,
            0.0131,
            0.931,
            0.8662,
            0.2918,
            0.5558,
            0.06286,
            0.1918),
  upper = c(NA, 3.038,
            8.433,
            2.586,
            0.7365,
            3.503,
            71.43,
            4.432,
            3.643,
            3.504,
            3.176,
            3.265,
            0.6728,
            1.952,
            2.614,
            23.17,
            2.051,
            1.863,
            21.61,
            0.709,
            2.453,
            40.65,
            0.8391,
            4.452,
            2.117,
            12.43,
            1.602,
            2.089,
            3.834,
            1.203,
            17.5,
            2.257,
            2.164)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -33L),
  class = "data.frame"
)

tabletext <- cbind(
  c("Tumortype", "HNSC",
    "CESC",
    "BLCA",
    "KIRC",
    "UCEC",
    "THYM",
    "MESO",
    "PAAD",
    "SARC",
    "LIHC",
    "BRCA",
    "UVM",
    "LUSC",
    "GBM",
    "CHOL",
    "LUAD",
    "OV",
    "THCA,",
    "KIRP",
    "LGG",
    "DLBC",
    "ACC",
    "UCS",
    "STAD",
    "KICH",
    "TGCT",
    "COAD",
    "SKCM",
    "ESCA",
    "PRAD",
    "PCPG",
    "READ"),
  c("p Value", "5.59E-08",
    "6.58E-07",
    "2.35E-05",
    "8.90E-05",
    "1.55E-04",
    "2.39E-04",
    "2.92E-04",
    "3.02E-04",
    "3.14E-04",
    "4.21E-04",
    "2.78E-03",
    "2.90E-03",
    "4.25E-03",
    "6.64E-03",
    "7.68E-03",
    "8.47E-03",
    "1.01E-02",
    "1.19E-02",
    "1.19E-02",
    "1.48E-02",
    "1.72E-02",
    "2.56E-02",
    "5.26E-02",
    "5.40E-02",
    "7.51E-02",
    "1.03E-01",
    "1.13E-01",
    "1.23E-01",
    "1.24E-01",
    "1.60E-01",
    "3.09E-01",
    "5.01E-01")
)

forestplot(tabletext,
           boxsize = 0.1,
           hrzl_lines = gpar(col = "#444444"),
           forestplot_data,new_page = TRUE,
           is.summary = c(TRUE,rep(FALSE,32)),
           #clip = c(0.1,1.0), 
           xlog = TRUE, 
           xlab = "Hazard ratio",
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           vertices = TRUE)

