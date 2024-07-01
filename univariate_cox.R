library(glmnet)
library(survival)
library(readxl) 
library(readr) 
library(maxstat)
library(org.Hs.eg.db)
library(survminer)

##Univariate cox regression for Riskmodel
## OS_5y mit den drei Endpunkten (OS; DSS; PFI)
## Daten einlesen
OS_5y_riskmodel <- read_excel("~/Desktop/Doktorarbeit/Data/LassoRiskscore_fpkm.xlsx", sheet = 1) 
OS_5y_riskmodel$Patient_ID <- OS_5y_riskmodel$id
OS_5y_riskmodel <- OS_5y_riskmodel %>%
  dplyr::select(Patient_ID, risk)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)
## Daten zusammenfügen für Analyse
OS_5y_riskmodeldata <- left_join(TCGA_clinicaldata, OS_5y_riskmodel, by = "Patient_ID", copy = TRUE)
OS_5y_riskmodeldata$OS_5y_event <- as.numeric(OS_5y_riskmodeldata$OS_5y_event)
OS_5y_riskmodeldata$DSS_5y_event <- as.numeric(OS_5y_riskmodeldata$DSS_5y_event)
OS_5y_riskmodeldata$PFI_5y_event <- as.numeric(OS_5y_riskmodeldata$PFI_5y_event)
## OS model mit OS Endpunkt
res.cox_OS_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ risk, data = OS_5y_riskmodeldata)
res.cox_OS_OS
summary(res.cox_OS_OS)
## OS model mit DSS Endpunkt
res.cox_OS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ risk, data = OS_5y_riskmodeldata)
res.cox_OS_DSS
summary(res.cox_OS_DSS)
## OS model mit PFI Endpunkt
res.cox_OS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ risk, data = OS_5y_riskmodeldata)
res.cox_OS_PFI
summary(res.cox_OS_PFI)


##Univariate cox regression for Riskmodel
## OS_5y mit den drei Endpunkten (OS; DSS; PFI)
## Daten einlesen
DSS_5y_riskmodel <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1) 
DSS_5y_riskmodel$Patient_ID <- DSS_5y_riskmodel$id
DSS_5y_riskmodel <- DSS_5y_riskmodel %>%
  dplyr::select(Patient_ID, risk)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)
## Daten zusammenfügen für Analyse
DSS_5y_riskmodeldata <- left_join(TCGA_clinicaldata, DSS_5y_riskmodel, by = "Patient_ID", copy = TRUE)
DSS_5y_riskmodeldata$OS_5y_event <- as.numeric(DSS_5y_riskmodeldata$OS_5y_event)
DSS_5y_riskmodeldata$DSS_5y_event <- as.numeric(DSS_5y_riskmodeldata$DSS_5y_event)
DSS_5y_riskmodeldata$PFI_5y_event <- as.numeric(DSS_5y_riskmodeldata$PFI_5y_event)
## DSS model mit OS Endpunkt
res.cox_DSS_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_OS
summary(res.cox_DSS_OS)
## DSS model mit DSS Endpunkt
res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)
## DSS model mit PFI Endpunkt
res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)


##Ergebnisse Forestplotten

library(forestplot)

forestplot_data <- structure(list(
  mean = c(NA, 0.4102, 0.4916, 0.7056, 0.3989, 0.3141, 0.5193),
  lower = c(NA, 0.2928, 0.3254, 0.5111, 0.2903, 0.2155, 0.3768),
  upper = c(NA, 0.5746, 0.7426, 0.9742, 0.5482, 0.4577, 0.7159)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -7L),
  class = "data.frame"
)

tabletext <- cbind(
  c("Riskmodel/Endpunkt", "OS_5y/OS_5y", "OS_5y/DSS_5y", "OS_5y/PFI_5y", "DSS_5y/OS_5y", "DSS_5y/DSS_5y", "DSS_5y/PFI_5y"),
  c("Hazard Ratio", "0.4102", "0.4916", "0.7056", "0.3989", "0.3141", "0.5193")
)

forestplot(tabletext, 
           forestplot_data,new_page = TRUE,
           is.summary = c(TRUE,rep(FALSE,6)),
           clip = c(0.1,1.0), 
           xlog = TRUE, 
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           vertices = TRUE)


##Ergebnisse Forestplotten

library(forestplot)

forestplot_data <- structure(list(
  mean = c(NA, 2.4378, 2.0342, 1.4172, 2.5069, 3.1837, 1.9257),
  lower = c(NA, 1.7403, 1.3466, 1.0265, 1.8242, 2.1848, 1.3968),
  upper = c(NA, 3.4153, 3.0731, 1.9566, 3.4447, 4.6404, 2.6539)),
  .Names = c("mean", "lower", "upper"),
  row.names = c(NA, -7L),
  class = "data.frame"
)

tabletext <- cbind(
  c("Risk model / Endpoint", "OS_5y/OS_5y", "OS_5y/DSS_5y", "OS_5y/PFI_5y", "DSS_5y/OS_5y", "DSS_5y/DSS_5y", "DSS_5y/PFI_5y"),
  c("p Value", "3E-08", "4E-04", "3E-02", "4E-08", "4E-09", "1E-04")
)

forestplot(tabletext,
           boxsize = 0.1,
           hrzl_lines = gpar(col = "#444444"),
           forestplot_data,new_page = TRUE,
           is.summary = c(TRUE,rep(FALSE,6)),
           #clip = c(0.1,1.0), 
           xlog = TRUE, 
           xlab = "Hazard ratio",
           col = fpColors(box = "royalblue",
                          line = "darkblue",
                          summary = "royalblue"),
           vertices = TRUE)

###Univariate and multivariate Cox regression analysis for overall survival of the TCGA-HNSC cohort
## den bereits vorher gemachten Datensatz benutzen und die variablen in binäre optionen ändern
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, Age_61 = case_when(Age > 61 ~ "No",
                                                  Age <= 61 ~ "Yes"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%  
  mutate(DSS_5y_riskmodeldata, cN = case_when(cN == "N0" ~ "cN0",
                                              cN != "N0" ~ "cN+"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%  
  mutate(DSS_5y_riskmodeldata, pN = case_when(pN == "N0" ~ "pN0",
                                              pN != "N0" ~ "pN+"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%  
  mutate(DSS_5y_riskmodeldata, cM = case_when(cM == "M0" ~ "cM0",
                                              cM != "M0" ~ "cM+"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%  
  mutate(DSS_5y_riskmodeldata, pM = case_when(pM == "M0" ~ "pM0",
                                              pM != "M0" ~ "pM+"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, cT = case_when(cT == "T1" ~ "cT1/2",
                                              cT == "T2" ~ "cT1/2",
                                              cT == "T3" ~ "cT3/4",
                                              cT == "TX" ~ "cT3/4",
                                              cT == "T4" ~ "cT3/4",
                                              cT == "T4a" ~ "cT3/4",
                                              cT == "T4b" ~ "cT3/4"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, pT = case_when(pT == "T0" ~ "pT1/2",
                                              pT == "T1" ~ "pT1/2",
                                              pT == "T2" ~ "pT1/2",
                                              pT == "T3" ~ "pT3/4",
                                              pT == "TX" ~ "pT3/4",
                                              pT == "T4" ~ "pT3/4",
                                              pT == "T4a" ~ "pT3/4",
                                              pT == "T4b" ~ "pT3/4"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, Grading = case_when(Grading == "G1" ~ "G1/2",
                                                   Grading == "G2" ~ "G1/2",
                                                   Grading == "G3" ~ "G3/4",
                                                   Grading == "G4" ~ "G3/4",
                                                   Grading == "GX" ~ "G3/4"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, margin_status = case_when(margin_status == "Negative" ~ "R0",
                                                         margin_status == "Positive" ~ "R+",
                                                         margin_status == "Close" ~ "R+"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, Gender = case_when(Gender == "Male" ~ "1",
                                                  Gender == "Female" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  #mutate(DSS_5y_riskmodeldata, Alcohol = case_when(Alcohol == "Yes" ~ "1",
                                                  # Alcohol == "No" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
 # mutate(DSS_5y_riskmodeldata, HPV16 = case_when(HPV16 == "Positive" ~ "1",
                                               #  HPV16 == "Negative" ~ "2"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, cN = case_when(cN == "cN0" ~ "1",
                                              cN == "cN+" ~ "2"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, pN = case_when(pN == "pN0" ~ "1",
                                              pN == "pN+" ~ "2"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, risk = case_when(risk == "low" ~ "1",
                                                risk == "high" ~ "2"))
DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
  mutate(DSS_5y_riskmodeldata, margin_status = case_when(margin_status == "R0" ~ "1",
                                                         margin_status == "R+" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
#  mutate(DSS_5y_riskmodeldata, pT = case_when(pT == "pT3/4" ~ "1",
#                                              pT == "pT1/2" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
#  mutate(DSS_5y_riskmodeldata, Grading = case_when(Grading == "G3/4" ~ "1",
#                                                   Grading == "G1/2" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
#  mutate(DSS_5y_riskmodeldata, ALI = case_when(ALI == "Yes" ~ "1",
#                                                   ALI == "No" ~ "2"))
#DSS_5y_riskmodeldata <- DSS_5y_riskmodeldata %>%
# mutate(DSS_5y_riskmodeldata, PNI = case_when(PNI == "Yes" ~ "1",
#                                                  PNI == "No" ~ "2"))


## univariate cox regression für jede variable

## GENDER

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Gender, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Gender, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Gender, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## ALTER

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Age_61, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Age_61, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Age_61, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## RAUCHEN

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Smoking, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Smoking, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Smoking, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## ALK

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Alcohol, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Alcohol, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Alcohol, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## HPV16

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ HPV16, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ HPV16, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ HPV16, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## TUMORGRÖßE

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ cT, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ cT, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ cT, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)



res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ pT, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ pT, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ pT, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## LYMPHKNOTEN

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ cN, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ cN, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ cN, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)


res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ pN, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ pN, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ pN, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## METASTASEN

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ cM, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ cM, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ cM, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)


res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ pM, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ pM, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ pM, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## GRADING

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Grading, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Grading, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Grading, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## MARGIN

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ margin_status, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ margin_status, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ margin_status, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## ALI

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ ALI, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ ALI, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ ALI, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## PNI

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ PNI, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ PNI, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ PNI, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## RISK

res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)

## Die Variablen die p<0.05 auf unabhängigkeit testen mittels multivariate cox regression
res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Gender + Age_61 + HPV16 + pT + pN + margin_status + ALI + PNI + risk, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ HPV16 + pT + pN + margin_status + PNI + risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ HPV16 + pT + pN + margin_status + PNI + risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)


## Die Variablen die n>400 haben auf unabhängigkeit checken
res.cox_OS <- coxph(Surv(OS_5y_months, OS_5y_event) ~ Gender + Age_61 + HPV16 + cT + pT + cN + pN + Smoking + Alcohol + Grading + margin_status + risk, data = DSS_5y_riskmodeldata)
res.cox_OS
summary(res.cox_OS)

res.cox_DSS_DSS <- coxph(Surv(DSS_5y_months, DSS_5y_event) ~ Gender + Age_61 + HPV16 + cT + pT + cN + pN + Smoking + Alcohol + Grading + margin_status + risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_DSS
summary(res.cox_DSS_DSS)

res.cox_DSS_PFI <- coxph(Surv(PFI_5y_months, PFI_5y_event) ~ Gender + Age_61 + HPV16 + cT + pT + cN + pN + Smoking + Alcohol + Grading + margin_status + risk, data = DSS_5y_riskmodeldata)
res.cox_DSS_PFI
summary(res.cox_DSS_PFI)


###Univariate and multivariate Cox regression analysis for overall survival of the CPTAC-HNSC cohort
## den bereits vorher gemachten Datensatz benutzen und die variablen in binäre optionen ändern

contr.grp.clinical <- read_delim("~/Desktop/Doktorarbeit/Data/HS_CPTAC_HNSCC_CLI.tsi.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
contr.grp.clinical <- contr.grp.clinical %>%
  dplyr::rename(`Patient_ID` = `case_id`)
contr.grp.clinical <- contr.grp.clinical %>%
  filter(contr.grp.clinical$Patient_ID %in% contr.grp.cox$Patient_ID)

contr.grp.cox <- contr.grp.all %>%
  dplyr::select(Patient_ID, risk, riskscore)

contr.grp.cox <- left_join(contr.grp.cox, contr.grp.clinical, by = "Patient_ID")
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, Age_61 = case_when(age > 61 ~ "No",
                                           age <= 61 ~ "Yes"))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, Smoking = case_when(smoke_age_start == 38 ~ "Yes",
                                            smoke_age_start == 37 ~ "Yes",
                                            smoke_age_start == 36 ~ "Yes",
                                            smoke_age_start == 35 ~ "Yes",
                                            smoke_age_start == 34 ~ "Yes",
                                            smoke_age_start == 33 ~ "Yes",
                                            smoke_age_start == 32 ~ "Yes",
                                            smoke_age_start == 31 ~ "Yes",
                                            smoke_age_start == 30 ~ "Yes",
                                            smoke_age_start == 29 ~ "Yes",
                                            smoke_age_start == 28 ~ "Yes",
                                            smoke_age_start == 27 ~ "Yes",
                                            smoke_age_start == 26 ~ "Yes",
                                            smoke_age_start == 25 ~ "Yes",
                                            smoke_age_start == 24 ~ "Yes",
                                            smoke_age_start == 23 ~ "Yes",
                                            smoke_age_start == 22 ~ "Yes",
                                            smoke_age_start == 21 ~ "Yes",
                                            smoke_age_start == 20 ~ "Yes",
                                            smoke_age_start == 19 ~ "Yes",
                                            smoke_age_start == 18 ~ "Yes",
                                            smoke_age_start == 17 ~ "Yes",
                                            smoke_age_start == 16 ~ "Yes",
                                            smoke_age_start == 15 ~ "Yes",
                                            smoke_age_start == 14 ~ "Yes",
                                            smoke_age_start == 13 ~ "Yes",
                                            smoke_age_start == 12 ~ "Yes",
                                            smoke_age_start == 11 ~ "Yes",
                                            smoke_age_start == 10 ~ "Yes",
                                            smoke_age_start == 9 ~ "Yes",
                                            smoke_age_start == 8 ~ "Yes",
                                            smoke_age_start == "Unknown" ~ "Yes",
                                            is.na(smoke_age_start) ~ "No"))

contr.grp.cox$alcohol_consum <- as.character(contr.grp.cox$alcohol_consum)
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, Alcohol = case_when(alcohol_consum == "Lifelong non-drinker" ~ "No",
                                           alcohol_consum == "Alcohol consumption equal to or less than 2 drinks per day for men and 1 drink or less per day for women" ~ "Yes",
                                           alcohol_consum == "Alcohol consumption more than 2 drinks per day for men and more than 1 drink per day for women" ~ "Yes",
                                           alcohol_consum == "Consumed alcohol in the past, but currently a non-drinker" ~ "Yes",
                                           #alcohol_consum == "Alcohol consumption history not available" ~ "NA"
                                           ))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, pT = case_when(       patho_staging_pt == "pT1" ~ "pT1/2",
                                              patho_staging_pt == "pT1b" ~ "pT1/2",
                                              patho_staging_pt == "pT2" ~ "pT1/2",
                                              patho_staging_pt == "T2" ~ "pT1/2",
                                              patho_staging_pt == "pT3" ~ "pT3/4",
                                              patho_staging_pt == "T4a" ~ "pT3/4",
                                              patho_staging_pt == "pT4" ~ "pT3/4",
                                              patho_staging_pt == "pT4a" ~ "pT3/4",
                                              patho_staging_pt == "pT4b" ~ "pT3/4"))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, pN = case_when(       patho_staging_pn == "pN0" ~ "pN0",
                                              patho_staging_pn != "pN0" ~ "pN+"))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, Grade = case_when(  histologic_grade == "G1 Well differentiated" ~ "G1/2",
                                            histologic_grade == "G2 Moderately differentiated" ~ "G1/2",
                                            histologic_grade == "G3 Poorly differentiated" ~ "G3"))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, gender = case_when(gender == "Male" ~ "1",
                                                  gender == "Female" ~ "2"))
contr.grp.cox <- contr.grp.cox %>%
  mutate(contr.grp.cox, risk = case_when(risk == "low" ~ "1",
                                                risk == "high" ~ "2"))
contr.grp.cox$overall_survival <- as.numeric(contr.grp.cox$overall_survival)
contr.grp.cox$overall_free_status <- as.numeric(contr.grp.cox$overall_free_status)
#univariaten test

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ gender, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ Age_61, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ Smoking, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ Alcohol, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ pT, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ pN, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ Grade, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ risk, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

## multi mit sig var

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ Grade + risk, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)

## mulit mit n var

res.cox_OS <- coxph(Surv(contr.grp.cox$overall_survival, contr.grp.cox$overall_free_status) ~ gender + Age_61+ Smoking + pT + pN +Grade + risk, data = contr.grp.cox)
res.cox_OS
summary(res.cox_OS)