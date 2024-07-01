library(glmnet)
library(survival)
library(readxl) 
library(maxstat)
library(org.Hs.eg.db)


#  DATEN EINLESEN
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
relevant_genelist_egfr <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)

# DATENSATZ VORBERITEN

TCGA_EGFR_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, OS_5y_event, OS_5y_months)

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS_5y_months
TCGA_EGFR_lasso$OS <- TCGA_EGFR_lasso$OS_5y_event
TCGA_EGFR_lasso$OS_5y_event <- NULL
TCGA_EGFR_lasso$OS_5y_months <- NULL
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  #dplyr::select(Patient_ID, OS.time, OS)


fpkm_data <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% relevant_genelist_egfr$GeneID_all) %>%
  as.data.frame(.)
colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Gene_ID <- rownames(fpkm_data)
#fpkm_data$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = fpkm_data$Gene_ID, column= "SYMBOL", keytype =  "ENTREZID")
#fpkm_data <- fpkm_data %>%
  #filter(!is.na(fpkm_data$Gene_Symbol))
#rownames(fpkm_data) <- fpkm_data$Gene_Symbol
#fpkm_data$Gene_ID <- NULL
#fpkm_data$Gene_Symbol <- NULL
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Patient_ID <- c(rownames(fpkm_data))
#rownames(fpkm_data) <- NULL

TCGA_EGFR_lasso <- TCGA_EGFR_lasso[match(rownames(fpkm_data),TCGA_EGFR_lasso$Patient_ID),]
TCGA_EGFR_lasso <- left_join(TCGA_EGFR_lasso, fpkm_data, by = "Patient_ID", copy = TRUE) #%>%
 # filter(TCGA_EGFR_lasso$Patient_ID %in% fpkm_data$Patient_ID)
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  #filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso <- as.data.frame(TCGA_EGFR_lasso)
#rownames(TCGA_EGFR_lasso) <- TCGA_EGFR_lasso$Patient_ID
#TCGA_EGFR_lasso$Patient_ID <- NULL

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/12
#TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/365

TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  filter(!is.na(TCGA_EGFR_lasso$OS))
TCGA_EGFR_lasso$OS <- as.numeric(TCGA_EGFR_lasso$OS)

x <- as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)])
y <- data.matrix(Surv(TCGA_EGFR_lasso$OS.time, TCGA_EGFR_lasso$OS))


fit <- glmnet(x, y, family = "cox")
plot(fit,  xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
#write.table

####TEST#####
test_risk <- TCGA_EGFR_lasso %>%
  dplyr::select(Patient_ID, OS.time, OS, '7535', '84215', '1948', '147372', '645432')
test_risk$riskscore <- test_risk$'7535' * (-0.00893014079168287) + test_risk$'84215' * -0.0245566035731716 + test_risk$'1948' * 0.00682598618190403 + test_risk$'147372' * 0.0527359814359829+ test_risk$'645432' * -0.484764320166299

#0.53818275*-0.00893014079168287 + 0.114368264 * -0.0245566035731716 + 24.0633075 * 0.00682598618190403 + 2.355539341 * 0.0527359814359829 + 0.054170154 * -0.484764320166299
######
riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
standard_data <- c("Patient_ID", "OS.time", "OS",lassoGene)
final_riskscore <- cbind(TCGA_EGFR_lasso[,outCol],riskScore=as.vector(riskScore))


#######

cvcoxlasso <- cv.glmnet(
  x = as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = y,
  family = 'cox')

plot(cvcoxlasso)
cvcoxlasso$lambda.1se

coef(cvcoxlasso)

plot(survival::survfit(cvcoxlasso, s = cvcoxlasso$lambda.1se, x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]), y = b))

coxlasso <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox')

plot(coxlasso)

lasso_fit <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox',
  alpha = 1,
  lambda = cv.coxlasso$lambda.1se
)
 plot(lasso_fit)


#######
lasso_result <- test_risk
rownames(lasso_result) <- lasso_result$Patient_ID
lasso_result$Patient_ID <- NULL

stat <- maxstat.test(Surv(lasso_result$OS.time,lasso_result$OS)~lasso_result$riskscore, data = lasso_result, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(lasso_result$riskscore > cutoff,"high","low"))
write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
            file="lasso-Risk-Cox-best-cut-5y_OS.txt",
            sep="\t",
            quote=F,
            row.names=F)



########### Lasso Cox mit DSS
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
relevant_genelist_egfr <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)

# DATENSATZ VORBERITEN

TCGA_EGFR_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, DSS, DSS_days)

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$DSS_days
TCGA_EGFR_lasso$OS <- TCGA_EGFR_lasso$DSS
TCGA_EGFR_lasso$DSS <- NULL
TCGA_EGFR_lasso$DSS_days <- NULL
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#dplyr::select(Patient_ID, OS.time, OS)


fpkm_data <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% relevant_genelist_egfr$GeneID_all) %>%
  as.data.frame(.)
colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))

rownames(fpkm_data) <- NULL

TCGA_EGFR_lasso <- TCGA_EGFR_lasso[match(rownames(fpkm_data),TCGA_EGFR_lasso$Patient_ID),]
TCGA_EGFR_lasso <- left_join(TCGA_EGFR_lasso, fpkm_data, by = "Patient_ID", copy = TRUE) #%>%

TCGA_EGFR_lasso <- as.data.frame(TCGA_EGFR_lasso)


TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/12


TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  filter(!is.na(TCGA_EGFR_lasso$OS))
TCGA_EGFR_lasso$OS <- as.numeric(TCGA_EGFR_lasso$OS)

x <- as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)])
y <- data.matrix(Surv(TCGA_EGFR_lasso$OS.time, TCGA_EGFR_lasso$OS))


fit <- glmnet(x, y, family = "cox")
plot(fit,  xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
#write.table

####TEST#####
test_risk <- TCGA_EGFR_lasso %>%
  dplyr::select(Patient_ID, OS.time, OS, '7535', '105373853', '1948', '147372', '645432')
test_risk$riskscore <- test_risk$'7535' * (-0.0555160154489259) + test_risk$'105373853' * 0.00181181168420555 + test_risk$'1948' * 0.00790511531356365 + test_risk$'147372' * 0.116634087076066+ test_risk$'645432' * -0.033711373761946

#0.53818275*-0.00893014079168287 + 0.114368264 * -0.0245566035731716 + 24.0633075 * 0.00682598618190403 + 2.355539341 * 0.0527359814359829 + 0.054170154 * -0.484764320166299
######
riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
standard_data <- c("Patient_ID", "OS.time", "OS",lassoGene)
final_riskscore <- cbind(TCGA_EGFR_lasso[,outCol],riskScore=as.vector(riskScore))


#######

cvcoxlasso <- cv.glmnet(
  x = as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = y,
  family = 'cox')

plot(cvcoxlasso)
cvcoxlasso$lambda.1se

coef(cvcoxlasso)

plot(survival::survfit(cvcoxlasso, s = cvcoxlasso$lambda.1se, x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]), y = b))

coxlasso <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox')

plot(coxlasso)

lasso_fit <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox',
  alpha = 1,
  lambda = cv.coxlasso$lambda.1se
)
plot(lasso_fit)


#######
lasso_result <- test_risk
rownames(lasso_result) <- lasso_result$Patient_ID
lasso_result$Patient_ID <- NULL

stat <- maxstat.test(Surv(lasso_result$OS.time,lasso_result$OS)~lasso_result$riskscore, data = lasso_result, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(lasso_result$riskscore > cutoff,"high","low"))
write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
            file="lasso-Risk-Cox-best-cut-DSS.txt",
            sep="\t",
            quote=F,
            row.names=F)
### DSS 5y
######
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
relevant_genelist_egfr <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)

# DATENSATZ VORBERITEN

TCGA_EGFR_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, DSS_5y_event, DSS_5y_months)

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$DSS_5y_months
TCGA_EGFR_lasso$OS <- TCGA_EGFR_lasso$DSS_5y_event
TCGA_EGFR_lasso$DSS_5y_event <- NULL
TCGA_EGFR_lasso$DSS_5y_months <- NULL
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#dplyr::select(Patient_ID, OS.time, OS)


fpkm_data <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% relevant_genelist_egfr$GeneID_all) %>%
  as.data.frame(.)
colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Gene_ID <- rownames(fpkm_data)
#fpkm_data$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = fpkm_data$Gene_ID, column= "SYMBOL", keytype =  "ENTREZID")
#fpkm_data <- fpkm_data %>%
#filter(!is.na(fpkm_data$Gene_Symbol))
#rownames(fpkm_data) <- fpkm_data$Gene_Symbol
#fpkm_data$Gene_ID <- NULL
#fpkm_data$Gene_Symbol <- NULL
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Patient_ID <- c(rownames(fpkm_data))
rownames(fpkm_data) <- NULL

TCGA_EGFR_lasso <- TCGA_EGFR_lasso[match(rownames(fpkm_data),TCGA_EGFR_lasso$Patient_ID),]
TCGA_EGFR_lasso <- left_join(TCGA_EGFR_lasso, fpkm_data, by = "Patient_ID", copy = TRUE) #%>%
# filter(TCGA_EGFR_lasso$Patient_ID %in% fpkm_data$Patient_ID)
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso <- as.data.frame(TCGA_EGFR_lasso)
#rownames(TCGA_EGFR_lasso) <- TCGA_EGFR_lasso$Patient_ID
#TCGA_EGFR_lasso$Patient_ID <- NULL

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/12
#TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/365

TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  filter(!is.na(TCGA_EGFR_lasso$OS))
TCGA_EGFR_lasso$OS <- as.numeric(TCGA_EGFR_lasso$OS)

x <- as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)])
y <- data.matrix(Surv(TCGA_EGFR_lasso$OS.time, TCGA_EGFR_lasso$OS))


fit <- glmnet(x, y, family = "cox")
plot(fit,  xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
#write.table

####TEST#####
test_risk <- TCGA_EGFR_lasso %>%
  dplyr::select(Patient_ID, OS.time, OS, '7535', '494514', '1948', '147372', '645432', '51806', '56936')
test_risk$riskscore <- test_risk$'7535' * (-0.0582152196382938) + test_risk$'494514' * 0.013184642798603 + test_risk$'1948' * 0.00470780696956776 + test_risk$'147372' * 0.11219886435049+ test_risk$'645432' * -0.267800991188706 + test_risk$'51806' * -7.50310759764005e-05 + test_risk$'56936' * -0.0844511106749959

#0.53818275*-0.00893014079168287 + 0.114368264 * -0.0245566035731716 + 24.0633075 * 0.00682598618190403 + 2.355539341 * 0.0527359814359829 + 0.054170154 * -0.484764320166299
######
riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
standard_data <- c("Patient_ID", "OS.time", "OS",lassoGene)
final_riskscore <- cbind(TCGA_EGFR_lasso[,outCol],riskScore=as.vector(riskScore))


#######

cvcoxlasso <- cv.glmnet(
  x = as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = y,
  family = 'cox')

plot(cvcoxlasso)
cvcoxlasso$lambda.1se

coef(cvcoxlasso)

plot(survival::survfit(cvcoxlasso, s = cvcoxlasso$lambda.1se, x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]), y = b))

coxlasso <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox')

plot(coxlasso)

lasso_fit <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox',
  alpha = 1,
  lambda = cv.coxlasso$lambda.1se
)
plot(lasso_fit)


#######
lasso_result <- test_risk
rownames(lasso_result) <- lasso_result$Patient_ID
lasso_result$Patient_ID <- NULL

stat <- maxstat.test(Surv(lasso_result$OS.time,lasso_result$OS)~lasso_result$riskscore, data = lasso_result, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(lasso_result$riskscore > cutoff,"high","low"))
write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
            file="lasso-Risk-Cox-best-cut-DSS-5y.txt",
            sep="\t",
            quote=F,
            row.names=F)

######PFI####
#Daten einlesen
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
relevant_genelist_egfr <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)

# DATENSATZ VORBERITEN

TCGA_EGFR_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, PFI, PFI_months)

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$PFI_months
TCGA_EGFR_lasso$OS <- TCGA_EGFR_lasso$PFI
TCGA_EGFR_lasso$PFI <- NULL
TCGA_EGFR_lasso$PFI_months <- NULL
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#dplyr::select(Patient_ID, OS.time, OS)


fpkm_data <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% relevant_genelist_egfr$GeneID_all) %>%
  as.data.frame(.)
colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Gene_ID <- rownames(fpkm_data)
#fpkm_data$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = fpkm_data$Gene_ID, column= "SYMBOL", keytype =  "ENTREZID")
#fpkm_data <- fpkm_data %>%
#filter(!is.na(fpkm_data$Gene_Symbol))
#rownames(fpkm_data) <- fpkm_data$Gene_Symbol
#fpkm_data$Gene_ID <- NULL
#fpkm_data$Gene_Symbol <- NULL
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Patient_ID <- c(rownames(fpkm_data))
rownames(fpkm_data) <- NULL

TCGA_EGFR_lasso <- TCGA_EGFR_lasso[match(rownames(fpkm_data),TCGA_EGFR_lasso$Patient_ID),]
TCGA_EGFR_lasso <- left_join(TCGA_EGFR_lasso, fpkm_data, by = "Patient_ID", copy = TRUE) #%>%
# filter(TCGA_EGFR_lasso$Patient_ID %in% fpkm_data$Patient_ID)
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso <- as.data.frame(TCGA_EGFR_lasso)
#rownames(TCGA_EGFR_lasso) <- TCGA_EGFR_lasso$Patient_ID
#TCGA_EGFR_lasso$Patient_ID <- NULL

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/12
#TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/365

TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso$OS <- as.numeric(TCGA_EGFR_lasso$OS)

x <- as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)])
y <- data.matrix(Surv(TCGA_EGFR_lasso$OS.time, TCGA_EGFR_lasso$OS))


fit <- glmnet(x, y, family = "cox")
plot(fit,  xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
#write.table

####TEST#####
test_risk <- TCGA_EGFR_lasso %>%
  dplyr::select(Patient_ID, OS.time, OS, '7535', '51316', '147372', '105373853', '3975')
test_risk$riskscore <- test_risk$'7535' * (-0.0152467008372717) + test_risk$'51316' * -0.000711082632583063 + test_risk$'147372' * 0.0743697053744941 + test_risk$'105373853' * 0.00860660545415523 + test_risk$'3975' * 0.0206256456812425

#0.53818275*-0.00893014079168287 + 0.114368264 * -0.0245566035731716 + 24.0633075 * 0.00682598618190403 + 2.355539341 * 0.0527359814359829 + 0.054170154 * -0.484764320166299
######
riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
standard_data <- c("Patient_ID", "OS.time", "OS",lassoGene)
final_riskscore <- cbind(TCGA_EGFR_lasso[,outCol],riskScore=as.vector(riskScore))


#######

cvcoxlasso <- cv.glmnet(
  x = as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = y,
  family = 'cox')

plot(cvcoxlasso)
cvcoxlasso$lambda.1se

coef(cvcoxlasso)

plot(survival::survfit(cvcoxlasso, s = cvcoxlasso$lambda.1se, x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]), y = b))

coxlasso <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox')

plot(coxlasso)

lasso_fit <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox',
  alpha = 1,
  lambda = cv.coxlasso$lambda.1se
)
plot(lasso_fit)


#######
lasso_result <- test_risk
rownames(lasso_result) <- lasso_result$Patient_ID
lasso_result$Patient_ID <- NULL

stat <- maxstat.test(Surv(lasso_result$OS.time,lasso_result$OS)~lasso_result$riskscore, data = lasso_result, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(lasso_result$riskscore > cutoff,"high","low"))
write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
            file="lasso-Risk-Cox-best-cut-pfi_5y.txt",
            sep="\t",
            quote=F,
            row.names=F)


##### PFI_5y####
#Daten einlesen
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
relevant_genelist_egfr <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 4)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)

# DATENSATZ VORBERITEN

TCGA_EGFR_lasso <- TCGA_clinicaldata %>%
  dplyr::select(Patient_ID, PFI_5y_event, PFI_5y_months)

TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$PFI_5y_months
TCGA_EGFR_lasso$OS <- TCGA_EGFR_lasso$PFI_5y_event
TCGA_EGFR_lasso$PFI_5y_event <- NULL
TCGA_EGFR_lasso$PFI_5y_months <- NULL
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#dplyr::select(Patient_ID, OS.time, OS)


fpkm_data <- rnaseq_data %>% 
  filter(rnaseq_data$Entrez_Gene_ID %in% relevant_genelist_egfr$GeneID_all) %>%
  as.data.frame(.)
colnames(fpkm_data) <- gsub("-01","",colnames(fpkm_data), fixed=TRUE)
rownames(fpkm_data) <- fpkm_data$Entrez_Gene_ID 
fpkm_data$Entrez_Gene_ID <- NULL
fpkm_data <- t(fpkm_data)
fpkm_data <- as.data.frame(fpkm_data)
fpkm_data$Patient_ID <- c(rownames(fpkm_data))
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Gene_ID <- rownames(fpkm_data)
#fpkm_data$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = fpkm_data$Gene_ID, column= "SYMBOL", keytype =  "ENTREZID")
#fpkm_data <- fpkm_data %>%
#filter(!is.na(fpkm_data$Gene_Symbol))
#rownames(fpkm_data) <- fpkm_data$Gene_Symbol
#fpkm_data$Gene_ID <- NULL
#fpkm_data$Gene_Symbol <- NULL
#fpkm_data <- t(fpkm_data)
#fpkm_data <- as.data.frame(fpkm_data)
#fpkm_data$Patient_ID <- c(rownames(fpkm_data))
rownames(fpkm_data) <- NULL

TCGA_EGFR_lasso <- TCGA_EGFR_lasso[match(rownames(fpkm_data),TCGA_EGFR_lasso$Patient_ID),]
TCGA_EGFR_lasso <- left_join(TCGA_EGFR_lasso, fpkm_data, by = "Patient_ID", copy = TRUE) #%>%
# filter(TCGA_EGFR_lasso$Patient_ID %in% fpkm_data$Patient_ID)
#TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
#filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso <- as.data.frame(TCGA_EGFR_lasso)
#rownames(TCGA_EGFR_lasso) <- TCGA_EGFR_lasso$Patient_ID
#TCGA_EGFR_lasso$Patient_ID <- NULL

#TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/12
#TCGA_EGFR_lasso$OS.time <- TCGA_EGFR_lasso$OS.time/365

TCGA_EGFR_lasso <- TCGA_EGFR_lasso %>%
  filter(!is.na(TCGA_EGFR_lasso$OS.time))
TCGA_EGFR_lasso$OS <- as.numeric(TCGA_EGFR_lasso$OS)

x <- as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)])
y <- data.matrix(Surv(TCGA_EGFR_lasso$OS.time, TCGA_EGFR_lasso$OS))


fit <- glmnet(x, y, family = "cox")
plot(fit,  xvar = "lambda", label = FALSE, lw=2)

cvfit <- cv.glmnet(x, y, family="cox")
abline(v=log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty="dashed")
plot(cvfit)

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]
geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)
#write.table

####TEST#####
test_risk <- TCGA_EGFR_lasso %>%
  dplyr::select(Patient_ID, OS.time, OS, '7535', '51316', '147372', '105373853', '3975')
test_risk$riskscore <- test_risk$'7535' * (-0.0152467008372717) + test_risk$'51316' * -0.000711082632583063 + test_risk$'147372' * 0.0743697053744941 + test_risk$'105373853' * 0.00860660545415523 + test_risk$'3975' * 0.0206256456812425

#0.53818275*-0.00893014079168287 + 0.114368264 * -0.0245566035731716 + 24.0633075 * 0.00682598618190403 + 2.355539341 * 0.0527359814359829 + 0.054170154 * -0.484764320166299
######
riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
standard_data <- c("Patient_ID", "OS.time", "OS",lassoGene)
final_riskscore <- cbind(TCGA_EGFR_lasso[,outCol],riskScore=as.vector(riskScore))


#######

cvcoxlasso <- cv.glmnet(
  x = as.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = y,
  family = 'cox')

plot(cvcoxlasso)
cvcoxlasso$lambda.1se

coef(cvcoxlasso)

plot(survival::survfit(cvcoxlasso, s = cvcoxlasso$lambda.1se, x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]), y = b))

coxlasso <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox')

plot(coxlasso)

lasso_fit <- glmnet(
  x = data.matrix(TCGA_EGFR_lasso[,4:ncol(TCGA_EGFR_lasso)]),
  y = b,
  family = 'cox',
  alpha = 1,
  lambda = cv.coxlasso$lambda.1se
)
plot(lasso_fit)


#######
lasso_result <- test_risk
rownames(lasso_result) <- lasso_result$Patient_ID
lasso_result$Patient_ID <- NULL

stat <- maxstat.test(Surv(lasso_result$OS.time,lasso_result$OS)~lasso_result$riskscore, data = lasso_result, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(lasso_result$riskscore > cutoff,"high","low"))
write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
            file="lasso-Risk-Cox-best-cut-pfi_5y.txt",
            sep="\t",
            quote=F,
            row.names=F)

