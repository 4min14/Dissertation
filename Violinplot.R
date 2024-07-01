library(openxlsx)       
library(ggplot2)
library(ggpubr)
library(forcats)
library(readr) 
library(readxl) 
library(easyGgplot2)


### Daten einlesen

rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
violinplot_fpkm <- rnaseq_data
violinplot_fpkm <- as.data.frame(violinplot_fpkm)
rownames(violinplot_fpkm) <- violinplot_fpkm$Entrez_Gene_ID
violinplot_fpkm$Entrez_Gene_ID <- NULL
colnames(violinplot_fpkm) <- gsub("-01","",colnames(violinplot_fpkm), fixed=TRUE)
violinplot_fpkm <- as.data.frame(violinplot_fpkm) %>%
  t(.)
violinplot_fpkm <- as.data.frame(violinplot_fpkm) %>%
  dplyr::select(`1956`)
violinplot_fpkm$Patient_ID <- rownames(violinplot_fpkm)
rownames(violinplot_fpkm) <- NULL
violinplot_fpkm$EGFR <- violinplot_fpkm$'1956'
violinplot_fpkm$'1956' <- NULL



violinplot_clusterdata <- clinical %>%
  dplyr::select(Cluster)
violinplot_clusterdata$Patient_ID <- rownames(violinplot_clusterdata)
rownames(violinplot_clusterdata) <- NULL


violinplot_data_egfr <- left_join(violinplot_clusterdata, violinplot_fpkm, by="Patient_ID") %>%
  dplyr::select(Cluster, EGFR)


ggplot(data = violinplot_data_egfr, aes(x=Cluster, fill=Cluster, y=EGFR))+geom_violin()+ylim(0,45)+geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("green","red1","purple","blue"))


anovatest <- aov(EGFR ~ Cluster, data = violinplot_data_egfr)
# Summary of the analysis
summary(anovatest)
TukeyHSD(anovatest)
###### T-Test (bei 2 gruppen)
#a=

#b=kk$Clusters
#t.test(a~b, var.equal = TRUE, alternative = "two.sided")
##### Violinplot für Progenydaten

violinplot_progeny <- progenydata
rownames(violinplot_progeny) <- violinplot_progeny$...1
violinplot_progeny <- as.data.frame(violinplot_progeny)
violinplot_progeny$...1 <- NULL
violinplot_progeny <- t(violinplot_progeny)
violinplot_progeny <- as.data.frame(violinplot_progeny)
violinplot_progeny$progenyscore <- violinplot_progeny$EGFR 
violinplot_progeny$EGFR <- NULL
violinplot_progeny$Patient_ID <- rownames(violinplot_progeny)
rownames(violinplot_progeny) <- NULL

violinplot_progeny <- left_join(violinplot_progeny, violinplot_clusterdata, by = "Patient_ID")
violinplot_progeny$Patient_ID <- NULL
violinplot_progeny <- violinplot_progeny %>%
  dplyr::filter(!is.na(violinplot_progeny$Cluster))

ggplot(data = violinplot_progeny, aes(x=Cluster, fill=Cluster, y=progenyscore))+geom_violin()+ylim(-2,2)+geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("green","red1","purple","blue"))

anovatest_prog <- aov(progenyscore ~ Cluster, data = violinplot_progeny)
# Summary of the analysis
summary(anovatest_prog)
TukeyHSD(anovatest_prog)

## Violinplot um die EGFR Verteilung in den Riskgruppen zu zeigen
# Daten einlesen
risks <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1)
rnaseq_data <- read_excel("~/Desktop/Doktorarbeit/Data/Heatmap_data_EGFR.xlsx", sheet = 1) #einfach die TCGA FPKM daten
violinplot_fpkm <- rnaseq_data
violinplot_fpkm <- as.data.frame(violinplot_fpkm)
rownames(violinplot_fpkm) <- violinplot_fpkm$Entrez_Gene_ID
violinplot_fpkm$Entrez_Gene_ID <- NULL
colnames(violinplot_fpkm) <- gsub("-01","",colnames(violinplot_fpkm), fixed=TRUE)
violinplot_fpkm <- as.data.frame(violinplot_fpkm) %>%
  t(.)
violinplot_fpkm <- as.data.frame(violinplot_fpkm) %>%
  dplyr::select(`1956`)
violinplot_fpkm$Patient_ID <- rownames(violinplot_fpkm)
rownames(violinplot_fpkm) <- NULL
violinplot_fpkm$EGFR <- violinplot_fpkm$'1956'
violinplot_fpkm$'1956' <- NULL


violinplot_risk.egfr <- left_join(risks, violinplot_fpkm, by = "Patient_ID") %>%
  dplyr::rename(`Risk` = `risk`)

# plotten

plot <- ggplot(data = violinplot_risk.egfr, aes(x=Risk, fill=Risk, y=EGFR))+geom_violin()+ylim(0,65)  +geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("black","grey"))

ggplot2.customize(plot, legendItemOrder = c("low", "high"))
# Remove plot legend
print(plot)

anovatest_prog <- aov(EGFR ~ Risk, data = violinplot_risk.egfr)
# Summary of the analysis
summary(anovatest_prog)
TukeyHSD(anovatest_prog)

## Violinplot um die EGFR PROGENy in den Riskgruppen zu zeigen
# Daten einlesen
risks <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1)
risks <- risks %>%
  dplyr::select(id, risk) %>%
  dplyr::rename(`Patient_ID` = `id`)
progenydata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA_HNSCC_progeny.xlsx", sheet = 1)
violinplot_progeny <- progenydata
rownames(violinplot_progeny) <- violinplot_progeny$...1
violinplot_progeny <- as.data.frame(violinplot_progeny)
violinplot_progeny$...1 <- NULL
violinplot_progeny <- t(violinplot_progeny)
violinplot_progeny <- as.data.frame(violinplot_progeny)
violinplot_progeny$progenyscore <- violinplot_progeny$EGFR 
violinplot_progeny$EGFR <- NULL
violinplot_progeny$Patient_ID <- rownames(violinplot_progeny)
rownames(violinplot_progeny) <- NULL


violinplot_risk.egfr <- left_join(risks, violinplot_progeny, by = "Patient_ID") %>%
  dplyr::rename(`Risk` = `risk`)
view(violinplot_risk.egfr)
# plotten

plot <- ggplot(data = violinplot_risk.egfr, aes(x=Risk, fill=Risk, y=progenyscore))+geom_violin()+ylim(-2,2)  +geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("black","grey"))

ggplot2.customize(plot, legendItemOrder = c("low", "high"))
# Remove plot legend
print(plot)

anovatest_prog <- aov(progenyscore ~ Risk, data = violinplot_risk.egfr)
# Summary of the analysis
summary(anovatest_prog)
TukeyHSD(anovatest_prog)

#####Violinplot um mutationcounts zu zeigen in riskgroups
###Daten einlesen
risks <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1)
risks$Patient_ID <- risks$id
risks <- risks %>%
  dplyr::select(Patient_ID, risk)
TCGA_clinicaldata <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1)
clincalandrisk <- left_join(TCGA_clinicaldata, risks, by = "Patient_ID")
violinplot_mutation <- clincalandrisk %>%
  dplyr::select(risk, "Mutation Count")
violinplot_mutation <- violinplot_mutation %>%
  dplyr::rename(`Risk` = `risk`)
violinplot_FGA <- clincalandrisk %>%
  dplyr::select(risk, "Fraction Genome Altered")
violinplot_FGA <- violinplot_FGA %>%
  dplyr::rename(`Risk` = `risk`)
#Nulls raus
violinplot_mutation <- violinplot_mutation %>%
  filter(!is.na(violinplot_mutation$risk))
violinplot_FGA <- violinplot_FGA %>%
  filter(!is.na(violinplot_FGA$risk))
violinplot_mutation <- violinplot_mutation %>%
  filter(!is.na(violinplot_mutation$`Mutation Count`))
violinplot_mutation$`Mutation Count` <- log2(violinplot_mutation$`Mutation Count`)
violinplot_FGA <- violinplot_FGA %>%
  filter(!is.na(violinplot_FGA$risk))
violinplot_FGA <- violinplot_FGA %>%
  filter(!is.na(violinplot_FGA$`Fraction Genome Altered`))
###für mutationcount
plot <- ggplot(data = violinplot_mutation, aes(x=Risk, fill=Risk, y=`Mutation Count`))+geom_violin()+ylim(0,12)  +geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("red1","blue"))+ylab("Mutation Count(log2)")

print(plot)
ggplot2.customize(plot, legendItemOrder = c("low", "high"))
#tests
t.test(violinplot_mutation$`Mutation Count`~violinplot_mutation$Risk)
anovatest_prog <- aov(`Mutation Count` ~ Risk, data = violinplot_mutation)
summary(anovatest_prog)
TukeyHSD(anovatest_prog)
###für fraction genome altered
plot <- ggplot(data = violinplot_FGA, aes(x=Risk, fill=Risk, y=`Fraction Genome Altered`))+geom_violin()+ylim(0,1)  +geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("red1","blue"))

print(plot)
ggplot2.customize(plot, legendItemOrder = c("low", "high"))
#tests
t.test(violinplot_FGA$`Fraction Genome Altered`~violinplot_FGA$Risk)
anovatest_prog <- aov(`Fraction Genome Altered` ~ Risk, data = violinplot_FGA)
summary(anovatest_prog)
TukeyHSD(anovatest_prog)

##Violinplot mit den Methylierungsdaten

meth_beta_mean <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA_HNSC_Methylation_mean.xlsx", sheet = 1)
rownames(meth_beta_mean) <-  meth_beta_mean$ID
meth_beta_mean <- as.data.frame(meth_beta_mean)
meth_beta_mean <- t(meth_beta_mean)
colnames(meth_beta_mean) <- gsub("-01","",colnames(meth_beta_mean), fixed=TRUE)
meth_beta_mean$Patient_ID <- rownames(meth_beta_mean)
meth_beta_mean <- meth_beta_mean %>%
  dplyr::select(Patient_ID, beta_mean)
risk_HNSC <- HNSC_risk %>%
  dplyr::rename(Patient_ID = bcr_patient_barcode) %>%
  dplyr::select(Patient_ID, risk)
violinplot_methbetamean <- left_join(meth_beta_mean, risk_HNSC, by = "Patient_ID")
violinplot_methbetamean <- violinplot_methbetamean %>%
  dplyr::rename(`Risk` = `risk`)
violinplot_methbetamean <- violinplot_methbetamean %>%
  filter(!is.na(violinplot_methbetamean$Risk))
violinplot_methbetamean <- violinplot_methbetamean %>%
  filter(!is.na(violinplot_methbetamean$beta_mean))
violinplot_methbetamean <- as.data.frame(violinplot_methbetamean) %>%
  dplyr::select(Risk, beta_mean)
violinplot_methbetamean$beta_mean <- as.numeric(violinplot_methbetamean$beta_mean)
###ploten
plot <- ggplot(data = violinplot_methbetamean, aes(x=Risk, fill=Risk, y=beta_mean))+geom_violin()+ylim(0,1)  +geom_boxplot(width=0.3)+ theme_classic()+scale_fill_manual(values = c("red1","blue"))+ylab("Global methylation (beta-mean)"                             )

print(plot)
ggplot2.customize(plot, legendItemOrder = c("low", "high"))
#tests
t.test(violinplot_methbetamean$beta_mean~violinplot_methbetamean$Risk)
anovatest_meth <- aov(beta_mean ~ Risk, data = violinplot_methbetamean)
summary(anovatest_meth)

write.table(riskybisky, file = "HNSC_risk.txt", quote = F, sep = " ")
TukeyHSD(anovatest_meth)