library(readr) 
library(readxl) 
library(tidyverse) 
library(DESeq2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap) #Heatmap
library(circlize) #colours
library(RColorBrewer) # colours

contr.grp <- read_delim("~/Desktop/Doktorarbeit/Data/HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Normal.cct.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
contr.grp.all <- read_delim("~/Desktop/Doktorarbeit/Data/HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Tumor.cct.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

contr.grp.clinical <- read_delim("~/Desktop/Doktorarbeit/Data/HS_CPTAC_HNSCC_CLI.tsi.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
#die gene für den riskscore reinhauen
list <- c(7535, 1948, 494514, 51806, 147372, 645432, 56936)
list <- as.data.frame(list)
list$ID <- list$list
list$ID <- as.character(list$ID)
list$Gene_Symbol <- mapIds(org.Hs.eg.db, keys = list$ID, column= "SYMBOL", keytype =  "ENTREZID")
list$list <- NULL



## die daten vorbeiretne
contr.grp.all <- as.data.frame(contr.grp.all)
rownames(contr.grp.all) <- contr.grp.all$Idx
contr.grp.all$Idx <- NULL
contr.grp.all <- t(contr.grp.all)
contr.grp.all <- as.data.frame(contr.grp.all)
contr.grp.all <- contr.grp.all %>%
  dplyr::select(ZAP70, EFNB2, TYMSOS, CALML5, CCBE1, ARRDC5, CCDC177)

contr.grp.clinical <- contr.grp.clinical %>%
  dplyr::rename(`Patient_ID` = `case_id`)
contr.grp.clinical <- contr.grp.clinical %>%
  dplyr::select(Patient_ID, overall_survival, overall_free_status)
## riskscore berechnen
contr.grp.all$riskscore <- contr.grp.all$ZAP70 * (-0.0582152196382938) + contr.grp.all$TYMSOS * 0.013184642798603 + contr.grp.all$EFNB2 * 0.00470780696956776 + contr.grp.all$CCBE1 * 0.11219886435049 + contr.grp.all$ARRDC5 * -0.267800991188706 + contr.grp.all$CALML5 * -7.50310759764005e-05 + contr.grp.all$CCDC177 * -0.0844511106749959
##klinische daten hinzufügen
contr.grp.all$Patient_ID <- rownames(contr.grp.all)
contr.grp.all <- left_join(contr.grp.all, contr.grp.clinical, by = "Patient_ID")
contr.grp.all$overall_survival <- as.numeric(contr.grp.all$overall_survival)
contr.grp.all$overall_free_status <- as.numeric(contr.grp.all$overall_free_status)


stat <- maxstat.test(Surv(contr.grp.all$overall_survival,contr.grp.all$overall_free_status)~contr.grp.all$riskscore, data = contr.grp.all, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)
risk <- as.vector(ifelse(contr.grp.all$riskscore > cutoff,"high","low"))
contr.grp.all$risk <- risk

#write.table(cbind(id=rownames(lasso_result),lasso_result,risk),
 #           file="lasso-Risk-Cox-best-cut-pfi_5y.txt",
  #          sep="\t",
    #        quote=F,
     #       row.names=F)

### kaplan meier mit den beiden gruppen
library(survival)
library(readxl) 


contr.grp.all$overall_survival <- contr.grp.all$overall_survival/365*12
surv_obj <- Surv(contr.grp.all$overall_survival,contr.grp.all$overall_free_status)

fit <- survfit(surv_obj~risk, data = contr.grp.all) 
summary(fit)

abline(h=0.5, col="red")

plot(fit, 
     conf.int=FALSE,
     xlab = "Survival in days",ylab = "Survival Rate",
     yscale = 100,
     las=1,
     lwd=2,
     col =c("red1","blue","green","purple"))

legend("bottomleft",                                      #man kann auch Koordinaten angeben wie (18, 0.95) oder "bottomleft"
       title= "RiskScore",
       legend=c("high","low"),
       col = c("red1","blue"),  #man kann auch die Linien zeichnen lassen mit fill/col=c("red1", "blue"...)
       lty=1,
       lwd=2,
       bty = "",
       cex=1)

###neues kaplan skript
ggsurvplot(fit, 
           title = "CPTAC HNSCC",
           xlab = "Survival in months",
           xlim = c(0, 50),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,        # show deaths
           risk.table = T,      # include number at risk table
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue"),
           conf.int = F,    # obv
           pval = T,
           legend.title = "Risk",
           legend.labs = c("high", "low"))

survdiff(Surv(contr.grp.all$overall_survival,contr.grp.all$overall_free_status) ~ risk ,data = contr.grp.all)  
