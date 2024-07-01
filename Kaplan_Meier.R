library(survival)
library(readxl) 


survival_sorv <- read_excel("~/Desktop/Doktorarbeit/Data/TCGA-HNSC_Master_clinical data_v1.xlsx", sheet = 1) 
survival_sorv <- left_join(survival_sorv, cluster, by = "Patient_ID")
survival_sorv <- survival_sorv %>%
  filter(!is.na(survival_sorv$Cluster)) %>%
 # filter(!is.na(survival_sorv$OS_5y_months)) %>%
  filter(!OS_5y_event == "z")
#survival_sorv <- survival_sorv %>%
#  filter(!is.na(survival_sorv$OS_5y_months))

survival_sorv$OS_5y_event <- as.numeric(survival_sorv$OS_5y_event)
surv_obj <- Surv(survival_sorv$OS_5y_months,survival_sorv$OS_5y_event)

fit_1 <- survfit(surv_obj~1, data = survival_sorv) 
summary(fit_1)

plot(fit_1, 
     conf.int=FALSE,
     xlab = "Survival in Days",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

fit <- survfit( Surv(OS_5y_months, OS_5y_event) ~ Cluster, data = survival_sorv )
ggsurvplot(fit, 
           title = "HNSC",
           xlab = "Survival in months",
           xlim = c(0, 60),
           break.x.by = 12,
           ylab = "Survival Rate",
           mark.time = T,        # show deaths
           risk.table = T,      # include number at risk table
           risk.table.title = "",
           risk.table.height = .25,
           palette =c("red1","blue", "green", "purple"),
           conf.int = F,    # obv
           pval = T,
           legend.title = "Cluster",
           legend.labs = c("B1", "B2b", "A1", "B2a"))

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "OS_5y_months",ylab = "Survival Rate",
     yscale = 100,
     las=1,
     lwd=2,
     col =c("red1","blue","green","purple"))

legend("bottomleft",                                      #man kann auch Koordinaten angeben wie (18, 0.95) oder "bottomleft"
       title= "Clusters",
       legend=c("B1","B2b","A1","B2a"),
       col = c("red1","blue","green","purple"),  #man kann auch die Linien zeichnen lassen mit fill/col=c("red1", "blue"...)
       lty=1,
       lwd=2,
       bty = "",
       cex=1)

survdiff(Surv(OS_months, OS) ~ Cluster,data = survival_sorv)  






library(survival)
library(readxl) 


kaplan_riskbased <- read_excel("~/Desktop/Doktorarbeit/Data/LassoRiskscore_fpkm.xlsx", sheet = 1) 
kaplan_riskbased$OS.time <- kaplan_riskbased$OS.time*12


surv_obj <- Surv(kaplan_riskbased$OS.time,kaplan_riskbased$OS)

fit_1 <- survfit(surv_obj~1, data = kaplan_riskbased)
summary(fit_1)



plot(fit_1, 
     conf.int=FALSE,
     xlab = "Survival in Months",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

abc <- survfit( Surv(OS.time, OS) ~ risk, data = kaplan_riskbased )

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "OS_5y_months",ylab = "Survival Rate",
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

survdiff(Surv(OS.time, OS) ~ risk ,data = kaplan_riskbased)  




library(survival)
library(readxl) 


kaplan_riskbased <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss.xlsx", sheet = 1) 
kaplan_riskbased$OS.time <- kaplan_riskbased$OS.time*30


surv_obj <- Surv(kaplan_riskbased$OS.time,kaplan_riskbased$OS)

fit_2 <- survfit(surv_obj~1, data = kaplan_riskbased) 
summary(fit_2)


plot(fit_2, 
     conf.int=FALSE,
     xlab = "Survival in days",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

abc <- survfit( Surv(OS.time, OS) ~ risk, data = kaplan_riskbased )

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "DSS_days",ylab = "Survival Rate",
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

survdiff(Surv(OS.time, OS) ~ risk ,data = kaplan_riskbased)  


library(survival)
library(readxl) 


kaplan_riskbased <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_dss_5y.xlsx", sheet = 1) 
#kaplan_riskbased$OS.time <- kaplan_riskbased$OS.time*30


surv_obj <- Surv(kaplan_riskbased$OS.time,kaplan_riskbased$OS)

fit_25 <- survfit(surv_obj~1, data = kaplan_riskbased) 
summary(fit_25)


plot(fit_25, 
     conf.int=FALSE,
     xlab = "Survival in days",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

abc <- survfit( Surv(OS.time, OS) ~ risk, data = kaplan_riskbased )

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "DSS_months",ylab = "Survival Rate",
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

survdiff(Surv(OS.time, OS) ~ risk ,data = kaplan_riskbased)  


library(survival)
library(readxl) 


kaplan_riskbased <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_pfi.xlsx", sheet = 1) 
#kaplan_riskbased$OS.time <- kaplan_riskbased$OS.time*30


surv_obj <- Surv(kaplan_riskbased$OS.time,kaplan_riskbased$OS)

fit_3 <- survfit(surv_obj~1, data = kaplan_riskbased) 
summary(fit_3)


plot(fit_3, 
     conf.int=FALSE,
     xlab = "Survival in days",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

abc <- survfit( Surv(OS.time, OS) ~ risk, data = kaplan_riskbased )

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "PFI_months",ylab = "Survival Rate",
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

survdiff(Surv(OS.time, OS) ~ risk ,data = kaplan_riskbased)  

library(survival)
library(readxl) 


kaplan_riskbased <- read_excel("~/Desktop/Doktorarbeit/Data/lassocox_egfr_pfi_5y.xlsx", sheet = 1) 
#kaplan_riskbased$OS.time <- kaplan_riskbased$OS.time*30


surv_obj <- Surv(kaplan_riskbased$OS.time,kaplan_riskbased$OS)

fit_35 <- survfit(surv_obj~1, data = kaplan_riskbased) 
summary(fit_35)


plot(fit_35, 
     conf.int=FALSE,
     xlab = "Survival in days",ylab = "Survival Rate",
     yscale = 100,  
     las=1,
     lwd=2)

abc <- survfit( Surv(OS.time, OS) ~ risk, data = kaplan_riskbased )

abline(h=0.5, col="red")

plot(abc, 
     conf.int=FALSE,
     xlab = "PFI_months",ylab = "Survival Rate",
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

survdiff(Surv(OS.time, OS) ~ risk ,data = kaplan_riskbased)  


