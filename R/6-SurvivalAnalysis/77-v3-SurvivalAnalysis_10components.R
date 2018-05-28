#version 3 (Apr/3/2018)
#NOT INCLUDED IN THE PAPER go to 79

#analize the data
rm(list = ls())
source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
#load("temp/77-v3-data.RData")
load("temp/77-v3-result.RData")
#required for the cox proportional hazard model
loadls("plyr survival missForest survAUC prodlim survminer")

loadData <- function()
{
  load("results/10March2018/77-v2-data.RData")
  w_muta = w

  load("temp/w_clinical.RData")
  w <- as.data.frame(cbind(w,row.names(w)))
  names(w) <- c( paste("C",1:10,sep=''), "Patient.ID")
  w_clini=w
  rm(w)

  survivalClinicalNMF_clini =  merge(survivalClinical, w_clini, by="Patient.ID")
  which(is.na(survivalClinicalNMF_clini))
  count(survivalClinicalNMF_clini, 'Overall.Survival.Status')

  save.image("temp/77-v3-data.RData")
}

#' Survival Analysis
survAna <- function(){

  #clinical variables
  coxFit <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
                    Diagnosis.Age+#ER.Status+HER2.Status+
                    Metastasis.Coded+#PR.Status+
                    Tumor..T1.Coded,
                  x=TRUE, y=TRUE, method = "breslow", data=survivalClinicalNTF)
  summary(coxFit)
  #Likelihood ratio test= 19.09  on 5 df,   p=0.001851
  #Wald test            = 25.37  on 5 df,   p=0.0001184
  #Score (logrank) test = 30.29  on 5 df,   p=1.292e-05

  #NFT components
  coxFit2 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ Comp.1+
                     Comp.2+Comp.3+Comp.4+Comp.5+Comp.6+Comp.7+Comp.8+Comp.9+Comp.10,
                   method = "breslow", data=survivalClinicalNTF)
  summary(coxFit2)
  #p=8.327e-15


  #clinical variables and components
  coxFit3 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
                     Diagnosis.Age+#ER.Status+
                     HER2.Status+Metastasis.Coded+#PR.Status+
                     Tumor..T1.Coded+
                     Comp.1+
                     Comp.2+Comp.3+Comp.4+Comp.5+Comp.6+Comp.7+Comp.8+Comp.9+Comp.10,
                   method = "breslow", data=survivalClinicalNTF)
  summary(coxFit3)
  #p=0



  #NMT components
  # coxFit4 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~C1 +
  #                    C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10,
  #                  method = "breslow", data=survivalClinicalNMF)
  # summary(coxFit4)
  #p=0 but last for ever because the elements are strings

  patientVsComponents <- as.data.frame(w_muta[,1:10])
  patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
  patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)

  patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                  C=patientVsComponents$whichmax)
  survivalClinicalNMF_2 = merge(x=survivalClinical, y=patientAfiliation,
                                by.x="Patient.ID", by.y="patient")

  coxFit4_v2 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~factor(C),
                   method = "breslow", data=survivalClinicalNMF_2)
  summary(coxFit4_v2)
  #Score (logrank) test = 8.41  on 10 df,   p=0.5886
  #Score (logrank) test = 10.45  on 9 df,   p=0.315




  #clinical variables and NMF components
  coxFit5 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~factor(C)+
                     Diagnosis.Age+#ER.Status+
                     HER2.Status+Metastasis.Coded+#PR.Status+
                     Tumor..T1.Coded,
                      method = "breslow", data=survivalClinicalNMF_2)
  summary(coxFit5)
  #Score (logrank) test = 42.37  on 18 df,   p=0.0009802


  #clinical variables and NMF components
  coxFit6 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     Diagnosis.Age+#ER.Status+HER2.Status+
                     Metastasis.Coded+#PR.Status+
                     Tumor..T1.Coded + as.numeric(C1) +
                     as.numeric(C1) + as.numeric(C2) + as.numeric(C3) +
                     as.numeric(C4) + as.numeric(C5) + as.numeric(C6) +
                     as.numeric(C7) + as.numeric(C7) + as.numeric(C9) + as.numeric(C10),
                   method = "breslow", data=survivalClinicalNMF_clini)
  summary(coxFit6)
  # Likelihood ratio test= 82.58  on 14 df,   p=9.36e-12
  # Wald test            = 74.67  on 14 df,   p=2.718e-10
  # Score (logrank) test = 95.17  on 14 df,   p=3.975e-14


  patientVsComponentsC <- as.data.frame(w_clini[,1:10])
  patientVsComponentsC[, "max"] <- apply(patientVsComponentsC, 1, max)
  patientVsComponentsC[, "whichmax"] <- apply(patientVsComponentsC, 1, which.max)

  patientAfiliationC <- data.frame(patient=row.names(patientVsComponentsC),
                                  C=patientVsComponentsC$whichmax)
  survivalClinicalNMF_clini_2 = merge(x=survivalClinical, y=patientAfiliationC,
                                by.x="Patient.ID", by.y="patient")

  #NMF CLINI with clinical
  coxFit6 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     factor(C)+
                     Diagnosis.Age+#ER.Status+HER2.Status+
                     Metastasis.Coded+#PR.Status+
                     Tumor..T1.Coded ,
                   method = "breslow", data=survivalClinicalNMF_clini_2)
  summary(coxFit6)
  # Likelihood ratio test= 35  on 14 df,   p=0.001472
  # Wald test            = 44.22  on 14 df,   p=5.448e-05
  # Score (logrank) test = 55.52  on 14 df,   p=7.038e-07

  #NMF CLINI without clinical
  coxFit7 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     factor(C),
                   method = "breslow", data=survivalClinicalNMF_clini_2)
  summary(coxFit7)
  # Likelihood ratio test= 17.72  on 9 df,   p=0.03854
  # Wald test            = 21.05  on 9 df,   p=0.01245
  # Score (logrank) test = 28.64  on 9 df,   p=0.0007445

  save.image("temp/77-v3-result.RData")
}


#PLOTS
plots <- function(){

  #Kaplan Meier population
  plot(survfit(coxFit), main = "clinical")
  plot(survfit(coxFit2), main = "NTF")
  plot(survfit(coxFit3), main = "NTF + clinical")
  plot(survfit(coxFit4_v2), main = "NMF")
  plot(survfit(coxFit5), main = "NMF + clinical")

  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
                    Diagnosis.Age+#ER.Status+
                    HER2.Status+Metastasis.Coded+#PR.Status+
                    Tumor..T1.Coded,
                 method = "breslow",
                 bandwidth = 100,
                 data=survivalClinicalNTF)
  #NMF
  km4 <- prodlim(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                   factor(component),
                   data=survivalClinicalNMF_2)
  plot(km4, logrank=TRUE)


  #Complex Plots
  #clinical variables
  mfit <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
                    Diagnosis.Age+#ER.Status+
                    HER2.Status+Metastasis.Coded+#PR.Status+
                    Tumor..T1.Coded,
                  data=survivalClinicalNTF)
  plot.SurvivalComplex (coxFit, mfit, "clinical variables")

  ####////////////
  #Using survminer library
  #Clinical
  ggsurvplot(mfit, conf.int = F, legend="none")
  ggrisktable(mfit, tables.theme = theme_cleantable(), color=NULL,
              fontsize = 1)

  #NFT components
  patientVsComponents2 <- as.data.frame(A[,1:10])
  patientVsComponents2[, "max"] <- apply(patientVsComponents2, 1, max)
  patientVsComponents2[, "whichmax"] <- apply(patientVsComponents2, 1, which.max)

  patientAfiliation2 <- data.frame(patient=row.names(patientVsComponents2),
                                  component=patientVsComponents2$whichmax)
  survivalClinicalNTF_2 = merge(x=survivalClinical, y=patientAfiliation2,
                                by.x="Patient.ID", by.y="patient")
  mfit2 <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     factor(component),
                   data=survivalClinicalNTF_2)
  ggsurvplot(mfit2, conf.int = F)
  #these plots look too messy

  #NMF
  mfit4 <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     factor(C),
                    data=survivalClinicalNMF_2)
  ggsurvplot(mfit4, conf.int = F)
  ggrisktable(mfit4, tables.theme = theme_cleantable(), color=NULL)
  #ZOOM and copy to word


  #' These plots look too messy so we should get the top three subtypes and plot those
  #' > table(survivalClinicalNMF_2$C)
  #
  # 1*  2  3  4 *5  6  7 *8  9 10
  # 68 50 41 53 54 53 50 58 38 42
  # > table(survivalClinicalNTF_2$component)
  #
  # 1   2   3   4   5   6   7   8   9  10
  # 9 138  34  41   7  14 109  99  37  19
  #'
  #' The component with largest population are 2,7,8
  #' > dim(survivalClinicalNTF_2red)
  # [1] 346  17
  #NTF top3
  survivalClinicalNTF_top3 <- survivalClinicalNTF_2[which(survivalClinicalNTF_2$component %in% c(2,7,8)),]
  mfit2_top3 <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                     factor(component),
                   data=survivalClinicalNTF_top3)
  ggsurvplot(mfit2_top3, conf.int = F)

  coxFit_top3 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                           factor(component),
                         data=survivalClinicalNTF_top3)
  summary(coxFit_top3)
  #these plots look too messy

  ###
  #NMF top 3
  survivalClinicalNMF_top3 <- survivalClinicalNMF_2[which(survivalClinicalNMF_2$C %in% c(1,5,8)),]
  nmffit_top3 <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                          factor(C),
                        data=survivalClinicalNMF_top3)
  ggsurvplot(nmffit_top3, conf.int = F)

  coxFitNMF_top3 <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                         factor(C),
                       data=survivalClinicalNMF_top3)
  summary(coxFitNMF_top3)


}

save.image("temp/77-v3.RData")

