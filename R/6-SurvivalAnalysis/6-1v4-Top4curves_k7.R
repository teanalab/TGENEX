#version 3 (May/28/2018)
#' Select the top 4 groups with best p-value
#' Pre: 79-Top3curvesFrom77_v2

rm(list = ls())

LoadMyData <- function()
{
  load("data/factors/factors_7.RData")
  load("data/NMF/NMF_PxC_R_7.RData")
  load("data/NMF/NMF_PxM_R_7.RData")
  load("data/survivalClinical-5-1_30.RData")
  kMax <- 7
  load("data/tensorClinical.RData")
  row.names(tensorClinical) <- tensorClinical$patient.bcr_patient_barcode
  tensorClinical <- tensorClinical[patients,-1]
  load("data/myLib.RData")
  source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
  load("data/patientAfiliation/patientAfiliation_7.RData")
  save.image(file="temp/6-1v4-data.RData")
}
load("temp/6-1v4-data.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)

# k = 7
# table(patientAfiliation$component)
# 1   2   3   4   5   6
# 4  10   4 302   7 129

# k = 12
# 1   2   3   4   5   6   7   9  10  11  12
# 283  13   3  13  32  21   7  46  13  20   5


survivalAnalysisNTF <- function(patiAfi_NTF){
  ### by score
  patiFact_NTF <- patientsF
  patiFact_NTF <- merge(as.data.frame(patiFact_NTF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  row.names(patiFact_NTF) <- patiFact_NTF$Row.names
  patiFact_NTF <- patiFact_NTF[,-1]
  #NFT components
  coxFit_fact_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiFact_NTF)
  summary(coxFit_fact_NTF)
  # Concordance= 0.669  (se = 0.054 )
  # Rsquare= 0.022   (max possible= 0.606 )
  # Likelihood ratio test= 10.11  on 7 df,   p=0.2
  # Wald test            = 14.17  on 7 df,   p=0.05
  # Score (logrank) test = 20.48  on 7 df,   p=0.005


  ### by principal verotype
  patiAfi_NTF <- patientAfiliation
  row.names(patiAfi_NTF) <- patiAfi_NTF$patient
  patiAfi_NTF <- merge(as.data.frame(patiAfi_NTF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiAfi_NTF <- patiAfi_NTF[,-1]
  row.names(patiAfi_NTF) <- patiAfi_NTF$patient
  patiAfi_NTF <- patiAfi_NTF[,-1]
  #NFT components
  coxFit_Afi_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~factor(component),
                      x=T, y=T, model=T, method = "breslow",
                      data=patiAfi_NTF)
  summary(coxFit_Afi_NTF)
  # Concordance= 0.59  (se = 0.044 )
  # Rsquare= 0.013   (max possible= 0.606 )
  # Likelihood ratio test= 6.08  on 5 df,   p=0.3
  # Wald test            = 7.22  on 5 df,   p=0.2
  # Score (logrank) test = 8.79  on 5 df,   p=0.1


  ## by both
  patiBoth_NTF <- merge(as.data.frame(patientsF), as.data.frame(patiAfi_NTF), by='row.names', all=TRUE)
  row.names(patiBoth_NTF) <- patiBoth_NTF$Row.names
  patiBoth_NTF <- patiBoth_NTF[,-1]
  patiBoth_NTF$component <- factor(patiBoth_NTF$component)
  #NFT components
  coxFit_both_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiBoth_NTF)
  summary(coxFit_both_NTF)
  # Concordance= 0.69  (se = 0.054 )
  # Rsquare= 0.035   (max possible= 0.606 )
  # Likelihood ratio test= 16.09  on 12 df,   p=0.2
  # Wald test            = 20.36  on 12 df,   p=0.06
  # Score (logrank) test = 29.85  on 12 df,   p=0.003

save.image("temp/6-1v4-survNTF.RData")
}
#load("temp/6-1v4-survNTF.RData")

survivalAnalysisNMF_Mut <- function(patiAfi_NTF){
  load("data/patientA_NMF_M/patientA_NMF_M_7.RData")
  ### by score
  #NMF mutation
  #normalize rows
  pati_NMF_Mut <- data.frame( t(apply(patientNMF_PxM_R_7,1,function(x){x/sum(x)})) )
  names(pati_NMF_Mut) <- paste('V',seq(1:kMax),sep='')
  pati_NMF_Mut <- merge(as.data.frame(pati_NMF_Mut), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  pati_NMF_Mut <- pati_NMF_Mut[,-1]
  coxFit_NMF_mut <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~.,
                          x=T, y=T, model=T, method = "breslow",
                          data=pati_NMF_Mut)
  summary(coxFit_NMF_mut)
  # Concordance= 0.615  (se = 0.054 )
  # Rsquare= 0.022   (max possible= 0.606 )
  # Likelihood ratio test= 9.98  on 6 df,   p=0.1
  # Wald test            = 11.43  on 6 df,   p=0.08
  # Score (logrank) test = 11.78  on 6 df,   p=0.07

  ### by principal verotype
  patiAfi_NMF_mut <- patientAfiliation
  row.names(patiAfi_NMF_mut) <- patiAfi_NMF_mut$patient
  patiAfi_NMF_mut <- merge(as.data.frame(patiAfi_NMF_mut), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiAfi_NMF_mut <- patiAfi_NMF_mut[,-1]
  row.names(patiAfi_NMF_mut) <- patiAfi_NMF_mut$patient
  patiAfi_NMF_mut <- patiAfi_NMF_mut[,-1]
  coxFit_Afi_NMF_mut <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~factor(component),
                          x=T, y=T, model=T, method = "breslow",
                          data=patiAfi_NMF_mut)
  summary(coxFit_Afi_NMF_mut)
  # Concordance= 0.632  (se = 0.052 )
  # Rsquare= 0.022   (max possible= 0.606 )
  # Likelihood ratio test= 10.06  on 6 df,   p=0.1
  # Wald test            = 10.79  on 6 df,   p=0.1
  # Score (logrank) test = 11.84  on 6 df,   p=0.07

  save.image("temp/6-1v4-survNMF_mut.RData")
}
#load("temp/6-1v4-survNMF_mut.RData")
# table(patientAfiliation$component)
# 1   2   3   4   5   6   7
# 31  85   9 137  26   9 159

survivalAnalysisNMF_Clini <- function(patiAfi_NTF){
  load("data/patientA_NMF_C/patientA_NMF_C_7.RData")
  ### by score

  #NMF clinical
  #normalize rows
  pati_NMF_Clin <- data.frame( t(apply(patientNMF_PxC_R_7,1,function(x){x/sum(x)})) )
  names(pati_NMF_Clin) <- paste('V',seq(1:kMax),sep='')
  pati_NMF_Clin <- merge(as.data.frame(pati_NMF_Clin), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  pati_NMF_Clin <- pati_NMF_Clin[,-1]

  coxFit_NMF_cli <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~.,
                          x=T, y=T, model=T, method = "breslow",
                          data=pati_NMF_Clin)
  summary(coxFit_NMF_cli)
  # Concordance= 0.678  (se = 0.054 )
  # Rsquare= 0.02   (max possible= 0.606 )
  # Likelihood ratio test= 9.33  on 6 df,   p=0.2
  # Wald test            = 9.55  on 6 df,   p=0.1
  # Score (logrank) test = 9.96  on 6 df,   p=0.1



  ### by principal verotype
  patiAfi_NMF_clin <- patientAfiliation
  row.names(patiAfi_NMF_clin) <- patiAfi_NMF_clin$patient
  patiAfi_NMF_clin <- merge(as.data.frame(patiAfi_NMF_clin), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiAfi_NMF_clin <- patiAfi_NMF_clin[,-1]
  row.names(patiAfi_NMF_clin) <- patiAfi_NMF_clin$patient
  patiAfi_NMF_clin <- patiAfi_NMF_clin[,-1]
  coxFit_Afi_NMF_clin <- coxph(Surv(time = Overall.Survival..Months.,
                                   event = Overall.Survival.Status)~factor(component),
                              x=T, y=T, model=T, method = "breslow",
                              data=patiAfi_NMF_clin)
  summary(coxFit_Afi_NMF_clin)
  # Concordance= 0.686  (se = 0.052 )
  # Rsquare= 0.026   (max possible= 0.606 )
  # Likelihood ratio test= 12.16  on 6 df,   p=0.06
  # Wald test            = 11.51  on 6 df,   p=0.07
  # Score (logrank) test = 12.73  on 6 df,   p=0.05

  save.image("temp/6-1v4-survNMF_clin.RData")
}
#load("temp/6-1v4-survNMF_clin.RData")
# table(patientAfiliation$component)
# 1   2   3   4   5   6   7
# 49 152  47  32  59  17 100

#' Survival Analysis
survAna <- function(){
  #row clinical variables
  tensorClini <- merge(as.data.frame(tensorClinical), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  tensorClini <- tensorClini[,-1]
  coxFit_tensorClini <- coxph(Surv(time = Overall.Survival..Months.,
                       event = Overall.Survival.Status)~.,
                  x=T, y=T, model=T, method = "breslow",
                  data=tensorClini)
  summary(coxFit_tensorClini)
  ###
  # Concordance= 0.959  (se = 0.058 )
  # Rsquare= 0.378   (max possible= 0.655 )
  # Likelihood ratio test= 151.1  on 97 df,   p=4e-04
  # Wald test            = 273.2  on 97 df,   p=<2e-16
  # Score (logrank) test = 187.9  on 97 df,   p=9e-08

  #WARNING: coxph ran out of iterations and did not converge because:
  with(tensorClini, table(Overall.Survival.Status, Converted.Stage))
  #solution:
  #group all stage III together
  #also group age
  #https://stats.stackexchange.com/questions/66591/coxph-ran-out-of-iterations-and-did-not-converge




  #random
  coxFit_random <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~.,
                          x=T, y=T, model=T, method = "breslow",
                          data=randPatfm)
  summary(coxFit_random)
  # Concordance= 0.638  (se = 0.054 )
  # Rsquare= 0.053   (max possible= 0.606 )
  # Likelihood ratio test= 25.03  on 29 df,   p=0.7
  # Wald test            = 24.08  on 29 df,   p=0.7
  # Score (logrank) test = 24.36  on 29 df,   p=0.7

  save.image("temp/6-1v4-surAna.RData")
}
#load("temp/6-1v4-surAna.RData")


###
plotPAM504paper <- function(){
  #PAM50
  row.names(survClinical) <- survClinical$patient.bcr_patient_barcode
  PAM50_surv_c <- survClinical[patients,c("patient.bcr_patient_barcode", "PAM50.Subtype")]
  PAM50_surv_c <- merge(as.data.frame(PAM50_surv_c), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  PAM50_surv_c <- PAM50_surv_c[,-1]
  coxFit_PAM50_surv_c <- coxph(Surv(time = Overall.Survival..Months.,
                                    event = Overall.Survival.Status)~factor(PAM50.Subtype),
                               x=T, y=T, model=T, method = "breslow",
                               data=PAM50_surv_c)
  summary(coxFit_PAM50_surv_c)
  # Concordance= 0.585  (se = 0.05 )
  # Rsquare= 0.019   (max possible= 0.606 )
  # Likelihood ratio test= 8.67  on 4 df,   p=0.07
  # Wald test            = 8.07  on 4 df,   p=0.09
  # Score (logrank) test = 8.6  on 4 df,   p=0.07

  # table(PAM50_surv_c$PAM50.Subtype)

  mfit_PAM50_surv_c <- survfit(Surv(time = Overall.Survival..Months.,
                       event = Overall.Survival.Status)~factor(PAM50.Subtype),
                  data=PAM50_surv_c)

  test <- summary(coxFit_PAM50_surv_c)
  ggsurvplot(mfit_PAM50_surv_c, conf.int = F, pval = test$sctest[3])

  PAM50_surv_c_no_normal <- which(PAM50_surv_c$PAM50.Subtype == "Normal-like")
  coxFit_PAM50_surv_c_no_normal <- coxph(Surv(time = Overall.Survival..Months.,
                                    event = Overall.Survival.Status)~factor(PAM50.Subtype),
                               x=T, y=T, model=T, method = "breslow",
                               data=PAM50_surv_c[-PAM50_surv_c_no_normal,])

  summary(coxFit_PAM50_surv_c_no_normal)

  mfit_PAM50_no_normal <- survfit(Surv(time = Overall.Survival..Months.,
                                    event = Overall.Survival.Status)~factor(PAM50.Subtype),
                               data=PAM50_surv_c[-PAM50_surv_c_no_normal,])

  test <- summary(coxFit_PAM50_surv_c_no_normal)
  ggsurvplot(mfit_PAM50_no_normal, conf.int = F, pval = test$sctest[3])

  table(PAM50_surv_c$PAM50.Subtype)

  setEPS()
  postscript("output4paper/KM_standard_subtypes.eps", width = 12, height = 9)
  #tiff(file = "output4paper/KM_standard_subtypes.tiff",  width = 12, height = 9, units = 'in', res=300)
  #pdf( file = "output4paper/KM_standard_subtypes.pdf",  onefile = TRUE, width = 12, height = 9)
  plot.Survival4paper(coxFit_PAM50_surv_c_no_normal, mfit_PAM50_no_normal, location = "topright",
                      colorsP = c("red","black","blue","#1E8449" ),
                      colorsL = c("#1E8449","blue","black","red" ),
                      labelClu = c("Luminal B (110)", "Luminal A (200)",  "HER2-enriched (54)", "Basal-like (85)") )
                      # font_size_times = 0.5)
  dev.off()
}


#PLOTS
plots <- function(){

  #Kaplan Meier population
  plot(survfit(coxFit), main = "clinical")
  plot(survfit(coxFit2), main = "NTF")
  plot(survfit(coxFit3), main = "NTF + clinical")
  plot(survfit(coxFit4_v2), main = "NMF")
  plot(survfit(coxFit_PAM50_surv_c), main = "PAM50")



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
                    Diagnosis.Age+#ER.Status+HER2.Status+
                    Metastasis.Coded+#PR.Status+
                    Tumor..T1.Coded,
                  data=survivalClinicalNTF)
  plot.SurvivalComplex (coxFit, mfit, "clinical variables")

  ####////////////
  #Using survminer library
  #Clinical
  ggsurvplot(mfit, conf.int = F, legend="none")
  ggrisktable(mfit, tables.theme = theme_cleantable(), color=NULL,
              fontsize = 1)

  #NTF components
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
  save.image("temp/79.3-result.RData")
}

removeUnwantedVars <- function(){
  X=ls()
  length(i <- grep("survival", X))
  X[i]
  length(j <- grep("plot", X))
  X[j]
  rm(list=ls()[-which(ls() %in% X[c(i,j)])])
}


#source
plotAnyTop <- function(threeGroups){ #survivalClinicalTop
  mfitTop <- survfit(Surv(time = Overall.Survival..Months.,
                          event = Overall.Survival.Status)~
                       factor(component),
                     data=survivalClinicalTop)

  coxFitTop <- coxph(Surv(time = Overall.Survival..Months.,
                          event = Overall.Survival.Status)~.,
                     data=survivalClinicalTop)
  print(summary(coxFitTop))
  list(mfitTop,coxFitTop)
}

plotAnyNTFforPaper <- function(threeGroups){
  # pdf( file = "temp/79-kaplan-NTF4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     colorsL = c("blue","black","red","green"),
  #                     labelClu = c("verotype 10","verotype 8", "verotype 2", "verotype 1") )
  # dev.off()
  # pdf( file = "temp/79-kaplan-NTFnl4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location="bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     labelClu = NULL)
  # dev.off()
  #pdf( file = "temp/79-kaplan-NTF4_num.pdf",  onefile = TRUE, width = 9, height = 7)
  tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
                      colorsP = c("red","black","blue","#1E8449"),
                      colorsL = c("#1E8449","blue","black","red"),
                      #labelClu = c("verotype 10","verotype 8", "verotype 2", "verotype 1") )
                      #1    9, 2  138, 8   99, 10   19
                      #labelClu = c("vero. 10 = 19","vero. 8 = 99", "vero. 2 = 138", "vero. 1 = 9") )
                      labelClu = c("subtype 1 (19)","subtype 2 (99)", "subtype 3 (138)", "subtype 4 (9)") )
  dev.off()
}

coxfromGroup <- function(Groups){
  survivalgroup <- surCli[which(as.numeric(surCli$component) %in% Groups),]
  newInstruction <- try(coxFit_G <- coxph(Surv(time = Overall.Survival..Months.,
                             Overall.Survival.Status)~.,
                        data=survivalgroup), silent=TRUE )
  if( inherits(newInstruction, "try-error") ){
    warning(paste("****** coxfromGroup ****** G = (", Groups,")", collapse = " "))
    return(Inf)
  }
  test <- summary(coxFit_G)
  test$sctest[3]
}

###


plotAnyNMFforPaper <- function(threeGroups){
  # pdf( file = "temp/79-kaplan-NMF4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     colorsL = c("blue","black","red","green"),
  #                     labelClu = c("verotype 10","verotype 8", "verotype 2", "verotype 1") )
  # dev.off()
  # pdf( file = "temp/79-kaplan-NMFnl4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location="bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     labelClu = NULL)
  # dev.off()
  #pdf( file = "temp/79-kaplan-NMF4_num.pdf",  onefile = TRUE, width = 9, height = 7)
  #png( file = "temp/79-kaplan-NMF4_num.png",  width = 680, height = 460)
  tiff(file = "temp/79-kaplan-NTF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
                      colorsP = c("red","black","blue","#1E8449"),
                      colorsL = c("#1E8449","blue","black","red"),
                      # x freq
                      # 1  1   68
                      # 2  7   50
                      # 3  8   58
                      # 4 10   42
                      #labelClu = c("group 10 (42)","group 8 (58)", "group 7 (50)", "group 1 (68)") )
                      labelClu = c("subtype 1 (42)","subtype 2 (58)", "subtype 3 (50)", "subtype 4 (68)") )
  dev.off()
}


plotAnyNMFforPaper_clini <- function(threeGroups){
  # pdf( file = "temp/79-kaplan-NMF4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     colorsL = c("blue","black","red","green"),
  #                     labelClu = c("verotype 10","verotype 8", "verotype 2", "verotype 1") )
  # dev.off()
  # pdf( file = "temp/79-kaplan-NMFnl4.pdf",  onefile = TRUE, width = 9, height = 7)
  # plot.Survival4paper(coxFitTop, mfitTop, location="bottomleft",
  #                     colorsP = c("red","black","blue","green"),
  #                     labelClu = NULL)
  # dev.off()
  #pdf( file = "temp/79-kaplan-NMF4_clini_num.pdf",  onefile = TRUE, width = 9.5, height = 7)
  #png(file = "temp/79-kaplan-NMF4_clini_num.png",  width = 680, height = 460)
  tiff(file = "temp/79-kaplan-NMF4_clini_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
                      colorsP = c("red","black","blue","#1E8449"),
                      colorsL = c("#1E8449","blue","black","red"),
                      # x freq
                      # 1  1   73
                      # 2  7   71
                      # 3  8   88
                      # 4 10   40
                      #labelClu = c("group 10 (40)","group 8 (88)", "group 7 (71)", "group 1 (73)") )
                      labelClu = c("subtype 1 (40)","subtype 2 (88)", "subtype 3 (71)", "subtype 4 (73)") )
  dev.off()
}


best_X_KM <- function(surCli){
  X <- 5
  surCli = patiAfi_NTF  #NTF
  #surCli = patiAfi_NMF_mut  #mutation
  #surCli = patiAfi_NMF_clin  #clinical
  #table(surCli$component)
  combi4 <- combn(kMax,X)
  allPvals4 <- lapply(c(1:dim(combi4)[2]),
                      function(x){
                        coxfromGroup(Groups=combi4[,x])})
  minpos <- which.min(unlist(allPvals4))
  bc4 = combi4[,minpos]
  minpos
  min(unlist(allPvals4))
  bc4
  count(surCli[which(surCli$component %in% bc4),]$component)

  survivalClinicalTop <- surCli[which(surCli$component %in% bc4),]
  testPlot4t <- plotAnyTop(bc4)
  mfitTop = testPlot4t[[1]]
  coxFitTop = testPlot4t[[2]]
  test <- summary(coxFitTop)
  ggsurvplot(mfitTop, conf.int = F,pval = test$sctest[3] )

  #p-val NMF_clin
  test <- summary(coxFit_Afi_NMF_clin)
  ggsurvplot(mfitTop, conf.int = F,pval = test$sctest[3] )


  plotAnyNTFforPaper(threeGroups)
}


#save.image("temp/79.RData")
