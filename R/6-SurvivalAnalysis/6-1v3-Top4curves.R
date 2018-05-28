#version 3 (May/28/2018)
#' Select the top 4 groups with best p-value
#' Pre: 79-Top3curvesFrom77_v2


rm(list = ls())
load("temp/6-1v3-data.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)


LoadMyData <- funtion()
{
  load("temp/5-1v3.RData")
  load("data/factors.RData")
  load("data/NMF_30_PxC.RData")
  load("data/NMF_30_PxM.RData")
  rm(list=ls()[-which(ls() %in% c("randPatfm","patientsF","weightsC","patients",
                                  "patientAfiliationNMF_PxC", "patientNMF_PxC", "patientNMF_PxM",
                                  "patientAfiliationNMF_PxM"))])
  load("data/survivalClinical-5-1_30.RData")
  load("data/tensorClinical.RData")
  row.names(tensorClinical) <- tensorClinical$patient.bcr_patient_barcode
  tensorClinical <- tensorClinical[patients,-1]
  load("data/myLib.RData")
  loadls("plyr survival missForest survAUC prodlim survminer", F)
  source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
  save.image(file="temp/6-1v3-data.RData")
}


#source
ntfPatientFactorMatrix <- function(){
  weightsComp <- data.frame(namesC = paste('V',seq(1:kMax),sep=''),
                            weight=weightsC[1,] )
  weightsComp <- weightsComp[order(weightsComp$weight,decreasing = F),]
  patiF <- patientsF
  for (i in seq_along(weightsC)) {
    patiF[,i] <- patiF[,i]*weightsC[i]
  }
  #normalize rows
  patiF <- t(apply(patiF,1,function(x){x/sum(x)}))
  patiF <- patiF[,as.character(weightsComp$namesC) ]
}

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


  patiF <- ntfPatientFactorMatrix()
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[,-1]
  #NFT components
  coxFit_NTF <- coxph(Surv(time = Overall.Survival..Months.,
                        event = Overall.Survival.Status)~.,
                   x=T, y=T, model=T, method = "breslow",
                   data=patiF)
  summary(coxFit_NTF)
  # Concordance= 0.726  (se = 0.054 )
  # Rsquare= 0.094   (max possible= 0.606 )
  # Likelihood ratio test= 44.84  on 29 df,   p=0.03
  # Wald test            = 46.14  on 29 df,   p=0.02
  # Score (logrank) test = 54.73  on 29 df,   p=0.003


  #NMF mutation
  #normalize rows
  pati_NMF_Mut <- data.frame( t(apply(patientNMF_PxM,1,function(x){x/sum(x)})) )
  names(pati_NMF_Mut) <- paste('V',seq(1:kMax),sep='')
  pati_NMF_Mut <- merge(as.data.frame(pati_NMF_Mut), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  pati_NMF_Mut <- pati_NMF_Mut[,-1]
  coxFit_NMF_mut <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~.,
                          x=T, y=T, model=T, method = "breslow",
                          data=pati_NMF_Mut)
  summary(coxFit_NMF_mut)
  # Concordance= 0.755  (se = 0.054 )
  # Rsquare= 0.085   (max possible= 0.606 )
  # Likelihood ratio test= 40.72  on 29 df,   p=0.07
  # Wald test            = 46.58  on 29 df,   p=0.02
  # Score (logrank) test = 54.54  on 29 df,   p=0.003


  #NMF clinical
  #normalize rows
  pati_NMF_Clin <- data.frame( t(apply(patientNMF_PxC,1,function(x){x/sum(x)})) )
  names(pati_NMF_Clin) <- paste('V',seq(1:kMax),sep='')
  pati_NMF_Clin <- merge(as.data.frame(pati_NMF_Clin), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  pati_NMF_Clin <- pati_NMF_Clin[,-1]

  coxFit_NMF_cli <- coxph(Surv(time = Overall.Survival..Months.,
                               event = Overall.Survival.Status)~.,
                          x=T, y=T, model=T, method = "breslow",
                          data=pati_NMF_Clin)
  # Warning message:
  #   In fitter(X, Y, strats, offset, init, control, weights = weights,  :
  #               Loglik converged before variable  12,28 ; beta may be infinite.
  summary(coxFit_NMF_cli)
  # Concordance= 0.75  (se = 0.054 )
  # Rsquare= 0.077   (max possible= 0.606 )
  # Likelihood ratio test= 36.35  on 27 df,   p=0.1
  # Wald test            = 37.17  on 27 df,   p=0.09
  # Score (logrank) test = 44.02  on 27 df,   p=0.02


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

  save.image("temp/6-1v3-result.RData")
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
plotAnyTopNTF <- function(threeGroups){
  mfitTop <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                          factor(component),
                        data=survivalClinicalTop)

  coxFitTop <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                         factor(component),
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

coxfromGroupNTF <- function(Groups){
  survivalgroup <- surCli[which(surCli$component %in% Groups),]
  coxFit_G <- coxph(Surv(time = Overall.Survival..Months.,
                         Overall.Survival.Status)~.,
                       data=survivalgroup)
  test <- summary(coxFit_G)
  test$sctest[3]
}

best4NTF <- function(){
  surCli = patiF
  # factor 'component'
  combi4 <- combn(kMax,4)
  allPvals4 <- lapply(c(1:dim(combi4)[2]),
                      function(x){
                        coxfromGroupNTF(combi4[,x])})
  which.min(unlist(allPvals4))
  #27
  min(unlist(allPvals4))
  bc4 = combi4[,27]
  bc4
  # >   min(unlist(allPvals4))
  # [1] 0.2405323
  # > bc4 = combi4[,27]
  # > bc4
  # [1]  1  2  8 10
  count(surCli[which(surCli$component %in% threeGroups),]$component)

  threeGroups = bc4
  survivalClinicalTop <- surCli[which(surCli$component %in% threeGroups),]
  testPlot4t <- plotAnyTopNTF(bc4)
  mfitTop = testPlot4t[[1]]
  coxFitTop = testPlot4t[[2]]
  ggsurvplot(mfitTop, conf.int = F)

  plotAnyNTFforPaper(threeGroups)
}

###

plotAnyTopNMF <- function(threeGroups){
  mfitTop <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                       factor(C),
                     data=survivalClinicalTop)

  coxFitTop <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                       factor(C),
                     data=survivalClinicalTop)
  print(summary(coxFitTop))
  list(mfitTop,coxFitTop)
}

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

coxfromGroupNMF <- function(Groups){
  survivalgroup <- surCli[which(surCli$C %in% Groups),]
  coxFit_G <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~
                      factor(C),
                    data=survivalgroup)
  test <- summary(coxFit_G)
  test$sctest[3]
}

best4NMF <- function(){
  surCli = survivalClinicalNMF_2
  # factor 'C'
  combi4 <- combn(10,4)
  allPvals4 <- lapply(c(1:dim(combi4)[2]),
                      function(x){
                        coxfromGroupNMF(combi4[,x])})
  which.min(unlist(allPvals4))
  #82
  min(unlist(allPvals4))
  bc4 = combi4[,82]
  bc4
  # >   min(unlist(allPvals4))
  # [1] 0.04928328
  # >   bc4 = combi4[,82]
  # >   bc4
  # [1]  1  7  8 10

  threeGroups = bc4
  count(surCli[which(surCli$C %in% threeGroups),]$C)


  survivalClinicalTop <- surCli[which(surCli$C %in% threeGroups),]
  testPlot4t <- plotAnyTopNMF(bc4)
  mfitTop = testPlot4t[[1]]
  coxFitTop = testPlot4t[[2]]
  ggsurvplot(mfitTop, conf.int = F)

  plotAnyNMFforPaper(threeGroups)

}

###

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

best4NMF_clini <- function(){
  surCli = survivalClinicalNMF_clini_2
  # factor 'C'
  combi4 <- combn(10,4)
  allPvals4 <- lapply(c(1:dim(combi4)[2]),
                      function(x){
                        coxfromGroupNMF(combi4[,x])})
  which.min(unlist(allPvals4[-c(87,141,140, 139, 138, 137, 136)]))
  #87
  min(unlist(allPvals4[-c(87,141,140, 139, 138, 137, 136)]))
  bc4 = combi4[,82]
  bc4
  # >   bc4
  # [1]  1  7  8 10

  threeGroups = bc4
  count(surCli[which(surCli$C %in% threeGroups),]$C)

  survivalClinicalTop <- surCli[which(surCli$C %in% threeGroups),]
  testPlot4t <- plotAnyTopNMF(bc4)
  #p=0.2073
  mfitTop = testPlot4t[[1]]
  coxFitTop = testPlot4t[[2]]
  ggsurvplot(mfitTop, conf.int = F)

  plotAnyNMFforPaper_clini(threeGroups)

}

#save.image("temp/79.RData")
