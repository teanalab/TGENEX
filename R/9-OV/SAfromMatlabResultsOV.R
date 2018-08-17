#version 2 (Aug/2/2018)
rm(list = ls())
load("loadls.RData")
loadls("plyr survival missForest survAUC prodlim survminer RTCGA.clinical", F)


LoadData <- function()
{
  clustersNBS <<- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_output_clustersNBS4.csv", sep = "\t", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  standardNMF <<- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_output_standardNMF_indClust_ctrl4.csv", sep = ",", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  patients <<-    read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_output_sample_id4.csv", sep = "\t", quote = "'",
                        header = FALSE, stringsAsFactors = FALSE)
  clustersNBS_cc <<- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_CC_clustersNBS4.csv", sep = "\t", quote = "",
                            header = FALSE, stringsAsFactors = FALSE)
  clusters_cc <<- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_CC_4.csv", sep = "\t", quote = "",
                               header = FALSE, stringsAsFactors = FALSE)
  ###
  #More details on Phenotype For OV.R
  survivalTCGA(OV.clinical) -> OV.survInfo
  ppp <<- intersect( tolower(as.character(OV.survInfo$bcr_patient_barcode)), as.character(unlist(patients)) )
  pheno_OV <<- OV.survInfo[which(tolower(as.character(OV.survInfo$bcr_patient_barcode)) %in% ppp),]

  names(pheno_OV) <<- c( "months", "sample", "survi")
  pheno_OV$sample <<- tolower(pheno_OV$sample)
  row.names(pheno_OV) <<- pheno_OV$sample

  #NMF R
  load(file="data/NMF_4_PxC.RData", envir = globalenv())
  load(file="data/NMF_k4_PxM.RData", envir = globalenv())


  #Rubik
  load(file="data/patientAfiliationRubik_PxC.Rd", envir = globalenv())

  ##plotting functions
  source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
  #save.image(file="temp-SA-data.RData")
}

LoadData()

#' Survival Analysis NBS
survAnaNBS <- function(){
  patiF <- data.frame(patients,clustersNBS)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2)]
  patiF <- patiF[,-3]

  toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  anyNA(patiF)


  patiF$months <- as.numeric(patiF$months)

  #NBS components
  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_NBS <- coxph(Surv(time = months,
                           event = survi)~factor(cluster),
                      method = "breslow", data=patiF,
                      x=T, y=T, model=T)
  x <- summary(coxFit_NBS)
  # Concordance= 0.588  (se = 0.024 )
  # Rsquare= 0.04   (max possible= 0.988 )
  # Likelihood ratio test= 14.27  on 3 df,   p=0.003
  # Wald test            = 12.25  on 3 df,   p=0.007
  # Score (logrank) test = 13.28  on 3 df,   p=0.004

  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_NBS, mfit, "NBS")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return (x)
}

pValNBS <-survAnaNBS()

#' Survival Analysis NMF
survAnaNMF <- function(){
  patiF <- data.frame(patients,standardNMF)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]


  toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  anyNA(patiF)


  patiF$months <- as.numeric(patiF$months)

  #NBS components
  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_NBS <- coxph(Surv(time = months,
                           event = survi)~factor(cluster),
                      method = "breslow", data=patiF,
                      x=T, y=T, model=T)
  x <- summary(coxFit_NBS)
  # Concordance= 0.513  (se = 0.011 )
  # Rsquare= 0.012   (max possible= 0.988 )
  # Likelihood ratio test= 4.25  on 3 df,   p=0.2
  # Wald test            = 5.38  on 3 df,   p=0.1
  # Score (logrank) test = 5.76  on 3 df,   p=0.1


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_NBS, mfit, "NBS")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}

pValNMF <-survAnaNMF()


#' Survival Analysis NBS cc
survAnaNBS_cc <- function(){
  patiF <- data.frame(patients,clustersNBS_cc)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]


  toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  anyNA(patiF)


  patiF$months <- as.numeric(patiF$months)

  #NBS components
  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_NBS <- coxph(Surv(time = months,
                           event = survi)~factor(cluster),
                      method = "breslow", data=patiF,
                      x=T, y=T, model=T)
  x <- summary(coxFit_NBS)
  # Concordance= 0.579  (se = 0.023 )
  # Rsquare= 0.042   (max possible= 0.988 )
  # Likelihood ratio test= 15.04  on 3 df,   p=0.002
  # Wald test            = 11.96  on 3 df,   p=0.008
  # Score (logrank) test = 13.27  on 3 df,   p=0.004


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_NBS, mfit, "NBS")

  fit <- survfit(Surv(months, survi) ~ factor(cluster),
                 data = patiF)
  # Visualize with survminer
  ggsurvplot(fit, data = patiF, risk.table = TRUE)


  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}
pValNBS_cc <-survAnaNBS_cc()


#' Survival Analysis cc
survAna_cc <- function(){
  patiF <- data.frame(patients,clusters_cc)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]


  toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  anyNA(patiF)


  patiF$months <- as.numeric(patiF$months)

  #NBS components
  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_NBS <- coxph(Surv(time = months,
                           event = survi)~factor(cluster),
                      method = "breslow", data=patiF,
                      x=T, y=T, model=T)
  x <- summary(coxFit_NBS)
  # Concordance= 0.501  (se = 0.011 )
  # Rsquare= 0.004   (max possible= 0.988 )
  # Likelihood ratio test= 1.57  on 3 df,   p=0.7
  # Wald test            = 1.92  on 3 df,   p=0.6
  # Score (logrank) test = 2.01  on 3 df,   p=0.6


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_NBS, mfit, "NBS")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}
pVal_cc <-survAna_cc()


#' Survival Analysis NMFr_m
survAna_NMFr_m <- function(){
  patiF <- data.frame(patientAfiliationNMF_PxM)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(patiF$patients)
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]

  anyNA(patiF)
  #toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  patiF$months <- as.numeric(patiF$months)

  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_ <- coxph(Surv(time = months,
                           event = survi)~factor(cluster),
                      method = "breslow", data=patiF,
                      x=T, y=T, model=T)
  x <- summary(coxFit_)
  # Concordance= 0.511  (se = 0.013 )
  # Rsquare= 0.004   (max possible= 0.988 )
  # Likelihood ratio test= 1.36  on 3 df,   p=0.7
  # Wald test            = 1.52  on 3 df,   p=0.7
  # Score (logrank) test = 1.55  on 3 df,   p=0.7


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_, mfit, "NBS")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}
pVsurvAna_NMFr_m <- survAna_NMFr_m()

#' Survival Analysis NMFr_c
survAna_NMFr_c <- function(){
  patiF <- data.frame(patientAfiliationNMF_PxC)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(patiF$patients)
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]

  anyNA(patiF)
  #toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  patiF$months <- as.numeric(patiF$months)

  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_ <- coxph(Surv(time = months,
                        event = survi)~factor(cluster),
                   method = "breslow", data=patiF,
                   x=T, y=T, model=T)
  x <- summary(coxFit_)
  # Concordance= 0.511  (se = 0.013 )
  # Rsquare= 0.004   (max possible= 0.988 )
  # Likelihood ratio test= 1.36  on 3 df,   p=0.7
  # Wald test            = 1.52  on 3 df,   p=0.7
  # Score (logrank) test = 1.55  on 3 df,   p=0.7


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_, mfit, "Survival")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}
pVsurvAna_NMFr_c <- survAna_NMFr_c()


#' Survival Analysis Rubik
survAna_Rubik <- function(){
  patiF <- data.frame(patientAfiliationRubik_PxC)
  names(patiF) <- c("patients", "cluster")
  row.names(patiF) <- as.character(patiF$patients)
  patiF <- patiF[ppp,]
  pheno_OV <- pheno_OV[ppp, ]
  patiF <- merge(as.data.frame(patiF), as.data.frame(pheno_OV), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2,5)]

  anyNA(patiF)
  #toDel <- which(is.na(patiF$months))
  #patiF <- patiF[-11,]

  patiF$months <- as.numeric(patiF$months)

  SurFunction <- Surv(time = patiF$months,
                      event = patiF$survi)~factor(patiF$cluster)
  coxFit_ <- coxph(Surv(time = months,
                        event = survi)~factor(cluster),
                   method = "breslow", data=patiF,
                   x=T, y=T, model=T)
  x <- summary(coxFit_)
  # Concordance= 0.545  (se = 0.024 )
  # Rsquare= 0.019   (max possible= 0.988 )
  # Likelihood ratio test= 6.62  on 3 df,   p=0.08
  # Wald test            = 7.78  on 3 df,   p=0.05
  # Score (logrank) test = 8.09  on 3 df,   p=0.04


  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(Surv(time = months,
                      event = survi)~factor(cluster),
                 method = "breslow", data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction)
  plot.SurvivalComplex (coxFit_, mfit, "Survival")

  # plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
  #                     colorsP = c("red","black","blue", "green"),
  #                     colorsL = c("blue","black","red", "green"),
  #                     labelClu = c("subtype 1","subtype 2", "subtype 3", "green") )
  return(x)
}
pVsurvAna_Rubik <- survAna_Rubik()


cat("results: \n",
    "\nNBS pvalue \t", pValNBS$sctest[3],
    "\nNMF pvalue \t", pValNMF$sctest[3],
    "\nNBS_cc pvalue \t", pValNBS_cc$sctest[3],
    "\nCC pvalue \t", pVal_cc$sctest[3],
    "\nsurvAna_NMFr_m \t", pVsurvAna_NMFr_m$sctest[3],
    "\nsurvAna_NMFr_c \t", pVsurvAna_NMFr_c$sctest[3],
    "\nsurvAna_Rubik \t", pVsurvAna_Rubik$sctest[3])


save.image("temp/SA.RData")
