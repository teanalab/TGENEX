#version 3 (May/28/2018)
rm(list = ls())
load("data/myLib.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)


LoadMyData <- funtion()
{
  #old data
  load("data/survivalClinical-5-1_30.RData")
  # load("data/tensorClinical.RData")
  # row.names(tensorClinical) <- tensorClinical$patient.bcr_patient_barcode
  # tensorClinical <- tensorClinical[patients,-1]


  binaMredu <-read.table(file="nbs_release_v0.2/BCRAdata/10-binaMredu.csv",sep = ";", quote = '"',
                         header = TRUE, stringsAsFactors = FALSE)
  clustersNBS <-read.table(file="nbs_release_v0.2/BCRAdata/output_clustersNBS.csv", sep = "\t", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  standardNMF <-read.table(file="nbs_release_v0.2/BCRAdata/output_standardNMF.csv", sep = ",", quote = "",
                           header = FALSE, stringsAsFactors = FALSE)
  genes_redu <-read.table(file="nbs_release_v0.2/BCRAdata/10-genes_redu.csv", sep = "\t", quote = '"',
                          header = TRUE, stringsAsFactors = FALSE)
  clinical <-read.table(file="nbs_release_v0.2/BCRAdata/10-clinical.csv",sep = ";", quote = '"',
                        header = TRUE, stringsAsFactors = FALSE)
  patients <-read.table(file="nbs_release_v0.2/BCRAdata/10-patients.csv", sep = "\t", quote = '"',
                        header = TRUE, stringsAsFactors = FALSE)

  survivalClinical <- survivalClinical[as.character(unlist(patients)),]

  #plotting functions
  source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")
  save.image(file="temp/6-8v1-data.RData")
}


#' Survival Analysis
survAna <- function(){
  patiF <- data.frame(patients,clustersNBS)
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-1]
  patiF <- patiF[,-1]
  #NBS components
  SurFunction <- Surv(time = Overall.Survival..Months.,
                      event = Overall.Survival.Status)~V1
  coxFit_NBS <- coxph(SurFunction,
                        x=T, y=T, model=T, method = "breslow",
                        data=patiF)
  summary(coxFit_NBS)
  # Concordance= 0.536  (se = 0.05 )
  # Rsquare= 0.001   (max possible= 0.606 )
  # Likelihood ratio test= 0.41  on 1 df,   p=0.5
  # Wald test            = 0.41  on 1 df,   p=0.5
  # Score (logrank) test = 0.41  on 1 df,   p=0.5

  #Kaplan Meier population
  plot(survfit(coxFit_NBS), main = "NBS")

  #Nonparametric estimation in event history analysis
  #clinical variables
  km0 <- prodlim(SurFunction,
                 method = "breslow",
                 data=patiF)
  #Complex Plots
  mfit <- survfit(SurFunction,
                  data=patiF)
  plot.SurvivalComplex (coxFit_NBS, mfit, "NBS")

  #
  plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
                      colorsP = c("red","black","blue"),
                      colorsL = c("blue","black","red"),
                      labelClu = c("subtype 1","subtype 2", "subtype 3") )
}


LoadMyDataSeveralKs <- funtion()
{
  clustersNBSks <- data.frame(patients, matrix(NA,nrow=length(patients), ncol = 12))

  clustersNBSks[,2] <- unlist(clustersNBS)

  for (i in c(3:13)) {
    fileName <- paste("nbs_release_v0.2/BCRAdata/output_clustersNBS_k",i+1,".csv",sep='')
    clustersNBSks[,i] <-read.table(file=fileName, sep = "\t", quote = "",
                               header = FALSE, stringsAsFactors = FALSE)
  }

  #verifying correct reading
  apply(clustersNBSks, 2, max)
  anyNA(clustersNBSks)

  save.image(file="temp/6-8v1-data2.RData")
}


survAnaDifferentKs <- function(){
  figureKs <- "temp/SA_NBS_ks"
  #dir.create(figureKs)

  for (i in c(2:13)) {
    rm(list=ls()[which(ls() %in% c("SurFunction", "coxFit_NBS", "mfit" ))])
    clustersNBSi <- clustersNBSks[,i]
    patiFi <- data.frame(patients,clustersNBSi)
    row.names(patiFi) <- as.character(unlist(patients))
    patiFi <- merge(as.data.frame(patiFi), as.data.frame(survivalClinical), by='row.names', all=TRUE)
    row.names(patiFi) <- patiFi[,1]
    patiFi <- patiFi[,-1]
    patiFi <- patiFi[,-1]
    #NBS components
    SurFunction <- Surv(time = patiFi$Overall.Survival..Months.,
                        event = patiFi$Overall.Survival.Status)~patiFi$clustersNBSi
    coxFit_NBS <- coxph(SurFunction,
                        x=T, y=T, model=T, method = "breslow")
    mfit <- survfit(SurFunction)

    pdf( file = paste0(figureKs,"/kaplan-NBS",(i+1),".pdf"),  onefile = TRUE, width = 9, height = 7)
    #tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
    plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
                        colorsP = c("red","black","blue"),
                        colorsL = c("blue","black","red"),
                        labelClu = c("subtype 1","subtype 2", "subtype 3") )
    dev.off()
  }
}

save.image("temp/6-8.RData")
