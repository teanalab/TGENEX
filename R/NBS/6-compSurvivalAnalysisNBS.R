rm(list = ls())
library(purrr)
library(cluster)
library(survival)
library(prodlim)

#nice colors
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
#loat plotting functions
source("R/6-plotSurvivalUpdate.R")


survAna <- function(k){
  patiF <- data.frame(patients,affiliationList[patients] )
  names(patiF) <- c("patients", "group")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[-c(1,2),]
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2)]
  
  
  #remove NA
  #patiF <- patiF[-c(1,2),]
  anyNA(patiF)
  remNA = which(is.na(patiF$Overall.Survival.Status))
  patiF <- patiF[-remNA,]
  
  
  #remove groups with one patient
  table(patiF$group)
  removeG <- which(table(patiF$group)<= 3)
  if (length(removeG) >= 1){
    patiF <- patiF[-(which(patiF$group %in% removeG)),]
    table(patiF$group)
  }
  
  patiF <<- patiF
  #levels(factor(patiF$group))
  
  # components
  SurFunction <- Surv(time = Overall.Survival..Months.,
                      event = Overall.Survival.Status)~factor(group)
  coxFit_NBS <- coxph(SurFunction,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiF)
  summary(coxFit_NBS)
  mfit <- survfit(SurFunction, data = patiF)
  
  
  ####
  #complete plot
  # km4 <- prodlim(Surv(time = Overall.Survival..Months., event = Overall.Survival.Status)~
  #                  factor(group),data=patiF)
  # plot(km4, logrank=TRUE )
  
  ####
  #Plot for paper
  main.k=k
  tr = table(patiF$group)
  lr = levels(as.factor(patiF$group))
  k = length(lr)
  pdf( file = paste0("temp/kaplan-bestNBS_k",main.k,".pdf"),  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.5)
  dev.off()
}


#Figure 4.B) best subtypes obtained with NBS
#
load("dataExample/survival/survivalClinical.RData")
k=6
load(file=paste("dataExample/survival/affiliationNBS_k",k,".RData",sep='' ) )
patients <- names(affiliationList)
row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)
survClinical <- survClinical[patients,]
survivalClinical <- survClinical[,c( "patient.vital_status", "OS_MONTHS")]
survivalClinical$patient.vital_status = survivalClinical$patient.vital_status != "alive"
names(survivalClinical) <- c( "Overall.Survival.Status", "Overall.Survival..Months.")
survAna(k)



#Plot all clusters ##
# load("dataExample/survival/survivalClinical.RData")
# survClinicalAll <- survClinical
# 
# for(k in c(2:15)){
#   load(file=paste("dataExample/survival/affiliationNBS_k",k,".RData",sep='' ) , globalenv())
#   patients <- names(affiliationList)
#   survClinicalAll -> survClinical
#   row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)
#   survClinical <- survClinical[patients,]
#   survivalClinical <- survClinical[,c( "patient.vital_status", "OS_MONTHS")]
#   survivalClinical$patient.vital_status = survivalClinical$patient.vital_status != "alive"
#   names(survivalClinical) <- c( "Overall.Survival.Status", "Overall.Survival..Months.")
#   
#   survAna(k)
# }
