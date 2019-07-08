rm(list = ls())
load("data/loadls.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)

library("purrr")
library("cluster")

rm(list = ls())
load("dataExample/survival/survivalClinical.RData")
load(file="temp/6survival.RData")

patients <- row.names(survivalClinical)
row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)
survClinical <- survClinical[patients,]

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


  patiF <- data.frame(patients,survClinical$PAM50.Subtype) 
  names(patiF) <- c("patients", "group")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[-c(1,2),]
  row.names(patiF) <- patiF[,1]

  #remove NA
  #patiF <- patiF[-c(1,2),]
  anyNA(patiF)
  #remNA = which(is.na(patiF$Overall.Survival.Status))
  #patiF <- patiF[-remNA,]


  #remove groups with one patient
  table(patiF$group)
  removeG <- which(table(patiF$group)<= 3)
  if (length(removeG) >= 1){
    patiF <- patiF[-(which(patiF$group %in% removeG)),]
    table(patiF$group)
  }

  # remove groups with Normal-like patients e.g. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4011983/
  # "Tumors categorized by PAM50 as Normal-like, a subtype with no corresponding clinicopathological category, were excluded from the kappa calculation." 
  patiF <- patiF[-(which(patiF$group == "Normal-like")),]
  
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
  km4 <- prodlim(Surv(time = Overall.Survival..Months., event = Overall.Survival.Status)~
                   factor(group),
                 data=patiF)
  plot(km4, logrank=TRUE)


  #Plot for paper
  main.i=i
  tr = table(factor(patiF$group))
  lr = levels(factor(patiF$group))
  i = length(lr)
  typesS = levels(factor(patiF$group))
  pdf( file = paste0("temp/KM_standard_subtypes.pdf"),  onefile = TRUE, width = 9, height = 7)
  # #tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  plot.Survival4paper(coxFit_NBS, mfit, location = "bottomleft",
                      colorsP = c25[1:i],
                      colorsL = c25[1:i],
                      labelClu = paste(typesS, " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray95")
  dev.off()
