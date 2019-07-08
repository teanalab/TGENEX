
rm(list = ls())
load("data/loadls.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)


load("dataExample/survival/survivalClinical.RData")
rm(list = ls()[-which(ls() %in% c("survClinical")  )])

save(survClinical, file="dataExample/survival/survClinical.RData")

anyNA(survClinical)
row.names(survClinical)
ppp <- as.character(survClinical$patient.bcr_patient_barcode)
survClinical[remNA,]

length(intersect(ppp,patients))


binaMredu <-read.table(file="dataExample/fromNBS/10-binaM.csv",sep = ";", quote = '"',
                       header = TRUE, stringsAsFactors = FALSE)
clustersNBS <-read.table(file="dataExample/fromNBS/clustersNBS_woCC_k3.csv", sep = "\t", quote = "",
                         header = FALSE, stringsAsFactors = FALSE)
standardNMF <-read.table(file="dataExample/fromNBS/standardNMF_k4.csv", sep = ",", quote = "",
                         header = FALSE, stringsAsFactors = FALSE)
genes_redu <-read.table(file="dataExample/fromNBS/genes.csv", sep = "\t", quote = '"',
                        header = TRUE, stringsAsFactors = FALSE)
clinical <-read.table(file="dataExample/fromNBS/10-binac.csv",sep = ";", quote = '"',
                      header = TRUE, stringsAsFactors = FALSE)
load("dataExample/Filtered/patients.RData")

row.names(survClinical) <- as.character(survClinical$patient.bcr_patient_barcode)
survivalClinical <- survClinical[,c( "patient.vital_status", "OS_MONTHS")]
survivalClinical$patient.vital_status = survivalClinical$patient.vital_status != "alive"
names(survivalClinical) <- c( "Overall.Survival.Status", "Overall.Survival..Months.")
survivalClinical <- survivalClinical[as.character(unlist(patients)),]
rm(list = ls()[which(ls() %in% c("survClinical")  )])

affiliationList = clustersNBS
affiliationList = unlist(affiliationList)
names(affiliationList) = as.character(unlist(patients))
save(affiliationList, file = "dataExample/afiliations/affiliationNBSk3.RData" )


#plotting functions
source("R/6-plotSurvivalUpdate.R")
#save.image(file="temp/6-8v1-data.RData")


#' Survival Analysis
survAna <- function(patients, group, survivalClinical){
  patiF <- data.frame(patients,clustersNBS)
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2)]
  
  #remove NA
  #patiF <- patiF[-c(1,2),]
  if(anyNA(patiF)){
    remNA = which(is.na(patiF$Overall.Survival.Status))
    patiF <- patiF[-remNA,]
  }

  
  # #Kaplan Meier all samples
  # plot(survfit(coxFit_NBS), main = "NBS")
  
  
  #NBS components
  SurFunction <- Surv(time = Overall.Survival..Months.,
                      event = Overall.Survival.Status)~V1
  coxFit_NBS <- coxph(SurFunction,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiF)
  summary(coxFit_NBS)
  test <- summary(coxFit_NBS)
  test$sctest[3]
  
  
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
  clustersNBSks <- data.frame(patients, matrix(NA, nrow=length(patients), ncol = 13))
  
  clustersNBSks[,2] <- unlist(clustersNBS)
  
  for (i in c(3:14)) {
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
  
  #colorblind palette
  # library("ggplot2")
  # library("scales")
  # library("ggthemes")
  # show_col(colorblind_pal()(8))
  # p <- ggplot(mtcars) + geom_point(aes(x = wt, y = mpg,
  #                                      colour = factor(gear))) + facet_wrap(~am)
  # p + theme_igray() + scale_colour_colorblind()
  # https://rdrr.io/cran/ggthemes/man/colorblind.html
  
  
  #better colors
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
  pie(rep(1,25), col=c25)
  
  for (i in c(2:14)) {
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
                        colorsP = c25[1:(i+1)],
                        colorsL = c25[1:(i+1)],
                        labelClu = paste("subtype ", c(1:(i+1)) ),
                        font_size_times = 1)
    dev.off()
  }
}

save.image("temp/6-8.RData")
