#Integrating and final survival

library(plyr)
library(survival)
library(classInt) #for sturges
source("R/6-plotSurvivalUpdate.R") #for survival plot

# parameters ##########
k=3
minimumGroupSize = 4
foldermRNA="R/clean_mRNA/"
#1. get list of diseases
diseases = list.files(path=foldermRNA)
diseases = gsub(".RData", "", diseases, fixed = TRUE)

#loadAllmRNA #####
# for (i in seq_along(diseases))
# {
#   disease = diseases[i]
#   load(paste0(foldermRNA, "/", disease,".RData"))
# }


plotDendoSurvival <- function(disease, outFolder = "kmplots_from_dendo/"){
  #load survival data
  load(paste0("survival/",disease,"/survival.Rd"))

  #load mRNA data
  load(paste0(foldermRNA, "/", disease,".RData"))
  eval(parse(text= paste0("TCGA_DB <- ",disease,".mRNA")))

  # find intersection between survival and mRNA
  patients <- row.names(TCGA_DB)
  p_s <-survivalData$bcr_patient_barcode
  patients <- intersect(patients, p_s)

  # update tables
  survivalData <- survivalData[patients,]
  TCGA_DB <- TCGA_DB[patients,]


  if ( anyNA(TCGA_DB) ){
    stop("ERROR: Table TCGA_DB have NA values for ", disease)
  }
  if ( anyNA(survivalData) ){
    stop("ERROR: Table survivalData have NA values for ", disease)
  }

  cat( "For", disease, "Number of genes", dim(TCGA_DB)[2], "Number of patients", dim(TCGA_DB)[1] , "\n")

  ## Row- and column-wise clustering
  hr <- hclust(as.dist(1-cor(t(TCGA_DB), method="pearson")), method="complete")
  # get groups
  survivalData$group <- as.numeric(cutree(hr, k = k))
  distSurvGroups <- table(survivalData$group)
  print(distSurvGroups)

  if (any(distSurvGroups < minimumGroupSize)){
    cat("ERROR: one or more groups are smaller than the treshold of ",
        minimumGroupSize, "adjusting k from", k, "to", k-1, "FOR disease", disease)
    # get groups
    survivalData$group <- as.numeric(cutree(hr, k = k-1))
    distSurvGroups <- table(survivalData$group)
    print(distSurvGroups)
    if (any(distSurvGroups < minimumGroupSize)){
      cat("ERROR AFTER ADJUSTMENT: one or more groups are smaller than the treshold of ",
          minimumGroupSize, "with k = ", k-1, "FOR disease", disease, "REMOVING SMALLER GROUP")

      survivalData$group <- as.numeric(cutree(hr, k = k))
      distSurvGroups <- table(survivalData$group)
      print(distSurvGroups)

      smallerGNr <- which(distSurvGroups < minimumGroupSize)
      survivalData <- survivalData[-which(survivalData$group == smallerGNr),]
      distSurvGroups <- table(survivalData$group)
      print(distSurvGroups)
      if ( any(distSurvGroups < minimumGroupSize) ){
        stop(paste("ERROR AFTER ADJUSTMENT: one or more groups are smaller than the treshold of ",
                   minimumGroupSize, "with k = ", k-1, "FOR disease", disease))
      }
    }
  }


  # survival analysis
  cat("\n Performing survival analysis for ", disease)
  # components
  SurFunction <- Surv(time = times,
                      event = patient.vital_status)~as.factor(group)
  coxFit <- coxph(SurFunction,
                  x=T, y=T, model=T, method = "breslow",
                  data=survivalData)
  cat(". Cox Log p-val = ",format.pval(summary(coxFit)$sctest[3], digits = 3) )
  mfit <- survfit(SurFunction, data = survivalData)

  # FINALLY, plot
  tr = table(survivalData$group)
  fileName = paste0(outFolder,disease,".pdf")
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "topright",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()
  #crop borders to fit paper
  system(paste("pdfcrop", fileName, fileName))

  save(survivalData, file = paste0(outFolder,disease,"_survival.RData"))

}





plotIntegraSurvival <- function(disease, outFolder = "kmplots_f/", outFolderSurv = "kmplots/", outFolderDendo = "kmplots_from_dendo/"){
  #load survival dendo
  load(paste0(outFolderDendo,disease,"_survival.RData"))
  s_dendo <- survivalData

  #load survival
  load(paste0(outFolderSurv,disease,"_survival.RData"))

  cat( "For", disease, "\n")

  distSurvGroups <- table(survivalData$group)
  print(distSurvGroups)

  if (any(distSurvGroups < minimumGroupSize)){
    cat("ERROR: one or more groups are smaller than the treshold of ",
        minimumGroupSize, "adjusting k from", k, "to", k-1, "FOR disease", disease, "REMOVING SMALLER GROUP")

    smallerGNr <- which(distSurvGroups < minimumGroupSize)
    survivalData <- survivalData[-which(survivalData$group == smallerGNr),]
    distSurvGroups <- table(survivalData$group)
    print(distSurvGroups)
    if ( any(distSurvGroups < minimumGroupSize) ){
      stop(paste("ERROR AFTER ADJUSTMENT: one or more groups are smaller than the treshold of ",
                 minimumGroupSize, "with k as shown in table FOR disease", disease))
    }
  }

  # survival analysis
  cat("\n Performing survival analysis for ", disease)
  # components
  SurFunction <- Surv(time = times,
                      event = patient.vital_status)~as.factor(group)
  coxFit <- coxph(SurFunction,
                  x=T, y=T, model=T, method = "breslow",
                  data=survivalData)
  cat(". Cox Log p-val = ",format.pval(summary(coxFit)$sctest[3], digits = 3) )
  mfit <- survfit(SurFunction, data = survivalData)

  # FINALLY, plot
  tr = table(survivalData$group)
  fileName = paste0(outFolder,disease,"_tr.pdf")
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "topright",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()

  #crop borders to fit paper
  #system(paste("pdfcrop", fileName, fileName))


  fileName = paste0(outFolder,disease,"_br.pdf")
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "bottomright",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()

  #crop borders to fit paper
  #system(paste("pdfcrop", fileName, fileName))


  fileName = paste0(outFolder,disease,"_bl.pdf")
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "bottomleft",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()

  #crop borders to fit paper
  #system(paste("pdfcrop", fileName, fileName))


  save(survivalData, file = paste0(outFolder,disease,"_survival.RData"))


}

