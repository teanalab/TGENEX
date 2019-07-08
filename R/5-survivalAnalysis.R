library("plyr")
library("survival")
library("missForest")
library("survAUC")
library("prodlim")
library("survminer")
library("gsubfn")


survivalAnalysis <- function(survivalClinical,r, k, threshold_groups=3, plots=F, 
                             plotName = paste0("temp/plotSurvCurves_r_",r,"k_",k,".pdf"),
                             location = "bottomleft")
{
  if(plots) {
    source("R/6-plotSurvivalUpdate.R")
  }
    
  patients <- names(affiliationList)
  patiF <- data.frame(patients, affiliationList[patients] )
  names(patiF) <- c("patients", "group")
  row.names(patiF) <- as.character(unlist(patients))
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  row.names(patiF) <- patiF[,1]
  patiF <- patiF[,-c(1,2)]
  
  #remove NA
  if(anyNA(patiF)) {
    remNA = which(is.na(patiF$Overall.Survival.Status))
    patiF <- patiF[-remNA,]
  }
  
  #remove groups with one patient
  #table(patiF$group)
  removeG <- which(table(patiF$group) <= threshold_groups)
  if (length(removeG) >= 1){
    patiF <- patiF[-(which(patiF$group %in% removeG)),]
    table(patiF$group)
  }
  
  tr = table(patiF$group)
  lr = levels(factor(patiF$group))
  patiF <- patiF
  nFactors <- as.numeric(length(lr))
  ok <- TRUE
  
  if (nFactors == 1){
    warning ( paste0("ONLY ONE GROUP LEFT FOR r = ",r," and k = ", k) )
    return (pval=99,ok,patiF,nFactors)
  }

  # components
  SurFunction <- Surv(time = Overall.Survival..Months.,
                      event = Overall.Survival.Status)~as.factor(group)
  coxFit_NBS <- coxph(SurFunction,
                      x=T, y=T, model=T, method = "breslow",
                      data=patiF)
  #summary(coxFit_NBS)
  mfit <- survfit(SurFunction, data = patiF)
  
  tryCatch(coxph(SurFunction, x=T, y=T, model=T, method = "breslow",
                 data=patiF),
           warning=function(w) {
             cat("problem values: r = ",r,", k = ",k,"\n")
             ok <- FALSE
           }
  )
  
  s <- summary(coxFit_NBS)
  
  if(plots) {
    ####
    # #complete plot
    # km4 <- prodlim(Surv(time = Overall.Survival..Months., event = Overall.Survival.Status)~
    #                 factor(group),
    #               data=patiF)
    # plot(km4, logrank=TRUE)
    
   
    #Plot for paper
    pdf( file = plotName,  onefile = TRUE, width = 9, height = 7)
    # #tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
    plot.Survival4paper(coxFit_NBS, mfit, location = location,
                        colorsP = c25[1:nFactors],
                        colorsL = c25[1:nFactors],
                        labelClu = paste("subtype ", seq(nFactors), " (", tr, ")", sep='') ,
                        font_size_times = 1.8, legendbg = "gray98")
    dev.off()
  }
  
  return (list(pval = s$sctest[3], ok,patiF,nFactors)) 
}
