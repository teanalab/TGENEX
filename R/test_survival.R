# run survival analisis on datasets

library(plyr)
library(survival)
library(classInt) #for sturges
source("R/6-plotSurvivalUpdate.R") #for survival plot

# parameters ##########
k=3
minimumGroupSize = 4
foldermRNA="R/clean_mRNA/"
outFolder = "kmplots/"
#1. get list of diseases
diseases = list.files(path=foldermRNA)
diseases = gsub(".RData", "", diseases, fixed = TRUE)


plotIdealSurvival <- function(disease){
  #load survival data
  load( paste0("survival/",disease,"/survival.Rd",sep='') )

  #load mRNA data
  load(paste0(foldermRNA, "/", disease,".RData"))
  eval(parse(text= paste0("TCGA_DB <- ",disease,".mRNA")))

  # find intersection between survival and mRNA
  patients <- row.names(TCGA_DB)
  p_s <-survivalData$bcr_patient_barcode
  patients <- intersect(patients, p_s)

  # update tables
  survivalData <- survivalData[patients,]


  # get survival groups from time
  aInt <- classIntervals(survivalData$times, n=k, style = "jenks")
  aCategory <- findCols(aInt)
  survivalData$group <- aCategory

  # survival analysis
  cat("\n Disease", disease)
  cat(". Groups", table(survivalData$group))
  # components
  SurFunction <- Surv(time = times,
                      event = patient.vital_status)~as.factor(group)
  coxFit <- coxph(SurFunction,
                  x=T, y=T, model=T, method = "breslow",
                  data=survivalData)
  cat(". Log p-val",format.pval(summary(coxFit)$sctest[3], digits = 3) )
  mfit <- survfit(SurFunction, data = survivalData)


  # FINALLY, plot
  tr = table(survivalData$group)
  fileName = paste0( outFolder,disease,"_tr.pdf" )
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "topright",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()

  #crop borders to fit paper
  #system( paste("pdfcrop", fileName, fileName) )

  fileName = paste0( outFolder,disease,"_br.pdf" )
  #Plot for paper
  pdf( file = fileName,  onefile = TRUE, width = 9, height = 7)
  plot.Survival4paper(coxFit, mfit, location = "bottomright",
                      colorsP = c25[1:k],
                      colorsL = c25[1:k],
                      labelClu = paste("subtype ", seq(k), " (", tr, ")", sep='') ,
                      font_size_times = 1.8, legendbg = "gray98")
  dev.off()

  #crop borders to fit paper
  #system( paste("pdfcrop", fileName, fileName) )

  save(survivalData, file = paste0(outFolder,disease,"_survival.RData"))

}



#Start trying with one disease
for (i in seq_along(diseases))
{
  disease = diseases[i]
  plotIdealSurvival(disease)
}





######## OUTPUT
# Disease BRCA. Groups 749 280 69. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/BRCA.pdf'.
#
# Disease COAD. Groups 396 44 17. Log p-val 2.41e-11PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/COAD.pdf'.
#
# Disease COADREAD. Groups 552 56 20. Log p-val 1.7e-13PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/COADREAD.pdf'.
#
# Disease GBMLGG. Groups 892 175 43. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/GBMLGG.pdf'.
#
# Disease KIPAN. Groups 487 306 144. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/KIPAN.pdf'.
#
# Disease KIRC. Groups 233 184 120. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/KIRC.pdf'.
#
# Disease KIRP. Groups 188 69 31. Log p-val 4.9e-14PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/KIRP.pdf'.
#
# Disease LGG. Groups 374 111 30. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/LGG.pdf'.
#
# Disease LUAD. Groups 388 107 8. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/LUAD.pdf'.
#
# Disease LUSC. Groups 329 123 42. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/LUSC.pdf'.
#
# Disease OV. Groups 279 221 76. Log p-val <2e-16PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/OV.pdf'.
#
# Disease READ. Groups 156 13 2. Log p-val 0.012PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/READ.pdf'.
#
# Disease UCEC. Groups 328 165 54. Log p-val 4.16e-15PDFCROP 1.38, 2012/11/02 - Copyright (c) 2002-2012 by Heiko Oberdiek.
# ==> 1 page written on `kmplots/UCEC.pdf'.
