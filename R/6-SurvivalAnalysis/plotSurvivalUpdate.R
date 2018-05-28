#To move to github


#' Plot complex Survival Curves
#'
#' @param groups
#' @param mainTitle
#' @param survData
#' @param location
#' @param colorsP
#' @param colorsL
#' @param labelClu
#' @param centerT
#'
#' @return
#' @export
#'
#' @examples
#'   coxFit <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
# Diagnosis.Age+#ER.Status+
#   HER2.Status+Metastasis.Coded+#PR.Status+
#   Tumor..T1.Coded,
# x=TRUE, y=TRUE, method = "breslow", data=survivalClinicalNTF)
#'   mfit <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ #Converted.Stage+
# Diagnosis.Age+#ER.Status+
#   HER2.Status+Metastasis.Coded+#PR.Status+
#   Tumor..T1.Coded,
# data=survivalClinicalNTF)
#' plot.SurvivalComplex (coxFit, mfit, "clinical variables")
plot.SurvivalComplex <- function(coxFit, mfit, mainTitle)
{
  require("survival")
  # location
  # colorsP = c("red","black","blue")
  # ,
  # colorsL = c("red","black","blue")
  # ,
  # labelClu = c("group1","group2", "group3")
  centerT = 0.4


  font_size_times = 1 #1.5
  par( mgp=c(1.8,0.3,0), #axis title, axis labels and axis line
       mai=c(2.0,2.0,1.0,1.0) ) #c(bottom, left, top, right)

  # coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups,
  #                 data = survData, ties="exact")
  # print(summary(coxFit))

  resCox <- mySumm(coxFit)

  pVal<-format.pval(resCox$sctest[3], digits = 3)

  plot(mfit, col=rainbow(resCox$sctest[2]),
       #lty = 1, #line type solid
       #lwd=5, #line tightness
       #main = mainTitle,
       cex.lab=font_size_times, #font size
       #xaxt="n", #do not plot x axe
       xlab = "Years",
       conf.int=F,
       #cex.main=(font_size_times+0.1),
       cex.sub =font_size_times,
       cex.axis =font_size_times,
       mark=3, #type of sensored  mark http://www.statmethods.net/advgraphs/parameters.html
       #cex = 2, #size of mark
       #xscale = 365.25, #values in years
       ylab="Survival probability",
       xaxs="i", #survival curve should touch the y-axis
       tck=0.01 #length of tick marks
  )
  title(mainTitle, line = 0.5,
        adj=centerT, #centers the titles (0 left - 1 right)
        cex.main=(font_size_times+0.3))


  #   axis(1, at = seq(0, max(survData$Survival), length.out = 15),
  #       labels=round(seq(0, 4.1, length.out = 15), digits=1), las=1) #customize x axe
  # legend(location,
  #        labelClu,
  #        col=colorsL,
  #        text.col=colorsL,
  #        text.font=2, #bold
  #        fill=colorsL, #boxes
  #        cex=font_size_times,
  #        title = paste(" Cox p-value = ", pVal),
  #        title.col="black")
}



#http://www.r-statistics.com/2013/07/creating-good-looking-survival-curves-the-ggsurv-function/
todo111 <-function()
{
  loadls("pec")
  data(GBSG2,package="pec")
  GBSG2$tgrade <- as.factor(as.character(GBSG2$tgrade))
  GBSG2$Age <- cut(GBSG2$age,c(0,40,60,81),labels=c("21-40","41-60","61-80"))

  loadlib("pec", bioc=FALSE)
  library(prodlim)
  library(survival)


  #http://staff.pubhealth.ku.dk/~tag/Teaching/share/R-tutorials/SurvivalAnalysis.html
  #library(devtools)
  ## with the command
  ## install_github("TagTeam/Publish")

  library(Publish)
  #WHAT IS THIS?
  quantile(prodlim(Hist(time,cens)~1,data=GBSG2,reverse=TRUE))



  km0 <- prodlim(Hist(time,cens)~1,data=GBSG2)
  plot(km0)

  km1 <- prodlim(Hist(time,cens)~tgrade,data=GBSG2)
  par(mar=c(7,7,5,5), mgp=c(3,1,0))
  plot(km1,
       # atrisk.labels=paste("Tumor grade: ",c("I","II","III"),": "),
       #atrisk.title="",

       xlab="Years",  # label for x-axis
       axis1.at=seq(0,2900,365.25), # time grid for x-axis
       axis1.labels=0:7, # time labels for x-axis
       legend.x="bottomleft", # positition of legend
       legend.cex=0.8, # font size of legend
       legend.title="Tumor Grade\n", #
       logrank=TRUE) # show log-rank p-value




  km0 <- prodlim(Hist(time,cens)~1,data=GBSG2)
  quantile(km0)

  library(pec)
  cox2 <- coxph(Surv(time,cens)~tgrade+age+tsize+pnodes,data=GBSG2)
  newdata <- data.frame(tgrade=c("I","II","III"),age=50,tsize=30,pnodes=8)
  plotPredictSurvProb(cox2,
                      sort(unique(GBSG2$time)),
                      newdata=newdata,
                      col=c("darkgreen","darkorange","red"),
                      legend.title="Tumor grade",
                      legend.legend=c("I","II","III"))
  mtext("Individualized survival curves from multiple Cox regression
        age=50, tumor size= 30, no. positive lymph nodes=8",line=1.5,cex=1.3)



  cox3 <- coxph(Surv(time,cens)~tgrade+age+strata(menostat)+tsize+pnodes+progrec+estrec,data=GBSG2)
  summary(cox3)


  par(cex=4)
  cox2 <- coxph(Surv(time,cens)~tgrade+age+tsize+pnodes,data=GBSG2)
  rt2 <- do.call("rbind",regressionTable(cox2))
  plotConfidence(x=rt2$HazardRatio,
                 lower=rt2$Lower,
                 upper=rt2$Upper,
                 xlim=c(0,5),
                 factor.reference.pos=1,
                 labels=c("Tumor grade I",
                          "                     II",
                          "                     III",
                          "Age (years)",
                          "Tumor size (mm)",
                          "No. positive\nlymph nodes"),
                 title.labels="Factor",
                 cex=2,
                 xlab.cex=1.3,
                 xlab="Hazard ratio")
}



#plot.survfit manual
#https://stat.ethz.ch/R-manual/R-devel/library/survival/html/plot.survfit.html
# @param groups is a factor
# @param mainTitle is a string
# @param survivalData is data.frame with headers "PatientID" "Survival"  "Death"
#> head(survivalGBM)
#                     PatientID Survival Death
#1 TCGA-02-0001-01C-01T-0179-07      358     1
#2 TCGA-02-0003-01A-01T-0179-07      144     1
#3 TCGA-02-0007-01A-01T-0179-07      705     1
#
# Example 1:
# pdf( file=paste(resultsFolder, "SNF_all_kaplan.pdf", sep="")  ,  onefile=TRUE)
#      plot.SurvivalK3 (groups, mainTitle, survivalData, "bottomright")
# dev.off()
#
#
# Example 2:
# mainTitle <- "b) Survival curve, SNF, selected genes"
# pdf(
#   file = paste(resultsFolder, "SNF_union_kaplan.pdf", sep = ""),  onefile =
#     TRUE, width = 9, height = 7
# )
# plot.SurvivalK3 (groups, mainTitle, survivalData,"bottomright",
#                  colorsP = c("red","blue","black"))
# dev.off()
#
#
# Example 3:
# labelsGroups <- c(paste("group 1 (",table(groups)[2], ")", sep="" ),
#                   paste("group 2 (",table(groups)[1], ")", sep="" ),
#                   paste("group 3 (",table(groups)[3], ")", sep="" ) )
# colorsLines <- c("black","red","blue")
# colorsLabels <- c("red","black","blue")
#  pdf( file = pdfName,  onefile = TRUE, width = 9, height = 7)
#  plot.SurvivalK3 (groups, mainTitle, survivalData,"bottomright",
#                   labelClu = labelsGroups,
#                   colorsP = colorsLines,
#                   colorsL = colorsLabels,
#                   centerT = centerT )
#  dev.off()
#
#
#TRY THIS: http://www.r-statistics.com/2013/07/creating-good-looking-survival-curves-the-ggsurv-function/
#
plot.SurvivalKN <- function(groups, mainTitle, survData = survivalData,
                            location, colorsP = c("red","black","blue"),
                            colorsL = c("red","black","blue"),
                            labelClu = c("group1","group2", "group3"),
                            centerT = 0.4){

  font_size_times = 1 #1.5
  require("survival")
  par( mgp=c(1.8,0.3,0), #axis title, axis labels and axis line
       mai=c(2.0,2.0,1.0,1.0) ) #c(bottom, left, top, right)

  coxFit <- coxph(Surv(time = Survival, event = Death) ~ groups,
                  data = survData, ties="exact")
  print(summary(coxFit))
  resCox <- mySumm(coxFit)
  mfit <- survfit(Surv(time = Survival, event = Death) ~ groups, data = survData)
  pVal<-format.pval(resCox$sctest[3], digits = 3)
  plot(mfit, col=rainbow(12),
       lty = 1, #line type solid
       lwd=5, #line tightness
       #main = mainTitle,
       cex.lab=font_size_times, #font size
       #xaxt="n", #do not plot x axe
       xlab = "Years",
       conf.int=F,
       #cex.main=(font_size_times+0.1),
       cex.sub =font_size_times,
       cex.axis =font_size_times,
       mark=3, #type of sensored  mark http://www.statmethods.net/advgraphs/parameters.html
       cex = 2, #size of mark
       xscale = 365.25, #values in years
       ylab="Survival probability",
       xaxs="i", #survival curve should touch the y-axis
       tck=0.01 #length of tick marks
  )
  title(mainTitle, line = 0.5,
        adj=centerT, #centers the titles (0 left - 1 right)
        cex.main=(font_size_times+0.3))


  #   axis(1, at = seq(0, max(survData$Survival), length.out = 15),
  #       labels=round(seq(0, 4.1, length.out = 15), digits=1), las=1) #customize x axe
  legend(location,
         labelClu,
         col=colorsL,
         text.col=colorsL,
         text.font=2, #bold
         fill=colorsL, #boxes
         cex=font_size_times,
         title = paste(" Cox p-value = ", pVal),
         title.col="black")
}

#@ example:
# plotCompleteKplan(groups, survivalData, resultsFolder=getwd(), pdfName1="/pdf1.pdf")
plotCompleteKplan <- function (groups, survivalData, resultsFolder=getwd(), pdfName1="/pdf1.pdf"){

  #to check sizes
  library(prodlim)
  km0 <- prodlim(Hist(time = Survival, event = Death) ~ groups,
                 data = survivalData)

pdf( file = paste(resultsFolder, pdfName1 , sep = ""),  onefile = TRUE,
      width = 9, height = 7)

  plot(km0
       ,
        atrisk.labels=paste("Tumor grade: ",c("I","II","III"),": "),
       atrisk.title="",

       xlab="Years",  # label for x-axis
       axis1.at=seq(0,2900,365.25), # time grid for x-axis
       axis1.labels=0:7, # time labels for x-axis
       legend.x="bottomleft", # positition of legend
       legend.cex=0.8, # font size of legend
       legend.title="\n", #
       logrank=FALSE
    ) # show log-rank p-value

  dev.off()
}







#'@example:
#'   coxFit <- coxph(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ Diagnosis.Age,
#'   data=survivalClinicalNTF)
#'   mfit <- survfit(Surv(time = Overall.Survival..Months., event = (Overall.Survival.Status=="DECEASED"))~ Diagnosis.Age,
#'   data=survivalClinicalNTF)
#'   pdf( file = pdfName,  onefile = TRUE, width = 9, height = 7)
#'  plot.Survival4paper (coxFit, mfit, labelClu = NULL)
#'  dev.off()
plot.Survival4paper <- function(coxFit, mfit, mainTitle=NULL,
                                location = "right",
                                colorsP = c("red","black","blue"),
                                colorsL = c("red","black","blue"),
                                labelClu = c("group1","group2", "group3"),
                                centerT = 0.4,
                                legendTitle = "p-value =",
                                font_size_times = 1.3)
{
  require("survival")
  par( mgp=c(1.8,0.3,0), #axis title, axis labels and axis line
       mai=c(2.0,2.0,1.0,1.0) ) #c(bottom, left, top, right)
  print(summary(coxFit))
  resCox <- summary(coxFit)
  pVal<-format.pval(resCox$sctest[3], digits = 3)
  plot(mfit,
       col=colorsP,
       lty = 1, #line type solid
       lwd=5, #line tightness
       #main = mainTitle,
       cex.lab=font_size_times, #font size
       #xaxt="n", #do not plot x axe
       xlab = "Years",
       conf.int=F,
       #cex.main=(font_size_times+0.1),
       cex.sub =font_size_times,
       cex.axis =font_size_times,
       mark=3, #type of sensored  mark http://www.statmethods.net/advgraphs/parameters.html
       cex = 2, #size of mark
       xscale = 12, #365.25, #values in years
       ylab="Survival probability",
       xaxs="i", #survival curve should touch the y-axis
       tck=0.01 #length of tick marks
  )

  if(!is.null(mainTitle))
  {
    title(mainTitle, line = 0.5,
        adj=centerT, #centers the titles (0 left - 1 right)
        cex.main=(font_size_times+0.3))
  }

  if(is.null(labelClu))
  {
    legend(location,
           legend="",
           text.font=2, #bold
           cex=font_size_times,
           title = paste(legendTitle, pVal) ,
           title.col="black")
  } else {
  #   axis(1, at = seq(0, max(survData$Survival), length.out = 15),
  #       labels=round(seq(0, 4.1, length.out = 15), digits=1), las=1) #customize x axe
  legend(location,
         labelClu,
         col=colorsL,
         text.col=colorsL,
         text.font=2, #bold
         fill=colorsL, #boxes
         cex=font_size_times,
         title = paste(legendTitle, pVal),
         title.col="black")
  }
}
