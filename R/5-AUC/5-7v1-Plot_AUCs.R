#version 1 (4/1/2018)
# Plot 2 AUCs comparisons

# If I show AUC SH I can show how inestable NMF is and how it depends of k

rm(list = ls())
load("temp/107-data.RData")
loadlib("ggplot2")

#Get all needed data
LoadMyData <- function()
{
  load(file="temp/5-1-AUC-random.RData")
  load(file="temp/103-AUC-NTF.RData")
  load(file="temp/104-AUC-NMF.RData")
  load(file="temp/105-NMF_lasso.RData")
  load(file="temp/106-NTF_lasso.RData")

  rm(list=ls()[-which(ls() %in% c("AUC_CD_all", "AUCvsK_103_NTF", "AUCvsK_104_NMF", "AUCvsK_105_NMF", "AUCvsK_106_NTF"))])

  load(file="temp/AUCvsK_104_NMF_clini.RData")

  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  sapply(libs, function(u) {source(u)})

  loadlib("ggplot2")


  save.image(file="temp/107-data.RData")
}


WithRandomAsBaseline <- function(){

  AUCvsK_103_NTF$method = "M2"
  AUCvsK_104_NMF$method = "M1"
  AUCvsK_104_NMF_clini$method = "M3"

  #Only first points
  data10 <- rbind(AUCvsK_103_NTF[2:10,],AUCvsK_104_NMF[2:10,], AUCvsK_104_NMF_clini[2:10,])
  data10 = data10[,-3]  #remove the third column
  data10 <- rbind(data10, data.frame(k=seq(2,10), auc_CD=AUC_CD_all[2:10], method="M4") )

  #With clinical data
  theme_set(theme_light(base_size = 18))
  g1 <- ggplot(data10, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    xlab("# of most prevalent subtypes") +
    ylab("AUC") +
    geom_point(size=5) +
    # Remove title for all legends
    theme(legend.title=element_blank())

g1


  ##Export PDF 700 x 380  for paper as Figure7.png
  #Export PDF 8 x 5 inches for email
}



#not included because the clinical data model is not including all the variables
FourLinesClinicalCoxAsBaseline <- function(){

  AUCvsK_103_NTF$method = "M2"
  AUCvsK_104_NMF$method = "M1"
  AUCvsK_104_NMF_clini$method = "M3"

  #Only first points
  data10 <- rbind(AUCvsK_103_NTF[2:10,],AUCvsK_104_NMF[2:10,], AUCvsK_104_NMF_clini[2:10,])

  #With clinical data
  theme_set(theme_light(base_size = 18))
  g1 <- ggplot(data10, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    xlab("# of most prevalent subtypes") +
    ylab("AUC") +
    geom_point(size=5) +
    # Remove title for all legends
    theme(legend.title=element_blank())

  clinicalAUC_baseline = AUC_CD_iauc-0.14
  g2 = g1 + geom_hline(linetype = 1, yintercept=clinicalAUC_baseline,
                       colour="#5086ff")

  data4Line=data.frame(X=c(2:10), Y=clinicalAUC_baseline, method="M4")
  g2 + geom_line(linetype = 1, data=data4Line, aes(x=X, y=Y),
                 lineend = "square", linejoin = "bevel",
                 linemitre = 1 )


  ##Export PDF 700 x 380  for paper as Figure7.png
  #Export PDF 8 x 5 inches for email
}


ThreeLines <- function(){
  clinicalAUC_baseline = AUC_CD_iauc-0.14

  ##1- Using Cox regression

  AUCvsK_103_NTF$method = "N Matrix F"
  AUCvsK_104_NMF$method = "N Tensor F"

  dataALL <- rbind(AUCvsK_103_NTF,AUCvsK_104_NMF)
  dataALL <- dataALL[-16,]
  dataALL <- dataALL[-1,]

  ggplot(dataALL, aes(x=k, y=auc_SH, colour=method, shape=method, fill=method)) +
    ggtitle("Cox regression and AUC SH") +
    geom_line(linetype="dashed") +
    geom_point(size=3)


  ggplot(dataALL, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    ggtitle("higher AUC at different number of components (k)") +
    geom_line(linetype="dashed") +
    xlab("k first components") +
    ylab("A U C") +
    geom_point(size=3)

  #Export PDF 8 x 5 inches
  # AUCvsK_103_NTF$method = "NMF"
  # AUCvsK_104_NMF$method = "NTF"


   AUCvsK_103_NTF$method = "M2"
   AUCvsK_104_NMF$method = "M1"


  #Only first points
  data10 <- rbind(AUCvsK_103_NTF[2:10,],AUCvsK_104_NMF[2:10,])
  ggplot(data10, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    #xlab("# of top k components") +
    xlab("k components") +
    ylab("AUC") +
    geom_point(size=3) +
    # Remove title for all legends
    theme(legend.title=element_blank())
  #Export PDF 8 x 5 inches for email


  #Only 2 AUCs for paper
  theme_set(theme_light(base_size =18))
  ggplot(data10, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    #xlab("# of components") +
    xlab("# of top k components") +
    ylab("AUC") +
    geom_point(size=5) +
    # Remove title for all legends
    theme(legend.title=element_blank())
  # with out clinical ^

  #With clinical data
  theme_set(theme_light(base_size =18))
  g1 <- ggplot(data10, aes(x=k, y=auc_CD, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    xlab("# of most prevalent subtypes") +
    ylab("AUC") +
    geom_point(size=5) +
    # Remove title for all legends
    theme(legend.title=element_blank())
  #THE ACTUAL PLOT
  g2 = g1 + geom_hline(linetype = 1, yintercept=clinicalAUC_baseline,
               colour="#5086ff")

  ##Export PDF 7 x 4 inches for paper as 107-AUC.pdf


  #To manually edit the legend
  data4Line=data.frame(X=c(2:10), Y=clinicalAUC_baseline, method="M3")
  g2 + geom_line(linetype = 1, data=data4Line, aes(x=X, y=Y),
                 lineend = "square", linejoin = "bevel",
                 linemitre = 1 )

  ##########
  }

####
# if I ever need to put the actual value
#g1 +  geom_text(x=3, y=clinicalAUC_baseline+0.02, label=clinicalAUC_baseline,
#            colour="darkgreen")
#nicer box text at http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software


##1- Using (L1) Lasso regression
lassoR <- function(){
  AUCvsK_106_NTF$method = "N Matrix F"
  AUCvsK_105_NMF$method = "N Tensor F"

  dataALL <- rbind(AUCvsK_106_NTF,AUCvsK_105_NMF)

  dataALL <- dataALL[-16,]
  dataALL <- dataALL[-1,]

  ggplot(dataALL, aes(x=k, y=auc, colour=method, shape=method, fill=method)) +
    ggtitle("               (L1) Lasso regression") +
    geom_line(linetype="dashed") +
    geom_point(size=3)
}



