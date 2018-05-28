#' version 1 (May/28/2018)
#' Plot Cstatistics comparisons

rm(list = ls())

LoadMyData <- function()
{
  load(file="output4paper/randomCstat.RData")
  load(file="output4paper/NTF_Cstat.RData")
  load(file="output4paper/NMF_PxC_Cstat.RData")
  load(file="output4paper/NMF_PxM_Cstat.RData")

  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  loadlib("ggplot2",F)

  kMax=10
  intervalP = c(2:kMax)

  x = as.data.frame(NTF_Cstat)
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NTF_Cs = data.frame(k=intervalP, auc= meanX[intervalP], method = "M1", sd= sdX[intervalP])

  x = as.data.frame(NMF_PxM_Cstat)
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NMF_genes = data.frame(k=c(2:kMax), auc= meanX[intervalP], method = "M2", sd= sdX[intervalP])

  x = as.data.frame(NMF_PxC_Cstat)
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NMF_clini = data.frame(k=c(2:kMax), auc= meanX[2:kMax], method = "M3", sd= sdX[intervalP])

  x = as.data.frame(randomCstat)
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  randomAUC_ = data.frame(k=intervalP, auc= meanX[intervalP], method = "M4", sd= sdX[intervalP])
#}

#WithRandomAsBaseline <- function(){
  data10 <- rbind(NTF_Cs,NMF_genes, NMF_clini, randomAUC_)

  #With clinical data
  theme_set(theme_light(base_size = 18))
  g1 <- ggplot(data10, aes(x=k, y=auc, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    xlab("# of most prevalent subtypes") +
    ylab("C-statistic") +
    geom_point(size=3) +
    geom_errorbar(aes(ymin=auc-sd, ymax=auc+sd), width=2,
                  position=position_dodge(0.05))+
    theme(legend.title=element_blank()) # Remove title for all legends
  g1

  ##Export png 700 x 380  for paper as Figure7.png
  #Export PDF 8 x 5 inches for email

  # Smaller letters
  print( g1+theme_classic() )
}


