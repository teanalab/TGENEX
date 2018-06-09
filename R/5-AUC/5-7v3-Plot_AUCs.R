# version 3 (May/30/2018)
# Plot AUCs comparisons

rm(list = ls())

# LoadMyData <- function()
# {
  load(file="output4paper/random_AUC.RData")
  load(file="output4paper/NTF_AUC.RData")
  load(file="output4paper/NMF_PxC_AUC.RData")
  load(file="output4paper/NMF_PxM_AUC.RData")
  load("data/myLib.RData")
  loadlib("ggplot2",F)

  kMax=length(NTF_AUC)
  interval = c(2:kMax)

  x = as.data.frame(NTF_AUC[interval])
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NTF_AUC_ = data.frame(k=interval, auc= meanX, method = "M1", sd= sdX)

  x = as.data.frame(NMF_PxM_AUC[interval])
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NMF_PxM_AUC_ = data.frame(k=interval, auc= meanX, method = "M2", sd= sdX)

  x = as.data.frame(NMF_PxC_AUC[interval])
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  NMF_PxC_AUC_ = data.frame(k=interval, auc= meanX, method = "M3", sd= sdX)

  x = as.data.frame(random_AUC[interval])
  meanX <- apply(x,2,mean)
  sdX <- apply(x,2,sd)
  randomAUC_ = data.frame(k=interval, auc= meanX, method = "M4", sd= sdX)
#}


#WithRandomAsBaseline <- function(){
  data10 <- rbind(NTF_AUC_,NMF_PxM_AUC_, NMF_PxC_AUC_, randomAUC_)

  ggsave("output4paper/AUC.eps", width = 8, units = "in")
  theme_set(theme_light(base_size = 18))
  g1 <- ggplot(data10, aes(x=k, y=auc, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    # xlab("number of breast cancer subtypes") +
    # xlab("number of sub-types") +
    ylab("AUC") +
    geom_point(size=3) +
    # Remove title for all legends
    theme(legend.title=element_blank())
  g1 + scale_x_continuous(name="number of sub-types", breaks = (2:kMax))
  #http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels

  ##Export png 700 x 380  for paper as Figure7.png
  #Export PDF 8 x 5 inches for email
#}

  which.max(NTF_AUC_$auc)
  max(NTF_AUC_$auc)
