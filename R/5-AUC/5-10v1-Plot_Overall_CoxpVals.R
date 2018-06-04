# version 1 (May/31/2018)
# Plot Cox Pval overall

rm(list = ls())

#LoadMyData <- function()
#{
  load(file="output4paper/random_logRank.RData")
  load(file="output4paper/NTF_logRank.RData")
  load(file="output4paper/NMF_PxC_logRank.RData")
  load(file="output4paper/NMF_PxM_logRank.RData")
  load("data/myLib.RData")
  loadlib("ggplot2",F)

  kMax=12
  interval = c(2:kMax)

  NTF_AUC_ = data.frame(k=interval, auc= NTF_logRank[interval], method = "M1")

  NMF_PxM_AUC_ = data.frame(k=interval, auc= NMF_PxM_logRank[interval], method = "M2")

  NMF_PxC_AUC_ = data.frame(k=interval, auc= NMF_PxC_logRank[interval], method = "M3")

  randomAUC_ = data.frame(k=interval, auc= random_logRank[interval], method = "M4")
#}


#WithRandomAsBaseline <- function(){
  data10 <- rbind(NTF_AUC_,NMF_PxM_AUC_, NMF_PxC_AUC_, randomAUC_)

  #With clinical data
  theme_set(theme_light(base_size = 18))
  g1 <- ggplot(data10, aes(x=k, y=auc, colour=method, shape=method, fill=method)) +
    geom_line(linetype="dashed") +
    xlab("number of breast cancer subtypes") +
    ylab("Cox Log rank pvalue") +
    geom_point(size=3) +
    #geom_errorbar(aes(ymin=auc-sd, ymax=auc+sd), position=position_dodge(0.05))+
    # Remove title for all legends
    theme(legend.title=element_blank())
  g1

  g2 = g1 + geom_hline(linetype = 1, yintercept=0.05,
                       colour="gray")

  g2

  NTF_AUC_$auc[which(NTF_AUC_$k==7)]
#0.1176766
  NMF_PxC_AUC_$auc[which(NMF_PxC_AUC_$k==7)]
#0.04752304
  NMF_PxM_AUC_$auc[which(NMF_PxC_AUC_$k==7)]
#0.065641

  NTF_AUC_$auc[which(NTF_AUC_$k==12)]
  #0.05053501


  ##Export png 700 x 500  for paper as Figure7.png
  #Export PDF 8 x 5 inches for email
#}

