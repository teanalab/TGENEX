#version 1 (Jun/5/2018)
rm(list = ls())

load("temp/6-3v0-Top5_NTF_k7.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)
source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")

#standard plot
ggsurvplot(mfitTop, conf.int = F,pval = test$sctest[3], risk.table = T )
table(surCli[which(as.numeric(surCli$component) %in% bc4),"component"])

#plot for paper
setEPS()
postscript("output4paper/KM_top5_NTF.eps", width = 12, height = 9)
plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
                    colorsP = c("red","black","blue","#1E8449", "orange"),
                    colorsL = c("blue","#1E8449","red","black", "orange"),
                    labelClu = c("subtype 1 (26)","subtype 2 (9)",
                                 "subtype 3 (9)", "subtype 4 (137)",
                                 "subtype 5 (159)") )
dev.off()

