#version 1 (Jun/5/2018)

rm(list = ls())
load("temp/6-3v2-Top5_NMF_cli_k7.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)
source("R/6-SurvivalAnalysis/plotSurvivalUpdate.R")

#standard plot
ggsurvplot(mfitTop, conf.int = F,pval = test$sctest[3], risk.table = T )

table(surCli[which(as.numeric(surCli$component) %in% bc4),"component"])

setEPS()
postscript("output4paper/KM_top5_NMF_cli.eps", width = 12, height = 9)
plot.Survival4paper(coxFit_Afi_NMF_clin, mfitTop, location = "bottomleft",
                    colorsP = c("red","black","blue","#1E8449", "orange"),
                    colorsL = c("blue","#1E8449","red","black", "orange"),
                    labelClu = c("subtype 1 (47)","subtype 2 (59)",
                                 "subtype 3 (49)", "subtype 4 (152)",
                                 "subtype 5 (100)") )
dev.off()
