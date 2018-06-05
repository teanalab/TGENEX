#version 1 (Jun/5/2018)

rm(list = ls())
load("temp/6-3v1-Top5_NMF_mut_k7.RData")
loadls("plyr survival missForest survAUC prodlim survminer", F)

#standard plot
ggsurvplot(mfitTop, conf.int = F,pval = test$sctest[3], risk.table = T )

table(surCli[which(as.numeric(surCli$component) %in% bc4),"component"])

setEPS()
postscript("output4paper/KM_top5_NMF_mut.eps", width = 12, height = 9)
plot.Survival4paper(coxFitTop, mfitTop, location = "bottomleft",
                    colorsP = c("red","black","blue","#1E8449", "orange"),
                    colorsL = c("blue","#1E8449","red","black", "orange"),
                    labelClu = c("subtype 1 (4)","subtype 2 (302)",
                                 "subtype 3 (4)", "subtype 4 (10)", "subtype 5 (7)"),
                    font_size_times = 0.5)
dev.off()
