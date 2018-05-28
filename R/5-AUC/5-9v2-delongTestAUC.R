#' version 2 (5/27/2018)
#' Testing the AUCs with the DeLong test from the pROC package

rm(list = ls())
load("temp/5-9v2-data.RData")
loadls("pROC",FALSE)

LoadMyData <- funtion()
{
  # get max AUC from NTF
  load(file="temp/2v4-AUC_NTF_k1.RData")

  outcometr1 <- Surv.rsp[,2]
  predicttr1 <- lp
  outcomete1 <- Surv.rsp.new[,2]
  predictte1 <- lpnew

  rm(list=ls()[-which(ls() %in% c("outcometr1","outcomete1",
                                  "predicttr1", "predictte1" ))])

  load(file="temp/3v4-AUC_NMF_PxM_k10.RData")

  outcometr2 <- Surv.rsp[,2]
  predicttr2 <- lp
  outcomete2 <- Surv.rsp.new[,2]
  predictte2 <- lpnew

  rm(list=ls()[-which(ls() %in% c("outcometr1","outcomete1",
                                  "predicttr1", "predictte1" ,
                                  "outcometr2", "outcomete2",
                                  "predicttr2", "predictte2" ))])

  load(file="temp/4v2-AUC_NMF_PxC_k9.RData")

  outcometr3 <- Surv.rsp[,2]
  predicttr3 <- lp
  outcomete3 <- Surv.rsp.new[,2]
  predictte3 <- lpnew

  rm(list=ls()[-which(ls() %in% c("outcometr1","outcomete1",
                                  "predicttr1", "predictte1" ,
                                  "outcometr2", "outcomete2",
                                  "predicttr2", "predictte2" ,
                                  "outcometr3", "outcomete3",
                                  "predicttr3", "predictte3" ))])

  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  loadls("pROC",F)

  save.image(file="temp/5-9v2-data.RData")
}


#' computeDeLongT
computeDeLongT <- function() {
  auc1 <- auc(outcometr1, predicttr1) #NFT
  auc2 <- auc(outcometr2, predicttr2) #NMF GENES

  roc.test(auc1,auc2,method="delong")

  auc1 <- auc(outcometr1, predicttr1) #NFT
  auc3 <- auc(outcomete3, predictte3) #NMF CLINI
  roc.test(auc1,auc3,method="delong")
}

#end
sess <- sessionInfo() #save session on variable
save.image(file="temp/5-9v2.RData")
