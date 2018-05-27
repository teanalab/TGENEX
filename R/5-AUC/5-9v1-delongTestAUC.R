#' version 1 (3/14/2018)
#' Testing the AUCs with the DeLong test from the pROC package



rm(list = ls())
load("temp/109-Data.RData")
loadls("pROC",FALSE)

#' Get all needed data
LoadMyData <- funtion()
{
  load(file="results/10March2018/103-AUC_NMF_i9.RData")

  outcometr1 <- Surv.rsp[,2]
  predicttr1 <- lp
  outcomete1 <- Surv.rsp.new[,2]
  predictte1 <- lpnew

  rm(list=ls()[-which(ls() %in% c("outcometr1","outcomete1",
                                  "predicttr1", "predictte1" ))])

  load(file="results/10March2018/104-AUC_NMF_i10.RData")

  outcometr2 <- Surv.rsp[,2]
  predicttr2 <- lp
  outcomete2 <- Surv.rsp.new[,2]
  predictte2 <- lpnew

  rm(list=ls()[-which(ls() %in% c("outcometr1","outcomete1",
                                  "predicttr1", "predictte1" ,
                                  "outcometr2", "outcomete2",
                                  "predicttr2", "predictte2" ))])

  load(file="temp/104-AUC_NMF_clini_i6.RData")

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
  onlyBaseLibs()
  loadls("pROC")

  save.image(file="temp/109-Data.RData")
}


#' Title
#'
#' @section step 1
#'
#' @examples
foo <- function() {
  auc1 <- auc(outcometr1, predicttr1) #NFT
  auc2 <- auc(outcometr2, predicttr2) #NMF GENES

  roc.test(auc1,auc2,method="delong")

  auc3 <- auc(outcomete3, predictte3) #NMF CLINI
  roc.test(auc3,auc1,method="delong")
}

#end
sess <- sessionInfo() #save session on variable
save.image(file="temp/x.RData")
