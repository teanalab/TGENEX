# version 4 (May/26/2018)
# AUC comparing NMF. macro-average
# normalize rows of patien Factor
# 80-20

rm(list = ls())
load(file="temp/5-3v4-data.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  load(file="data/NMF_PxM.RData")
  load("data/survivalClinical-5-1.RData")

  rm(list=ls()[-which(ls() %in% c("patientNMF_PxM","survivalClinical",
                                  "patients","kMax", "numberOfPatiens"))])

  #load my libraries
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  #required for the cox proportional hazard model
  loadls("plyr survival Rcpp missForest survAUC perry",FALSE)

  save.image(file="temp/5-3v4-data.RData")
}

#line by line
NMF_PxM_AUC <- function()
{
  numSplits = 10
  AUC_CD_all <- rep(0,kMax)


  #normalize rows
  patiF <- data.frame( t(apply(patientNMF_PxM,1,function(x){x/sum(x)})) )
  names(patiF) <- paste('V',seq(1:10),sep='')
  patiF <- cbind.data.frame(patiF,survivalClinical)

  for (k in seq(2,kMax)){
    AUC_CD_K <- rep(0,numSplits)

    #we used 10 random splits
    # 80% of the sample size for training
    set.seed(k*1234)  # set seed for reproducibility
    #testing sample
    smp_size <- floor(0.2 * numberOfPatiens)

    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
      test_ind <- obsInTest[[4]][,i]
      train <- patiF[-test_ind, ]
      test <- patiF[test_ind, ]

      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~V1+V2+V3+V4+V5+V6+V7+V8+V9+V10,
                      x=TRUE, y=TRUE, data=train)

      #AUC
      lp <- predict(coxFit)
      lpnew <- predict(coxFit, newdata=test)
      Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                       event = train$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                           event = test$Overall.Survival.Status)
      maxMonth <- max(survivalClinical$Overall.Survival..Months.)
      #times	The vector of time points at which AUC is evaluated.
      times <- seq(from = maxMonth/100, to = maxMonth, by = maxMonth/100)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
    }

    #macro-average
    AUC_CD_all[k] = mean(AUC_CD_K)
  }
  AUC_CD_all

  NMF_PxM_AUC <- AUC_CD_all
  save(NMF_PxM_AUC, file = "output4paper/NMF_PxM_AUC.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-3v4.RData")
