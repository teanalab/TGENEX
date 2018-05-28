# version 4 (May/26/2018)
# matrix clinical versus patients
# AUC comparing NMF. macro-average
# normalize rows of patien Factor
# 80-20

rm(list = ls())
load(file="temp/5-4v2-data.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  #load(file="data/NMF_PxC.RData")
  load(file="data/NMF_30_PxC.RData")
  #load("data/survivalClinical-5-1.RData")
  load("data/survivalClinical-5-1_30.RData")

  rm(list=ls()[-which(ls() %in% c("patientNMF_PxC","survivalClinical",
                                  "patients","kMax", "numberOfPatiens"))])

  #load my libraries
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  #required for the cox proportional hazard model
  loadls("plyr survival Rcpp missForest survAUC perry",F)

  save.image(file="temp/5-4v2-data.RData")
}

#line by line
NMF_PxM_AUC <- function()
{
  numSplits = 10
  AUC_CD_all <- rep(0,kMax)


  #normalize rows
  patiF <- data.frame( t(apply(patientNMF_PxC,1,function(x){x/sum(x)})) )
  names(patiF) <- paste('V',seq(1:10),sep='')
  patiF <- cbind.data.frame(patiF,survivalClinical)

  for (k in seq(kMax)){
    AUC_CD_K <- rep(0,numSplits)
    set.seed(k*1234)  # set seed for reproducibility
    smp_size <- floor(0.2 * numberOfPatiens) #testing sample. 80% of the sample size for training

    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
      test_ind <- obsInTest[[4]][,i]
      train <- patiF[-test_ind,c(1:k,kMax+1,kMax+2)]
      test <- patiF[test_ind,c(1:k,kMax+1,kMax+2)]

      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
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
    #save.image(paste("temp/4v2-AUC_NMF_PxC_k",k,".RData",sep=''))
    AUC_CD_all[k] = mean(AUC_CD_K)
  }
  AUC_CD_all

  NMF_PxC_AUC <- AUC_CD_all
  #save(NMF_PxC_AUC, file = "output4paper/NMF_PxC_AUC.RData")
}

sess <- sessionInfo() #save session on variable
#save.image("temp/5-4v2.RData")
