#' version 3 (May/28/2018)
#' Cstatistic
#' matrix clinical versus patients
#' AUC comparing NMF. macro-average
#' normalize rows of patien Factor
#' 80-20

rm(list = ls())
load(file="temp/5-4v2-data.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  load(file="data/NMF_30_PxC.RData")
  load("data/survivalClinical-5-1_30.RData")

  rm(list=ls()[-which(ls() %in% c("patientNMF_PxC","survivalClinical",
                                  "patients","kMax", "numberOfPatiens"))])

  #load my libraries
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  #required for the cox proportional hazard model
  loadls("plyr survival Rcpp missForest survAUC perry",F)

  save.image(file="temp/5-4v3-data.RData")
}

#line by line
NMF_PxM_AUC <- function()
{
  numSplits = 10
  AUC_CD_all <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  smp_size <- floor(0.2 * numberOfPatiens) #testing sample size

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in c(2:kMax)){
    load(file= paste("data/NMF/NMF_PxC_R_",k,".RData", sep='') )
    #normalize rows
    patiF <- get(paste("patientNMF_PxC_R_",k,sep=''))
    patiF <- data.frame( t(apply(patiF,1,function(x){x/sum(x)})) )
    names(patiF) <- paste('V',seq(1:k),sep='')
    patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
    patiF <- patiF[,-1]

    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)

    set.seed(k*1234)  # set seed for reproducibility
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
      test_ind <- obsInTest$subsets[,i]
      train <- patiF[-test_ind, ]
      test <- patiF[test_ind, ]

      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      model=T, x=TRUE, y=TRUE, data=train)
      #AUC
      lp <- predict(coxFit)
      lpnew <- predict(coxFit, newdata=test)
      Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                       event = train$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                           event = test$Overall.Survival.Status)
      Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      Cstat_K[i] <- Cstat
    }
    save.image(paste("temp/4v4-AUC_NMF_PxC_k",k,".RData",sep=''))
    AUC_CD_all[k] = list(AUC_CD_K)
    Cstat_allK[k] <- list(Cstat_K)
  }

  NMF_PxC_AUC <- AUC_CD_all
  NMF_PxC_Cstat <- Cstat_allK
  save(NMF_PxC_AUC, file = "output4paper/NMF_PxC_AUC.RData")
  save(NMF_PxC_Cstat, file = "output4paper/NMF_PxC_Cstat.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-4v4.RData")
