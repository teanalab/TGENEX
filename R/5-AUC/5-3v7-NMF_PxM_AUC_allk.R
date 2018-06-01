#' version 7 (May/31/2018)
#' principal verotype as covariate

rm(list = ls())
load(file="temp/5-3v6-data.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  load("data/survivalClinical-5-1_30.RData")
  rm(list=ls()[-which(ls() %in% c("survivalClinical",
                                  "patients","kMax", "numberOfPatiens"))])
  load("data/myLib.RData")
  loadls("plyr survival Rcpp missForest survAUC perry",FALSE)
  save.image(file="temp/5-3v6-data.RData")
}

#line by line
NMF_PxM_AUC <- function()
{
  numSplits = 10
  AUC_CD_all <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  LogrankP_allK_train <- rep(list(),kMax)
  LogrankP_allK <- rep(0,kMax)
  smp_size <- floor(0.2 * numberOfPatiens) #testing sample size

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in c(2:kMax)){
    load( file=paste("data/patientA_NMF_M/patientA_NMF_M_", k,".RData",sep='') )
    row.names(patientAfiliation) <- patientAfiliation$patient
    patiF <- merge(as.data.frame(patientAfiliation), as.data.frame(survivalClinical), by='row.names', all=TRUE)
    patiF <- patiF[,-1]
    patiF$component <- factor(patiF$component)



    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)
    LogrankP_K <- rep(0,numSplits)

    set.seed(k*1234)  # set seed for reproducibility
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
      test_ind <- obsInTest$subsets[,i]
      train <- patiF[-test_ind, ]
      test <- patiF[test_ind, ]

      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~component,
                      model=T, x=T, y=T, data=train)
      logp_val <- summary(coxFit)
      logp_val <- logp_val$sctest[3]

      #AUC
      lp <- predict(coxFit, type= "expected")
      lpnew <- predict(coxFit, newdata=test, type= "expected")
      Surv.rsp <- Surv(time = train$Overall.Survival..Months.,
                       event = train$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = test$Overall.Survival..Months.,
                           event = test$Overall.Survival.Status)
      #Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      #Cstat_K[i] <- Cstat
      LogrankP_K[i] <- logp_val
    }

    coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                         event = Overall.Survival.Status)~component,
                    model=T, x=T, y=T, data=patiF)
    logp_val <- summary(coxFit)
    logp_val <- logp_val$sctest[3]
    LogrankP_allK[k] <- logp_val

    save.image(paste("temp/3v6-AUC_NMF_PxM_k",k,".RData",sep=''))
    AUC_CD_all[k] = list(AUC_CD_K)
    #Cstat_allK[k] <- list(Cstat_K)
    LogrankP_allK_train[k] <- list(LogrankP_K)
  }
  NMF_PxM_AUC <- AUC_CD_all
  #NMF_PxM_Cstat <- Cstat_allK
  NMF_PxM_longrp <- LogrankP_allK
  save(NMF_PxM_AUC, file = "output4paper/NMF_PxM_AUC.RData")
  #save(NMF_PxM_Cstat, file = "output4paper/NMF_PxM_Cstat.RData")
  save(NMF_PxM_longrp, file = "output4paper/NMF_PxM_longrp.RData")

}

sess <- sessionInfo() #save session on variable
save.image("temp/5-3v6.RData")
