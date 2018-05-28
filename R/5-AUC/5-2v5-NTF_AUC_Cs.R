#' version 4 (May/26/2018)
#' Cstatistic
#' AUC comparing NTF. macro-average
#' normalize rows of patien Factor
#' 80-20

rm(list = ls())
load(file="temp/5-2v5-data.RData")
loadls("plyr survival Rcpp survAUC perry",F)

LoadMyData <- function(){
  load(file="data/factors.RData")
  load("data/survivalClinical-5-1_30.RData")
  # load(file="data/R_10/factors.RData")
  # load("data/R_10/survivalClinical-5-1_.RData")

  rm(list=ls()[-which(ls() %in% c("patientsF","survivalClinical","weightsC",
                                  "patients","kMax", "numberOfPatiens"))])

  #load my libraries
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/de2ecbae950d5d4b8cb1d7ba7bde5301c34186c2/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  #required for the cox proportional hazard model
  #missForest
  loadls("plyr survival Rcpp survAUC perry",F)

  save.image(file="temp/5-2v5-data.RData")
}

#source
ntfPatientFactorMatrix <- function(){
  weightsComp <- data.frame(namesC = paste('V',seq(1:kMax),sep=''),
                            weight=weightsC[1,] )
  weightsComp <- weightsComp[order(weightsComp$weight,decreasing = F),]
  patiF <- patientsF
  for (i in seq_along(weightsC)) {
    patiF[,i] <- patiF[,i]*weightsC[i]
  }
  #normalize rows
  patiF <- t(apply(patiF,1,function(x){x/sum(x)}))
  patiF <- patiF[,as.character(weightsComp$namesC) ]
}

NTF_AUC <- function()
{
  numSplits = 10
  AUC_CD_all <- rep(0,kMax)
  Cstat_allK <- rep(list(),kMax)
  smp_size <- floor(0.2 * numberOfPatiens) #testing sample size

  patiF <- ntfPatientFactorMatrix()
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[,-1]


  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in seq(kMax)){
    AUC_CD_K <- rep(0,numSplits)
    Cstat_K <- rep(0,numSplits)

    set.seed(k*1234)  # set seed for reproducibility
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = numSplits))

    for (i in seq(1,numSplits)) {
      test_ind <- obsInTest$subsets[,i]
      train <- patiF[-test_ind,c(c(kMax-k:kMax),kMax+1,kMax+2)]
      test <- patiF[test_ind,c(c(kMax-k:kMax),kMax+1,kMax+2)]

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
      Cstat <- BeggC(Surv.rsp, Surv.rsp.new, lp, lpnew)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, utimes)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K[i] <- AUC_CD$iauc
      Cstat_K[i] <- Cstat
    }
    #save.image(paste("temp/5-2v5-AUC_NTF_k",k,".RData",sep=''))
    AUC_CD_all[k] = mean(AUC_CD_K)
    Cstat_allK[k] <- list(Cstat_K)
  }

  NTF_AUC <- AUC_CD_all
  NTF_Cstat <- Cstat_allK
  save(NTF_AUC, file = "output4paper/NTF_AUC.RData")
  save(NTF_Cstat, file = "output4paper/NTF_Cstat.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-2v5.RData")
