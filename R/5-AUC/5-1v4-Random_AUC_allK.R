#' version 4 (May/30/2018)
#' AUC with random clustering

rm(list = ls())
load(file="data/survivalClinical-5-1_30.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)

LoadMyData <- function(){
  #load my libraries
  loadlib<-c("Packages.R")
  loadlib<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", loadlib,sep='')
  sapply(loadlib, function(u) {source(u)})
  #required for the cox proportional hazard model
  loadls("plyr survival Rcpp missForest survAUC",FALSE)

  load("data/survClinical.RData")
  load("data/patients.RData")


  survivalClinical <- survClinical
  survivalClinical[,"Patient.ID"] <- survivalClinical$patient.bcr_patient_barcode
  survivalClinical[,"Overall.Survival.Status"] <- (survivalClinical$patient.vital_status=="dead")
  survivalClinical[,"Overall.Survival..Months."] <- survivalClinical$OS_MONTHS

  survivalClinical <-  survivalClinical[which(survivalClinical$Patient.ID %in% patients),]
  numberOfPatiens <- dim(survivalClinical)[[1]]
  kMax <- 30
  row.names(survivalClinical) <- survivalClinical$Patient.ID
  survivalClinical <- survivalClinical[,c("Overall.Survival.Status","Overall.Survival..Months.")]
  save.image(file="data/survivalClinical-5-1_30.RData")
}

#line by line
randomAUCs <- function()
{
  kMax = 40
  AUC_CD_allK <- rep(list(),kMax)
  Cstat_allK <- rep(list(),kMax)
  nSplits = 100
  percenTesting = 0.2

  set.seed(100)

  #times	The vector of time points at which AUC is evaluated.
  utimes <- unique( survivalClinical[survivalClinical$Overall.Survival.Status == 1,"Overall.Survival..Months."] )
  utimes <- utimes[ order(utimes) ]

  for (k in c(2:kMax) ){
    varName = paste("randomF_R_",k,sep='')
    load(paste("data/random_norm/",varName,".RData",sep=''))
    #random factor
    assign("randPatfm", get(varName) )

    AUC_CD_K <- rep(0,nSplits)
    Cstat_K <- rep(0,nSplits)

    set.seed(k*1234)  # set seed for reproducibility
    smp_size <- floor(0.2 * numberOfPatiens)  #testing sample

    ## data folds for K-fold cross-validation
    # perrySplits(20, foldControl(K = 5))
    # perrySplits(20, foldControl(K = 5, R = 10))
    ## random data splits
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = nSplits))
    ## bootstrap samples
    # obsInTest <- perrySplits(numberOfPatiens, bootControl())

    for (i in seq(1,nSplits)) {
      set.seed(i)
      test_ind <- obsInTest$subsets[,i]
      train <- randPatfm[-test_ind, ]
      test <- randPatfm[test_ind, ]

      #cox proportional hazard model
      coxFit <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      model=TRUE, data=train)
      ## AUC
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
    save.image(paste("temp/5-1v3-AUC_Random_k",k,".RData",sep=''))
    AUC_CD_allK[k] = list(AUC_CD_K)
    Cstat_allK[k] <- list(Cstat_K)
  }
  randomAUC <- AUC_CD_allK
  randomCstat <- Cstat_allK
  save(randomAUC, file = "output4paper/randomAUC.RData")
  save(randomCstat, file = "output4paper/randomCstat.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-1v4.RData")
