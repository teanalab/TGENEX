# version 2 (May/26/2018)
# AUC with random clustering
# assing random values to the components
# the sum of components per each patient should sum to 1

rm(list = ls())
load(file="data/survivalClinical-5-1.RData")
loadls("plyr survival Rcpp missForest survAUC perry",FALSE)

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
  kMax <- 10
  row.names(survivalClinical) <- survivalClinical$Patient.ID

  save.image(file="data/survivalClinical-5-1.RData")
}

#source
randomPatientFactorMatrix <- function(){
  randPatfm <- matrix(0,numberOfPatiens,kMax)
  for (i in seq(numberOfPatiens)){
    newRow <- runif(10)
    newRow <- newRow/sum(newRow)
    randPatfm[i,] <- newRow
  }
  randPatfm <- as.data.frame(randPatfm)
}

#line by line
randomAUCs <- function()
{
  AUC_CD_allK <- rep(0,kMax)
  nRepeats = 10


  set.seed(100)
  #random factor
  randPatfm <-randomPatientFactorMatrix()
  randPatfm <- cbind.data.frame(randPatfm,survivalClinical)

  for (k in seq(2,kMax)){
    AUC_CD_K <- rep(0,nRepeats)

    #Previously we used 10 random splits
    ## One 50% of the sample size
    #Now, we do bootstrap on 80%
    set.seed(k*1234)  # set seed for reproducibility
    #testing sample
    smp_size <- floor(0.2 * numberOfPatiens)

    ## data folds for K-fold cross-validation
    # perrySplits(20, foldControl(K = 5))
    # perrySplits(20, foldControl(K = 5, R = 10))



    ## random data splits
    obsInTest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = nRepeats))


    ## bootstrap samples
    # obsInTest <- perrySplits(numberOfPatiens, bootControl())

    for (i in seq(1,nRepeats)) {
      set.seed(i)
      test_ind <- obsInTest[[4]][,i]
      train <- randPatfm[-test_ind, ]
      test <- randPatfm[test_ind, ]

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
    AUC_CD_allK[k] = mean(AUC_CD_K)
  }

  randomAUC <- AUC_CD_allK
  save(randomAUC, file = "output4paper/randomAUC.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-1v2.RData")
