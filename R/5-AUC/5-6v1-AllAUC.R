# version 1 (May/27/2018)

rm(list = ls())
#load(file="data/survivalClinical-5-1.RData")
load(file="data/survivalClinical-5-1_30.RData")
load(file="data/factors.RData")
loadls("plyr survival Rcpp missForest survAUC perry",F)


#source
randomPatientFactorMatrix <- function(){
  randPatfm <- matrix(0,numberOfPatiens,kMax)
  for (i in seq(numberOfPatiens)){
    newRow <- runif(kMax)
    newRow <- newRow/sum(newRow)
    randPatfm[i,] <- newRow
  }
  randPatfm <- as.data.frame(randPatfm)
}

ntfPatientFactorMatrix <- function(){
  weightsComp <- data.frame(namesC = paste('V',seq(1:kMax),sep=''),
                            weight=weightsC[1,] )
  weightsComp[order(weightsComp$weight,decreasing = T),]
  patiF <- patientsF
  for (i in seq_along(weightsC)) {
    patiF[,i] <- patiF[,i]*weightsC[i]
  }
  #normalize rows
  patiF <- t(apply(patiF,1,function(x){x/sum(x)}))
  patiF
}

#line by line
getAllAUCs <- function()
{
  AUC_CD_allK_Random <- list()
  AUC_CD_allK_NTF <- list()
  nRepeats = 10
  percentesting = 0.2

  set.seed(100)
  #random factor
  randPatfm <-randomPatientFactorMatrix()
  randPatfm <- cbind.data.frame(randPatfm,survivalClinical)

  #ntf
  patiF <- ntfPatientFactorMatrix()
  patiF <- merge(as.data.frame(patiF), as.data.frame(survivalClinical), by='row.names', all=TRUE)
  patiF <- patiF[,-1]

  #Previously we used 10 random splits
  ## One 50% of the sample size
  #Now, we do bootstrap on 80%
  #testRandoming sample
  smp_size <- floor(percentesting * numberOfPatiens)

  ## data folds for K-fold cross-validation
  # perrySplits(20, foldControl(K = 5))
  # perrySplits(20, foldControl(K = 5, R = 10))
  ## random data splits
  set.seed(1234)  # set seed for reproducibility
  obsIntest <- perrySplits(numberOfPatiens, splitControl(m = smp_size, R = nRepeats))
  ## bootstrap samples
  # obsIntest <- perrySplits(numberOfPatiens, bootControl())


  for (k in seq(2,kMax)){
    AUC_CD_K_Random <- rep(0,nRepeats)
    AUC_CD_K_NTF <- rep(0,nRepeats)

    for (i in seq(1,nRepeats)) {
      test_ind <- obsIntest$subsets[,i]
      trainRandom <- randPatfm[-test_ind,c(1:k,kMax+1,kMax+2)]
      testRandom <- randPatfm[test_ind,c(1:k,kMax+1,kMax+2)]

      trainNTF <- patiF[-test_ind,c(kMax-k:kMax,kMax+1,kMax+2)]
      testNTF <- patiF[test_ind,c(kMax-k:kMax,kMax+1,kMax+2)]

      #random
      ## cox proportional hazard models
      coxFitRandom <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      model=TRUE, data=trainRandom)

      ## AUC
      lp <- predict(coxFitRandom)
      lpnew <- predict(coxFitRandom, newdata=testRandom)
      Surv.rsp <- Surv(time = trainRandom$Overall.Survival..Months.,
                       event = trainRandom$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = testRandom$Overall.Survival..Months.,
                           event = testRandom$Overall.Survival.Status)
      #times	The vector of time points at which AUC is evaluated.
      maxMonth <- max(survivalClinical$Overall.Survival..Months.)
      times <- seq(from = maxMonth/100, to = maxMonth, by = maxMonth/100)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K_Random[i] <- AUC_CD$iauc



      #NTF
      ## cox proportional hazard models
      coxFitNTF <- coxph(Surv(time = Overall.Survival..Months.,
                           event = Overall.Survival.Status)~.,
                      model=TRUE, data=trainNTF)
      #AUC
      lp <- predict(coxFitNTF)
      lpnew <- predict(coxFitNTF, newdata=testNTF)
      Surv.rsp <- Surv(time = trainNTF$Overall.Survival..Months.,
                       event = trainNTF$Overall.Survival.Status)
      Surv.rsp.new <- Surv(time = testNTF$Overall.Survival..Months.,
                           event = testNTF$Overall.Survival.Status)
      maxMonth <- max(survivalClinical$Overall.Survival..Months.)
      #times	The vector of time points at which AUC is evaluated.
      times <- seq(from = maxMonth/100, to = maxMonth, by = maxMonth/100)
      AUC_CD <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, times)
      #plot(AUC_CD, main = paste("CD", AUC_CD$iauc))
      AUC_CD_K_NTF[i] <- AUC_CD$iauc
    }
    AUC_CD_allK_Random[[k]] = AUC_CD_K_Random
    AUC_CD_allK_NTF[[k]] = AUC_CD_K_NTF
  }

  save(AUC_CD_allK_Random, AUC_CD_allK_NTF, file = "output4paper/All-AUC.RData")
}

sess <- sessionInfo() #save session on variable
save.image("temp/5-6v1.RData")
