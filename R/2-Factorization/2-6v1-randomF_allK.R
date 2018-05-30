#' version 1 (May/30/2018)

rm(list = ls())
load(file="data/survivalClinical-5-1_30.RData")

randomPatientFactorMatrix <- function(k){
  set.seed(100)
  randPatfm <- matrix(runif(numberOfPatiens*k),numberOfPatiens,k)
  #normalize rows
  randPatfm <- t(apply(randPatfm,1,function(x){x/sum(x)}))
  randPatfm <- as.data.frame(randPatfm)
  randPatfm
}


#line by line
randomAUCs <- function()
{
  for (k in c(31:40)){
    #random factor
    randPatfm <-randomPatientFactorMatrix(k)
    set.seed(123*k)
    indexR <- sample(numberOfPatiens)
    randPatfm <- cbind.data.frame(randPatfm,survivalClinical[indexR,])
    patVarName <- paste("randomF_R_",k, sep='')
    assign(patVarName, randPatfm)
    save(list=c(patVarName), file = paste("data/random_norm/",patVarName,".RData",sep='') )
  }
}


sess <- sessionInfo() #save session on variable
save.image("temp/2-6v1.RData")
