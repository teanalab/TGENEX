#' version 1 (May/30/2018)

rm(list = ls())
load(file="data/survivalClinical-5-1_30.RData")



#line by line
randomAUCs <- function()
{
  set.seed(1)
  for (k in c(1:30)){
    #random factor
    randPatfm <- data.frame(patients= patients, component=sample(k, numberOfPatiens, replace = T) )
    patVarName <- paste("randomF_R_",k, sep='')
    assign(patVarName, randPatfm)
    save(list=c(patVarName), file = paste("data/random/",patVarName,".RData",sep='') )
  }
}


sess <- sessionInfo() #save session on variable
save.image("temp/2-6v1.RData")
