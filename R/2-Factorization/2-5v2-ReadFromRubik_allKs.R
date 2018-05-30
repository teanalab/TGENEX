#version 2 (May/30/2018)
#read results from rubik


rm(list = ls())
load(file= "data/boolClinical.RData")
load(file= "data/boolMutation.RData")


for ( k in c(31:40) )
{
  varName = paste("rubikOutput/R_",k,"/u3.csv",sep='')
  clinicalF<-read.table(file=varName, sep = ",", quote = '"',
                       header = FALSE, stringsAsFactors = FALSE)
  row.names(clinicalF) <- names(boolClinical)

  varName = paste("rubikOutput/R_",k,"/u1.csv",sep='')
  mutationF<-read.table(file=varName, sep = ",", quote = '"',
                        header = FALSE, stringsAsFactors = FALSE)
  row.names(mutationF) <- names(boolMutation)

  varName = paste("rubikOutput/R_",k,"/u2.csv",sep='')
  patientsF<-read.table(file=varName, sep = ",", quote = '"',
                        header = FALSE, stringsAsFactors = FALSE)
  row.names(patientsF) <- patients

  varName = paste("rubikOutput/R_",k,"/lambda.csv",sep='')
  weightsC <-read.table(file=varName, sep = ",", quote = '"',
                        header = FALSE, stringsAsFactors = FALSE)
  weightsC <- t(weightsC)

  save(clinicalF,mutationF,patientsF,patients,weightsC,
       file = paste("data/factors_",k,".RData", sep='') )
}
