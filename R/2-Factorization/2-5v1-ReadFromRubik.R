#version 1 (May/26/2018)
#read results from rubik


clinicalF<-read.table(file="rubikOutput/u3.csv", sep = ",", quote = '"',
                     header = FALSE, stringsAsFactors = FALSE)
load(file= "data/boolClinical.RData")
row.names(clinicalF) <- names(boolClinical)

mutationF<-read.table(file="rubikOutput/u1.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
load(file= "data/boolMutation.RData")
row.names(mutationF) <- names(boolMutation)

patientsF<-read.table(file="rubikOutput/u2.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
row.names(patientsF) <- patients


weightsC <-read.table(file="rubikOutput/lambda.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)

weightsC <- t(weightsC)


save(clinicalF,mutationF,patientsF,patients,weightsC, file="data/factors.RData")

