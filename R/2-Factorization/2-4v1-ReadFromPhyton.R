#version 1 (May/25/2018)
#read results from python


clinicalF<-read.table(file="pythonOutput/clinicalMay.csv", sep = ",", quote = '"',
                     header = FALSE, stringsAsFactors = FALSE)
load(file= "data/boolClinical.RData")
row.names(clinicalF) <- names(boolClinical)

mutationF<-read.table(file="pythonOutput/genesMay.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
load(file= "data/boolMutation.RData")
row.names(mutationF) <- names(boolMutation)

patientsF<-read.table(file="pythonOutput/patientsMay.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
row.names(patientsF) <- patients

#lambda is all equal to .1. DO NOT USE THIS PACKAGE

