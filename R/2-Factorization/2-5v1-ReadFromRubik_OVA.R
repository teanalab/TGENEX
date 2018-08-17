#version 1 (May/26/2018)
#read results from rubik
clinicalF<-read.table(file="rubikOutput/u3.csv", sep = ",", quote = '"',
                     header = FALSE, stringsAsFactors = FALSE)
load(file= "data/boolClinical.RData")
row.names(clinicalF) <- names(boolClinical)

mutationF<-read.table(file="rubikOutput/u1.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
load(file= "data/binaMutation.Rd")
row.names(mutationF) <- names(binaMutation)

patientsF<-read.table(file="rubikOutput/u2.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
load(file= "data/patients.RData")
row.names(patientsF) <- patients

weightsC <-read.table(file="rubikOutput/lambda.csv", sep = ",", quote = '"',
                      header = FALSE, stringsAsFactors = FALSE)
weightsC <- t(weightsC)

patientVsK <- as.data.frame(patientsF)
patientVsK[, "max"] <- apply(patientVsK, 1, max)
patientVsK[, "whichmax"] <- apply(patientVsK, 1, which.max)
patientAfiliationRubik_PxC <- data.frame(patient=patients,
                                       component=patientVsK$whichmax)


save(clinicalF,mutationF,patientsF,patients,weightsC, file="data/factors.RData")
save(patientAfiliationRubik_PxC, file="data/patientAfiliationRubik_PxC.Rd")
