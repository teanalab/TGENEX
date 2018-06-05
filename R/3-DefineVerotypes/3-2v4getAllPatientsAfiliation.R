# version 4 (May/31/2018)
# all verotypes

rm(list = ls())

writePatientsAfiliationNTF <- function(k){
  for (k in c(2:30) )
  {
    load(paste("data/factors/factors_",k,".RData",sep='') )
    #principal verotype from NTF assign
    patiF <- patientsF
    #for (i in seq_along(weightsC)) {
    #  patiF[,i] <- patiF[,i]*weightsC[i]
    #}
    patientVsComponents <- as.data.frame(patiF)
    patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
    patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
    patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                    component=patientVsComponents$whichmax)
    save(patientAfiliation, file=paste("data/patientAfiliation/patientAfiliation_",
                                       k,".RData",sep='') )
  }
}

writePatientsAfiliationRandoms <- function(k){
  for (k in c(2:40) )
  {
    randPatfm <-randomPatientFactorMatrix(k)
    set.seed(123*k)
    indexR <- sample(numberOfPatiens)
    randPatfm <- cbind.data.frame(randPatfm,survivalClinical[indexR,])
    patVarName <- paste("randomF_R_",k, sep='')
    assign(patVarName, randPatfm)

    patientsF <- get(paste("randomF_R_",k,sep=''))
    patientsF <- patientsF[,c(1:k)]
    patientVsComponents <- as.data.frame(patientsF)
    patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
    patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
    patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                    component=patientVsComponents$whichmax)
    save(patientAfiliation, file=paste("data/patientA_random/patientA_random_",
                                       k,".RData",sep='') )
  }
}

writePatientsAfiliationNMF_PxM <- function(){
  for (k in c(2:30) )
  {
    load(file= paste("data/NMF/NMF_PxM_R_",k,".RData", sep='') )

    #principal verotype from NTF assign
    patiF <- get( paste("patientNMF_PxM_R_",k,sep='') )
    patientVsComponents <- as.data.frame(patiF)
    patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
    patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
    patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                    component=patientVsComponents$whichmax)
    save(patientAfiliation, file=paste("data/patientA_NMF_M/patientA_NMF_M_",
                                       k,".RData",sep='') )
  }
}

writePatientsAfiliationNMF_PxC <- function(){
  for (k in c(2:30) )
  {
    load(file= paste("data/NMF/NMF_PxC_R_",k,".RData", sep='') )
    #principal verotype from NTF assign
    patiF <- get( paste("patientNMF_PxC_R_",k,sep='') )
    patientVsComponents <- as.data.frame(patiF)
    patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
    patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
    patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                    component=patientVsComponents$whichmax)
    save(patientAfiliation, file=paste("data/patientA_NMF_C/patientA_NMF_C_",
                                       k,".RData",sep='') )
  }
}


#end
sess <- sessionInfo() #save session on variable
save.image(file="temp/3-1v4.RData")
