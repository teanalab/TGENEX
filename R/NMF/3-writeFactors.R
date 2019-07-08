#read rubik output

rRange = c(3:78)


for(rankNMF in rRange)
{
  
  load(paste0("temp/NMF_",rankNMF,"_PxM.RData"))
  patientsFactors <- patientNMF_PxM
  mydata <- na.omit(patientsFactors)
  mydata <- scale(mydata) # standarize variables
  
  save(mydata,file=paste0("temp/patientFactors_r",rankNMF,".RData"))
}
