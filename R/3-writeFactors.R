#read rubik output

rRange = c(3:90)

load("data/patients.RData")

for(r in rRange)
{
  patientsFactors <- read.table(file=paste0("MATLAB-RUBIK/rubikOutput/R_",
                                            r,"/u2.csv"),
                                sep = ",", quote = "'",
                                header = FALSE, stringsAsFactors = FALSE,
                                skipNul = FALSE, blank.lines.skip = FALSE)
  mydata <- na.omit(patientsFactors)
  row.names(mydata) <- patients
  
  mydata <- scale(mydata) # standarize variables
  
  save(mydata,file=paste0("temp/patientFactors_r",r,".RData"))
}
