# version 3 (Jun/3/2018)

rm(list = ls())
LoadMyData <- function()
{
  k <- 7
  folder <- paste("temp/NMF_clini/k",k, sep ='')
  load( paste('data/NMF/NMF_PxC_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
}


writeTopClini <- function(){
  TopN = 10

  clini_table <- get(paste('mutationNMF_PxC_R_',k,sep=''))
  for(i in seq(k)) {
    clini_i <- data.frame(gene = colnames(clini_table), geneScore = clini_table[i,],
                         stringsAsFactors = FALSE)
    clini_i <- clini_i[order(clini_i$geneScore,decreasing = T),]
    write.table(clini_i[1:TopN,], paste(folder,'/Top',TopN,'ClinicalxComponent_',i,'_.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}


writeTopPatient <- function(){
  TopN = 50

  patient_table <- get(paste('patientNMF_PxM_R_',k,sep=''))
  for(i in seq(k)) {
    patient_i <- data.frame(patient = row.names(patient_table),
                         score_patient = patient_table[,i],
                         stringsAsFactors = FALSE)
    patient_i <- patient_i[order(patient_i$score_patient, decreasing = T),]
    write.table(patient_i[1:TopN,], paste('temp/k',k,'/Top',TopN,'patientcalxComponent_',i,'_.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}


