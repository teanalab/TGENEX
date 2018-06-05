# version 3 (Jun/3/2018)

rm(list = ls())
LoadMyData <- function()
{
  k <- 12
  folder <- paste("temp/NMF_mut/k",k, sep ='')
  load( paste('data/NMF/NMF_PxM_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
}


writeTopGenes <- function(){
  TopN = 100

  mut_table <- get(paste('mutationNMF_PxM_R_',k,sep=''))
  for(i in seq(k)) {
    Gene_i <- data.frame(gene = colnames(mut_table), geneScore = mut_table[i,],
                         stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = T),]
    write.table(Gene_i[1:TopN,], paste(folder,'/Top',TopN,'GenesxComponent_',i,'_.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}


writeTopPatient <- function(){
  TopN = 10

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


