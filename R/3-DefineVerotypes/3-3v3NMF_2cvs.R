# version 3 (Jun/3/2018)

rm(list = ls())

files2cvs_mutation7 <- function(){
  k <- 7
  folder <- paste("temp/NMF_mut/k",k, sep ='')
  load( paste('data/NMF/NMF_PxM_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
  write.table(mutationNMF_PxM_R_7, paste(folder,'/mutationNMF_PxM_R_7.csv',sep=''),
              sep=',', col.names = T, row.names = T )

  write.table(patientNMF_PxM_R_7, paste(folder,'/patientNMF_PxM_R_7.csv',sep=''),
              sep=',', col.names = T, row.names = F )
}

files2cvs_mutation12 <- function(){
  k <- 12
  folder <- paste("temp/NMF_mut/k",k, sep ='')
  load( paste('data/NMF/NMF_PxM_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
  write.table(mutationNMF_PxM_R_12, paste(folder,'/mutationNMF_PxM_R_12.csv',sep=''),
              sep=',', col.names = T, row.names = T )

  write.table(patientNMF_PxM_R_12, paste(folder,'/patientNMF_PxM_R_12.csv',sep=''),
              sep=',', col.names = T, row.names = F )
}


files2cvs_clinical7 <- function(){
  k <- 7
  folder <- paste("temp/NMF_clini/k",k, sep ='')
  load( paste('data/NMF/NMF_PxC_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
  write.table(mutationNMF_PxC_R_7, paste(folder,'/mutationNMF_PxC_R_7.csv',sep=''),
              sep=',', col.names = T, row.names = T )

  write.table(patientNMF_PxC_R_7, paste(folder,'/patientNMF_PxC_R_7.csv',sep=''),
              sep=',', col.names = T, row.names = F )
}


files2cvs_clinical12 <- function(){
  k <- 12
  folder <- paste("temp/NMF_clini/k",k, sep ='')
  load( paste('data/NMF/NMF_PxC_R_',k,'.RData', sep ='') ,envir = .GlobalEnv)

  if( !file.exists(folder) ) dir.create(folder)
  write.table(mutationNMF_PxC_R_12, paste(folder,'/mutationNMF_PxC_R_12.csv',sep=''),
              sep=',', col.names = T, row.names = T )

  write.table(patientNMF_PxC_R_12, paste(folder,'/patientNMF_PxC_R_12.csv',sep=''),
              sep=',', col.names = T, row.names = F )
}

