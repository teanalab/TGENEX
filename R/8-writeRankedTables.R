# NON ZERO ---------------

#' Projection of the factor matrices
projection <- function(patientsF, mutationF, clinicalF, k = 7) {
  A <- patientsF
  B <- mutationF
  C <- clinicalF
  
  # these are published on the paper table 2
  percenPatients <- apply(A,2,function(x) {length(which(x != 0.0))/dim(A)[1]})
  percenGenes <- apply(B,2,function(x) {length(which(x != 0.0))/dim(B)[1]})
  percenClini <- apply(C,2,function(x) {length(which(x != 0.0))/dim(C)[1]})
  
  table2paper <- rbind(percenPatients,percenGenes,percenClini)
  table2paper <- t(table2paper)
  table2paper <- round(table2paper, 5)
  table2paper <- cbind(component=paste("Component", c(1:k)),
                       table2paper)
  
  #write.csv2(table2paper, file = "output4paper/table2.csv",row.names = FALSE)
  
  
  #numbers
  nPatients <- apply(A,2,function(x) {length(which(x != 0.0))})
  nGenes <- apply(B,2,function(x) {length(which(x != 0.0))})
  nClini <- apply(C,2,function(x) {length(which(x != 0.0))})
  
  setClass("AComponent",
           representation(
             patientVectorP="numeric",
             generVectorP="numeric",
             cinicalVectorP="numeric",
             patientVectorNames="character",
             generVectorNames="character",
             cinicalVectorNames="character",
             genesDF = "data.frame"))
  
  ListComponents <- list()
  
  for( i in 1:k) {
    ListComponents[[i]] <- new("AComponent",
                               patientVectorP=A[which(A[,i] != 0.0),i],
                               generVectorP=B[which(B[,i] != 0.0),i],
                               cinicalVectorP=C[which(C[,i] != 0.0),i],
                               patientVectorNames=row.names(A)[which(A[,i] != 0.0)],
                               generVectorNames=row.names(B)[which(B[,i] != 0.0)],
                               cinicalVectorNames=row.names(C)[which(C[,i] != 0.0)],
                               genesDF = data.frame(geneName = row.names(B)[which(B[,i] != 0.0)],
                                                    geneScore = B[which(B[,i] != 0.0),i], stringsAsFactors = F)
    )
  }
  
  #validate the lengths match the last component
  length(unlist(ListComponents[[i]]@patientVectorP))
  length(unlist(ListComponents[[i]]@generVectorP))
  length(unlist(ListComponents[[i]]@cinicalVectorP))
  length(unlist(ListComponents[[i]]@patientVectorNames))
  length(unlist(ListComponents[[i]]@generVectorNames))
  length(unlist(ListComponents[[i]]@cinicalVectorNames))
  ListComponents
}



writeTopGenes <- function(ListComponents, TopN = 100, r){
  for(i in seq_along(ListComponents)) {
    Gene_i <- data.frame(ListComponents[[i]]@genesDF, stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = T),]
    write.table(Gene_i[1:TopN,], paste('temp/Top',TopN,'NonZero_GenesxComponent_',i,'_r_',r,'.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}


writeAllGenes <- function(ListComponents){
  for(i in seq_along(ListComponents)) {
    Gene_i <- data.frame(ListComponents[[i]]@genesDF, stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = T),]
    write.table(Gene_i, paste('temp/GenesxComponent_',i,'_.csv',sep=''), sep=',', col.names = T, row.names = F )
  }
}

writeTopClinical <- function(ListComponents, TopN = 10, r){
  for(i in seq_along(ListComponents)) {
    Clini_i <- data.frame(clinival_var = ListComponents[[i]]@cinicalVectorNames,
                          score_clini = ListComponents[[i]]@cinicalVectorP,
                          stringsAsFactors = FALSE)
    Clini_i <- Clini_i[order(Clini_i$score_clini,decreasing = T),]
    write.table(Clini_i[1:TopN,], paste('temp/Top',TopN,'NonZero_ClinicalxComponent_',i,'_r_',r,'.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}



writeNonZeroClinical_scores <- function(){
  rm(list=ls()[-which(ls() %in% c("ListComponents", "k"))])
  
  for(i in 1:k) {
    Clini_i <- unlist(ListComponents[[i]]@cinicalVectorP)
    write.table( data.frame(cli=t(Clini_i)), 'temp/k7/NonZeroClinical_scores.csv', append= T, sep=',', col.names = F )
  }
}

writeNonZeroClinical <- function(){
  rm(list=ls()[-which(ls() %in% c("ListComponents", "k"))])
  
  for(i in 1:k) {
    Clini_i <- unlist(ListComponents[[i]]@cinicalVectorNames)
    write.table( data.frame(cli=t(Clini_i)), 'temp/k7/NonZeroClinical_names.csv', append= T, sep=',', col.names = F )
  }
}

# DISTANCE TO CENTROIDS ---------------


writeTopGenesDC <- function(distancesGenes, TopN = 100, k){
  for(i in seq(1,k)) {
    Gene_i <- data.frame(geneName=row.names(distancesGenes),
                         geneScore=distancesGenes[,i], stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = F),]
    write.table(Gene_i[1:TopN,], paste('temp/Top',TopN,'GenesxComponent_',i,'_k_',k,'.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}

writeTopClinicalDC <- function(distancesClinical, TopN = 10, k){
  for(i in seq(1,k) ) {
    Clini_i <- data.frame(clinival_var = row.names(distancesClinical),
                          score_clini = distancesClinical[,i], stringsAsFactors = FALSE)
    Clini_i <- Clini_i[order(Clini_i$score_clini, decreasing = F),]
    write.table(Clini_i[1:TopN,], paste('temp/Top',TopN,'ClinicalxComponent_',i,'_k_',k,'.csv',sep=''),
                sep=',', col.names = T, row.names = F )
  }
}

