# version 3 (May/26/2018)

rm(list = ls())
LoadMyData <- function()
{
  #Get A
  load("data/factors.RData",envir = .GlobalEnv)
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
  sapply(libs, function(u) {source(u)})
}

LoadMyData()

#' Projection of the factor matrices
projection <- function() {
  k=10
  A <- patientsF
  B <- mutationF
  C <- clinicalF

  # these are published on the paper table 2
  percenPatients <- apply(A,2,function(x) {length(which(x>0.0))/dim(A)[1]})
  percenGenes <- apply(B,2,function(x) {length(which(x>0.0))/dim(B)[1]})
  percenClini <- apply(C,2,function(x) {length(which(x>0.0))/dim(C)[1]})

  table2paper <- rbind(percenPatients,percenGenes,percenClini)
  table2paper <- t(table2paper)
  table2paper <- round(table2paper, 5)
  table2paper <- cbind(component=paste("Component", c(1:10)), table2paper)

  write.csv2(table2paper, file = "output4paper/table2.csv",row.names = FALSE)


  #numbers
  nPatients <- apply(A,2,function(x) {length(which(x>0.0))})
  nGenes <- apply(B,2,function(x) {length(which(x>0.0))})
  nClini <- apply(C,2,function(x) {length(which(x>0.0))})

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
                               patientVectorP=A[which(A[,i]>0.0),i],
                               generVectorP=B[which(B[,i]>0.0),i],
                               cinicalVectorP=C[which(C[,i]>0.0),i],
                               patientVectorNames=row.names(A)[which(A[,i]>0.0)],
                               generVectorNames=row.names(B)[which(B[,i]>0.0)],
                               cinicalVectorNames=row.names(C)[which(C[,i]>0.0)],
                               genesDF = data.frame(geneName = row.names(B)[which(B[,i]>0.0)],
                                                    geneScore = B[which(B[,i]>0.0),i], stringsAsFactors = F)
                               )
  }

  #validate the lengths match the last component
  length(unlist(ListComponents[[i]]@patientVectorP))
  length(unlist(ListComponents[[i]]@generVectorP))
  length(unlist(ListComponents[[i]]@cinicalVectorP))
  length(unlist(ListComponents[[i]]@patientVectorNames))
  length(unlist(ListComponents[[i]]@generVectorNames))
  length(unlist(ListComponents[[i]]@cinicalVectorNames))

  save(ListComponents,file="data/ListObjectComponents-3-1.RData")
  # load("data/ListObjectComponents-3-1.RData")
}


writeNonZeroClinical <- function(){
  rm(list=ls()[-which(ls() %in% c("ListComponents", "k"))])

  for(i in 1:k) {
    Clini_i <- unlist(ListComponents[[i]]@cinicalVectorNames)
    write.table( data.frame(cli=t(Clini_i)), 'output4paper/NonZeroClinical.csv', append= T, sep=',', col.names = F )
  }
}

writeTopGenes <- function(){
  rm(list=ls()[-which(ls() %in% c("ListComponents", "k"))])

  TopN = 100

  for(i in 1:k) {
    Gene_i <- data.frame(ListComponents[[i]]@genesDF, stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore),]
    ListComponents[[i]]@genesDF <- data.frame(Gene_i[1:TopN,], stringsAsFactors = FALSE)
    write.table(t(as.character(Gene_i$geneName[1:TopN])),
                 'output4paper/TopGenesXComponent.csv', append= T, sep=',', col.names = F )
  }

  save(ListComponents,file="data/ListObjectComponentsTop-3-1.RData")
}


writePatientsMembership <- function(){
  load('data/factors.RData')
  rm(list=ls()[-which(ls() %in% c("patientsF"))])

  #principal verotype from NTF assign
  A <- patientsF
  patientVsComponents <- as.data.frame(A)
  patientVsComponents[, "max"] <- apply(patientVsComponents, 1, max)
  patientVsComponents[, "whichmax"] <- apply(patientVsComponents, 1, which.max)
  patientAfiliation <- data.frame(patient=row.names(patientVsComponents),
                                  component=patientVsComponents$whichmax)

  save(patientAfiliation, file="data/patientAfiliation.RData")
}

#end
sess <- sessionInfo() #save session on variable
save.image(file="temp/3-1v2.RData")
# load(file="temp/3-1v2.RData")
