#version 2 (Feb-6-2018)


#analize the data
rm(list = ls())
load("data/ListObjectComponents-3-1.RData")

writeTopGenes <- function(){
  TopN = 500

  for(i in seq_along(ListComponents)) {
    Gene_i <- data.frame(ListComponents[[i]]@genesDF, stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = T),]
    write.table(Gene_i[1:TopN,], paste('temp/Top',TopN,'GenesxComponent_',i,'_.csv',sep=''), sep=',', col.names = T, row.names = F )
  }
}

writeAllGenes <- function(){
  for(i in seq_along(ListComponents)) {
    Gene_i <- data.frame(ListComponents[[i]]@genesDF, stringsAsFactors = FALSE)
    Gene_i <- Gene_i[order(Gene_i$geneScore,decreasing = T),]
    write.table(Gene_i, paste('temp/GenesxComponent_',i,'_.csv',sep=''), sep=',', col.names = T, row.names = F )
  }
}



