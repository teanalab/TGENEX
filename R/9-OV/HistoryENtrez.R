rm(list = ls())
load(file="temp/boolMutation_OVA.RData")


#write.csv2(data.frame(names(boolMutation)),
#           file = "temp/genes.csv",row.names = FALSE)
#get genes from https://www.genenames.org/cgi-bin/symbol_checker
gene4genenames <-read.table(file="temp/1-geneFromGeneNames.txt", sep = "\t",
                            quote = "",
                            header = TRUE, stringsAsFactors = FALSE,
                            skipNul = TRUE, #skip = 1,
                            blank.lines.skip = TRUE,
                            fill=TRUE)
gene4genenames <- gene4genenames[!duplicated(gene4genenames$Input),]
genes2del <- which(gene4genenames$Match.type == "Unmatched")
gene4genenames <- gene4genenames[-genes2del,]


#from my gist
library("org.Hs.eg.db")
Sym2Entrez <- get("org.Hs.egSYMBOL2EG")
genes2map <- gene4genenames$Approved.symbol
genes2mapI <- gene4genenames$Input

genesEntrez <- rep(0, length(genes2map))
for (Gi in seq_along(genes2map) )
{
  newInstruction <- try( genesEntrez[Gi] <- as.list(Sym2Entrez[genes2map[Gi]])[[1]],
                         silent=TRUE )
  if( inherits(newInstruction, "try-error") )
  {
    try( genesEntrez[Gi] <- as.list(Sym2Entrez[genes2mapI[Gi]])[[1]] ,
         silent=FALSE )
  }
}

genesEntrezL <- unlist(genesEntrez, use.names = FALSE)
#make sure no value is equal to 0
#range(genesEntrezL)

boolMutationRedu <- boolMutation[,genes2mapI]
names(boolMutationRedu) <- genes2map #assign approve symbols
binaMutation <- as.data.frame(lapply(boolMutationRedu, as.numeric), stringsAsFactors = FALSE)
write.csv2(t(binaMutation), file = "temp/1-binaMredu_OVA.csv",row.names = FALSE)
write.csv2(genesEntrezL, file = "temp/10-genes_entrez_OVA.csv",row.names = FALSE)
write.csv2(genes2map, file = "temp/10-genes_redu_OVA.csv",row.names = FALSE)

write.csv2(data.frame(row.names(boolMutation)), file = "temp/1-patients.csv",row.names = FALSE)


