
#' Read a mutation file and construct the boolean matrix
#'
#' @return
#' @export
#'
#' @examples #TGCA data
#' boolMutation <- createBoolMutation(fileName = "/Users/diamac/10-Research/4Current/CLINGEN/12Datasets/ovarian/data_mutations_extended.txt")
#'     length(genes)
#'     length(patients)
#'     dim(boolMutation)
#'     #description boolean table
#'     genes <- names(boolMutation)
#'     patients <- row.names(boolMutation)
#'     save(boolMutation,patients, file="temp/boolMutation_OVA.RData")
createBoolMutation <- function(fileName, minNumMutation = NA, boolOrBin="bool"){
  mutation<-read.table(file = fileName, sep = "\t", quote = "",
                       header = TRUE, stringsAsFactors = FALSE, skip = 1, skipNul = TRUE, blank.lines.skip = TRUE)

  #uuid = NO A GOOD IDENTIFIER
  #Tumor_Sample_Barcode  = GOOD IDENTIFIER, BUT WE NEED TO CUT
  #TCGA-06-5411-01A-01D-1696-08
  mutationBarcodes <- substring(tolower(mutation$Tumor_Sample_Barcode), 1, 12)
  mutation$mutationBarcodes <- mutationBarcodes
  patients <- unique(mutationBarcodes)

  genes <- unique(as.character(mutation$Hugo_Symbol))


  if(boolOrBin == "bool")
  {
    boolMutation <- data.frame(matrix(FALSE,length(patients),length(genes)),
                               row.names = patients)
    trueVal = TRUE
  } else{
    boolMutation <- data.frame(matrix(0,length(patients),length(genes)),
                               row.names = patients)
    trueVal = 1
  }
  names(boolMutation) <- genes

  for(i in seq_len(dim(mutation)[1])){
    Apatient <- as.character(mutation$mutationBarcodes[i])
    if(Apatient %in% patients){
      #only non-silent mutations
      TMuta <- unique(mutation$Variant_Classification[i])
      if(TMuta != "Silent"){
        APgene <- mutation$Hugo_Symbol[i]
        boolMutation[Apatient,APgene] <- trueVal
      }
    }
  }

  ####
  #plot distributions
  #1 Number of mutations per patient
  mut4patient <- apply(boolMutation,1,sum)


 if( !is.na(minNumMutation) ){
   boolMutation <- boolMutation[-which(mut4patient<=10),]
   #after deletion
   mut4patient <- apply(boolMutation,1,sum)
 }


  # pdf(file = "temp/1-2v3-histomut4pati.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  # hist(mut4patient, main = "Histogram mutations per patient (mutation count)" ,
  #      ylab = "number of patients" ,
  #      xlab = "number of genes with mutations",
  #      breaks = 200)
  # dev.off()

  #Number of patients with mutation x
  pat4genes <- apply(boolMutation,2,sum)
  #Genes to remove
  if ( length(which(pat4genes<=0)) > 0 ){
    boolMutation <- boolMutation[,-which(pat4genes<=0)]
  }



  #    ########## REMOVING PATIENTS
  # #again after delete 0s
  # pat4genes <- apply(boolMutation,2,sum)
  # range(pat4genes)
  #   #NOTE: Genes to remove with 1
  #   length(genes)
  #   length(which(pat4genes==1))
  #   boolMutation <- boolMutation[,-which(pat4genes==1)]
  #   genes <- unique(names(boolMutation))
  #   length(genes)
  #
  #   #again after delete 0s
  #   pat4genes <- apply(boolMutation,2,sum)
  #   range(pat4genes)
  #
  #   #NOTE: Genes to remove with 2
  #   length(genes)
  #   length(which(pat4genes==2))
  #   boolMutation <- boolMutation[,-which(pat4genes==2)]
  #   genes <- unique(names(boolMutation))
  #   length(genes)
  #
  #   #again after delete 0s
  #   pat4genes <- apply(boolMutation,2,sum)
  #   range(pat4genes)
  #
  #   # #NOTE: Genes to remove with 3
  #   # length(genes)
  #   # length(which(pat4genes==3))
  #   # boolMutation <- boolMutation[,-which(pat4genes==3)]
  #   # genes <- unique(names(boolMutation))
  #   # length(genes)
  #   #
  #   # #again after delete 0s
  #   # pat4genes <- apply(boolMutation,2,sum)
  #   # range(pat4genes)
  #   #
  #   # #NOTE: Genes to remove with 4
  #   # length(genes)
  #   # length(which(pat4genes==4))
  #   # boolMutation <- boolMutation[,-which(pat4genes==4)]
  #   # genes <- unique(names(boolMutation))
  #   # length(genes)
  #   #
  #   # #again after delete 4s
  #   # pat4genes <- apply(boolMutation,2,sum)
  #   # range(pat4genes)
  #
  #
  #   pdf(file = "temp/1-2v3-Fig4.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  #   hist(pat4genes, main = "Histogram patient per gene (patient count)",
  #        ylab = "Number of genes",
  #        xlab = "Number of mutations", breaks=70)
  #   dev.off()
  #
  #
  return (boolMutation)
}
