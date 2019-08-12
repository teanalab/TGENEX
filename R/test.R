load("~/DDfileSystem/10-Research/1-Current/CLINGEN/9source/CLIGEN/git_teana/CLIGEN/allmRNA.RData")

initialScriptBRCAonly <-function(){
  #all with BRCA olny
  row.names( BRCA.mRNA$bcr_patient_barcode )
  BRCA.mRNA = BRCA.mRNA[, -1]
  BRCA.mRNA[0:4,0:4]


  library(gplots)

  heatmap.2(t(BRCA.mRNA)[1:10,1:10]) # Shortcut to final result


  y = BRCA.mRNA


  #remove NA
  anyNA(y)
  sum(is.na(y))

  #1695
  #delete columns because patients are more valuable
  delCols <- apply(y, 1, anyNA )
  sum(delCols)

  remNA = which(delCols)
  y <- y[,-remNA]

  anyNA(y)

  ## Row- and column-wise clustering
  hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete")
  hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")

  library(dendextend)
  library(circlize)

  # create a dendrogram
  dend <- as.dendrogram(hr)

  # modify the dendrogram to have some colors in the branches and labels
  dend <- dend %>%
    color_branches(k=4) %>%
    color_labels

  # plot the radial plot
  par(mar = rep(0,4))
  # circlize_dendrogram(dend, dend_track_height = 0.8)
  circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
  #circlize_dendrogram(dend, labels=F, labels_track_height = NA, dend_track_height = .4)


  pdf( file = "BC_dendo.pdf",  onefile = TRUE, width = 9, height = 9)
  #pdf( file = paste0(figureKs,"/kaplan-NBS",(i+1),".pdf"),  onefile = TRUE, width = 9, height = 7)
  #tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  circlize_dendrogram(dend, labels=F, labels_track_height = NA, dend_track_height = .4)
  dev.off()
}



#' Title
#'
#' @param TCGA_DB
#' @param fileName
#'
#' @return
#' @export
#'
#' @examples
#' library(gplots)
#' library(dendextend)
#' library(circlize)
#' TCGA_DB = COAD.mRNA
#' fileName = "COAD_dendo.pdf"
#' saveDendograms(TCGA_DB,fileName)
saveDendograms <- function( TCGA_DB, fileName ){
  row.names(TCGA_DB) <- TCGA_DB$bcr_patient_barcode
  TCGA_DB = TCGA_DB[, -1]

  #remove NA
  if(anyNA(TCGA_DB)){
    nCols = dim(TCGA_DB)[2]
    cat(sum(is.na(TCGA_DB)), " NA values found for ", fileName, "\n")
    #delete genes because patients are more valuable
    delCols <- apply(TCGA_DB, 2, anyNA )
    cat(sum(delCols), " genes to be deleted\n")
    remNA = which(delCols)
    TCGA_DB <- TCGA_DB[,-remNA]
    cat( (dim(TCGA_DB)[2] -  nCols), "genes were deleted", "\n")
  }
  if ( anyNA(TCGA_DB) ){
    stop("ERROR: Table stll have NA values for ", fileName)
  }

  cat( "For", fileName, "Number of genes", dim(TCGA_DB)[2], "Number of patients", dim(TCGA_DB)[1] , "\n")

  ## Row- and column-wise clustering
  hr <- hclust(as.dist(1-cor(t(TCGA_DB), method="pearson")), method="complete")
  #hc <- hclust(as.dist(1-cor(TCGA_DB, method="spearman")), method="complete")
  # create a dendrogram
  dend <- as.dendrogram(hr)
  # modify the dendrogram to have some colors in the branches and labels
  dend <- dend %>%
    color_branches(k=4) %>%
    color_labels
  # plot the radial plot
  par(mar = rep(0,4))
  # circlize_dendrogram(dend, dend_track_height = 0.8)
  #circlize_dendrogram(dend, labels_track_height = NA, dend_track_height = .4)
  circlize_dendrogram(dend, labels=F, labels_track_height = NA, dend_track_height = .4)

  pdf( file = fileName,  onefile = TRUE, width = 9, height = 9)
  #pdf( file = paste0(figureKs,"/kaplan-NBS",(i+1),".pdf"),  onefile = TRUE, width = 9, height = 7)
  #tiff(file = "temp/79-kaplan-NMF4_num.tiff",  width = 9, height = 7, units = 'in', res=300)
  circlize_dendrogram(dend, labels=F, labels_track_height = NA, dend_track_height = .4)
  dev.off()
}




library(gplots)
library(dendextend)
library(circlize)



TCGA_DB = BRCA.mRNA
fileName = "BRCA.pdf"
saveDendograms(TCGA_DB,fileName)


TCGA_DB = COAD.mRNA
fileName = "COAD.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = COADREAD.mRNA
fileName = "COADREAD.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = GBMLGG.mRNA
fileName = "GBMLGG.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = KIPAN.mRNA
fileName = "KIPAN.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = KIRC.mRNA
fileName = "KIRC.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = KIRP.mRNA
fileName = "KIRP.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = LGG.mRNA
fileName = "LGG.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = LUAD.mRNA
fileName = "LUAD.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = LUSC.mRNA
fileName = "LUSC.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = OV.mRNA
fileName = "OV.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = READ.mRNA
fileName = "READ.pdf"
saveDendograms(TCGA_DB,fileName)

TCGA_DB = UCEC.mRNA
fileName = "UCEC.pdf"
saveDendograms(TCGA_DB,fileName)



outputAug11 <- function(){

  > TCGA_DB = BRCA.mRNA

  > fileName = "BRCA.pdf"

  > saveDendograms(TCGA_DB,fileName)
  1695  NA values found for  BRCA.pdf
  534  genes to be deleted
  -534 genes were deleted
  For BRCA.pdf Number of genes 17280 Number of patients 590
  RStudioGD
  2

  > TCGA_DB = COAD.mRNA

  > fileName = "COAD.pdf"

  > saveDendograms(TCGA_DB,fileName)
  445  NA values found for  COAD.pdf
  274  genes to be deleted
  -274 genes were deleted
  For COAD.pdf Number of genes 17540 Number of patients 172
  RStudioGD
  2

  > TCGA_DB = COADREAD.mRNA

  > fileName = "COADREAD.pdf"

  > saveDendograms(TCGA_DB,fileName)
  546  NA values found for  COADREAD.pdf
  309  genes to be deleted
  -309 genes were deleted
  For COADREAD.pdf Number of genes 17505 Number of patients 244
  RStudioGD
  2

  > TCGA_DB = GBMLGG.mRNA

  > fileName = "GBMLGG.pdf"

  > saveDendograms(TCGA_DB,fileName)
  186  NA values found for  GBMLGG.pdf
  151  genes to be deleted
  -151 genes were deleted
  For GBMLGG.pdf Number of genes 17663 Number of patients 27
  RStudioGD
  2

  > TCGA_DB = KIPAN.mRNA

  > fileName = "KIPAN.pdf"

  > saveDendograms(TCGA_DB,fileName)
  783  NA values found for  KIPAN.pdf
  554  genes to be deleted
  -554 genes were deleted
  For KIPAN.pdf Number of genes 17260 Number of patients 88
  RStudioGD
  2

  > TCGA_DB = KIRC.mRNA

  > fileName = "KIRC.pdf"

  > saveDendograms(TCGA_DB,fileName)
  373  NA values found for  KIRC.pdf
  292  genes to be deleted
  -292 genes were deleted
  For KIRC.pdf Number of genes 17522 Number of patients 72
  RStudioGD
  2

  > TCGA_DB = KIRP.mRNA

  > fileName = "KIRP.pdf"

  > saveDendograms(TCGA_DB,fileName)
  410  NA values found for  KIRP.pdf
  344  genes to be deleted
  -344 genes were deleted
  For KIRP.pdf Number of genes 17470 Number of patients 16
  RStudioGD
  2

  > TCGA_DB = LGG.mRNA

  > fileName = "LGG.pdf"

  > saveDendograms(TCGA_DB,fileName)
  186  NA values found for  LGG.pdf
  151  genes to be deleted
  -151 genes were deleted
  For LGG.pdf Number of genes 17663 Number of patients 27
  RStudioGD
  2

  > TCGA_DB = LUAD.mRNA

  > fileName = "LUAD.pdf"

  > saveDendograms(TCGA_DB,fileName)
  517  NA values found for  LUAD.pdf
  418  genes to be deleted
  -418 genes were deleted
  For LUAD.pdf Number of genes 17396 Number of patients 32
  RStudioGD
  2

  > TCGA_DB = LUSC.mRNA

  > fileName = "LUSC.pdf"

  > saveDendograms(TCGA_DB,fileName)
  532  NA values found for  LUSC.pdf
  304  genes to be deleted
  -304 genes were deleted
  For LUSC.pdf Number of genes 17510 Number of patients 154
  RStudioGD
  2

  > TCGA_DB = OV.mRNA

  > fileName = "OV.pdf"

  > saveDendograms(TCGA_DB,fileName)
  755  NA values found for  OV.pdf
  375  genes to be deleted
  -375 genes were deleted
  For OV.pdf Number of genes 17439 Number of patients 561
  RStudioGD
  2

  > TCGA_DB = READ.mRNA

  > fileName = "READ.pdf"

  > saveDendograms(TCGA_DB,fileName)
  101  NA values found for  READ.pdf
  85  genes to be deleted
  -85 genes were deleted
  For READ.pdf Number of genes 17729 Number of patients 72
  RStudioGD
  2

  > TCGA_DB = UCEC.mRNA

  > fileName = "UCEC.pdf"

  > saveDendograms(TCGA_DB,fileName)
  556  NA values found for  UCEC.pdf
  377  genes to be deleted
  -377 genes were deleted
  For UCEC.pdf Number of genes 17437 Number of patients 54
  RStudioGD
  2
}
