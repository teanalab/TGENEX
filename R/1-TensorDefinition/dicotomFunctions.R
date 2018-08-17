
#plots for clinical variables
checkCol <- function(column, main=""){
  barplot(table(column), main=main,
          legend.text = c("> Num Vars: ", length(unique(column)),
                          "> VALUES: " ,unique(column) ),
          args.legend = list(x ='topright', bty='n'))
}

#plots for clinical variables
analyzeClinicVariables <- function(){
  pdf(file = "temp/10-Fig5_ClinicalVariablesRaw.pdf" ,  onefile = TRUE, pagecentre = FALSE, compress = FALSE)
  NumVals <- c()
  for(i in seq_along(names(clinical))){
    checkCol(clinical[,i],i)
    NumVals <- c(NumVals, length(unique(clinical[,i])) )
  }
  dev.off()

  # output the number of values per each column
  # to see if the variable is dichotomous or not
  NumVals
}


#1.1.6.	Dichotomization of clinical variables with more than 2 categories
#and blank values
deleteBlanks <- function(columnId)
{
  matCols <- as.data.frame(clinical[, columnId])
  names(matCols) <- names(clinical)[columnId]
  matCols[is.na(matCols)]<-''
  matCols <- lapply(matCols, factor)

  ##Recode categories to columns
  mati <- model.matrix(~ . + 0, data=matCols,
                       contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
  mati <- as.data.frame(mati)

  if(dim(mati)[1] != dim(clinical)[1])
    stop("ERROR. rows lost")

  #delete blank cols = col with shorter name
  nchars <- sapply(names(mati),nchar)
  ToDel <- which(nchars == min(nchars))
  if(length(ToDel) == 0){
    stop("No col to delete")
  } else if(length(ToDel) > 1){
    stop("more than one col to delete")
  } else {
    mati <- mati[,-ToDel]
  }
  mati
}


#Dichomotization of columns with more than 2 categories
#and [Not available] values
dichomoNotAvail <- function(columnId)
{
  matCols <- as.data.frame(clinical[, columnId])
  names(matCols) <- names(clinical)[columnId]
  matCols <- lapply(matCols, factor)

  ##Recode categories to columns
  mati <- model.matrix(~ . + 0, data=matCols,
                       contrasts.arg = lapply(matCols, contrasts, contrasts=FALSE))
  mati <- as.data.frame(mati)

  #delete [Not Available] cols
  ToDel <- grep("\\[", names(mati))
  if(length(ToDel) !=0){
    mati <- mati[,-ToDel]
  }
  mati
}

#Dichomotization of columns with more than 2 categories
dichoSeveralCols <-function(cols)
{
  dichoCols <- as.data.frame(clinical[patients, cols])
  dichoCols <- lapply(dichoCols, factor)

  ##Recode categories to columns
  matDich <- model.matrix(~ . + 0, data=dichoCols,
                          contrasts.arg = lapply(dichoCols, contrasts, contrasts=FALSE))
  matDich <- as.data.frame(matDich)
  matDich <- as.data.frame(lapply(matDich, as.logical))
  row.names(matDich) <- patients
  matDich
}


#Dichomotization of discrete columns with several categories
classesSturges <- function(colId)
{
  require(classInt)
  aInt <- classIntervals(as.numeric(as.character(clinical[patients,colId])), style = "jenks")
  aCategory <- findCols(aInt)

  ##Recode categories to columns
  boolClinical <- data.frame(aCategory)
  boolClinical <- lapply(boolClinical, factor)
  matA <- model.matrix(~ . + 0, data=boolClinical,
                       contrasts.arg = lapply(boolClinical, contrasts, contrasts=FALSE))
  matA <- as.data.frame(matA)

  interNames <-  print(aInt)
  interNames <- names(interNames)

  names(matA) <- paste(names(clinical)[colId],interNames, sep='')

  checkCol(aCategory, main = names(clinical)[colId])
  matA
}

