rm(list = ls())
k<-2

# LoadFiles --------
load(paste("dataExample/Filtered/bcCPk",k,".Rd",sep=''), envir = globalenv())

#' Projection of the factor matrices
# projection -------
A <- bcCP$A  #patients x k

hists <- function(A,k)
{
  for(i in seq(1:k)){
    hist(A[,i], breaks = 100, density = T)
    invisible(readline(prompt="Press [enter] to continue"))
  }
}

#hists(A,k)

backUpA <- A
backUpA -> A

#remove outliers
outliers <- which(A[,1] < (-4))
length(outliers)
row.names(A)[outliers] #"tcga-e2-a10c"
A <- A[-outliers,]

outliers <- which(A[,2] > (2.8))
length(outliers)
row.names(A)[outliers] #"tcga-a1-a0so"
A <- A[-outliers,]

A1 <- A
A <- backUpA
#hists(A1,k)

#A1<-A #do not remove outliers

# Prepare Data
mydata <- na.omit(A1) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables

save(mydata,file=paste("temp/subtypes_k",k,".RData",sep=''))

mydata <- na.omit(A) # listwise deletion of missing
mydata <- scale(mydata) # standardize variables

save(mydata,file=paste("temp/subtypesWithOutL_k",k,".RData",sep=''))
