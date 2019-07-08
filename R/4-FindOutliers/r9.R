rm(list = ls())
k<-9

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
outliers <- which(A[,1] < (-4.5))
length(outliers)
row.names(A)[outliers] #"tcga-c8-a1hm"
A <- A[-outliers,]

outliers <- which(A[,3] < (-1))
length(outliers)
row.names(A)[outliers] #"tcga-e2-a10c"
A <- A[-outliers,]

outliers <- which(A[,4] < (-1))
length(outliers)
row.names(A)[outliers] #"tcga-a8-a09g" "tcga-bh-a0dz"
A <- A[-outliers,]

outliers <- which(A[,5] > (2))
length(outliers)
row.names(A)[outliers] # "tcga-a2-a0t0"
A <- A[-outliers,]

outliers <- which(A[,6] < (-2))
length(outliers)
#row.names(A)[outliers] # "tcga-a2-a0t0"
#A <- A[-outliers,]

outliers <- which(A[,7] > (2))
length(outliers)
row.names(A)[outliers] #"tcga-bh-a0hp"
A <- A[-outliers,]

outliers <- which(A[,8] > (1))
length(outliers)
row.names(A)[outliers] #"tcga-ar-a0tx"
A <- A[-outliers,]

outliers <- which(A[,9] < (-1.5))
length(outliers)
row.names(A)[outliers] #"tcga-a8-a08l" "tcga-an-a0xw" "tcga-bh-a0bz"
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

