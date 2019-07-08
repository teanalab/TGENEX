#' version 1 (May/30/2018)

rm(list = ls())
k<-10

# LoadFiles --------
load(paste("tempServer/bcCPk",k,".Rd",sep=''), envir = globalenv())

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
  outliers <- which(A[,1] < (-5))
  length(outliers)
  row.names(A)[outliers] #"tcga-c8-a1hm"
  A <- A[-outliers,]

  outliers <- which(A[,3] < (-2))
  length(outliers)
  row.names(A)[outliers] #"tcga-e2-a10c"
  A <- A[-outliers,]

  outliers <- which(A[,4] < (-1.5))
  length(outliers)
  row.names(A)[outliers] #"tcga-a8-a09g" "tcga-bh-a0dz"
  A <- A[-outliers,]

  outliers <- which(A[,6] > (1))
  length(outliers)
  row.names(A)[outliers] # "tcga-a2-a0t0"
  A <- A[-outliers,]

  outliers <- which(A[,8] > (1))
  length(outliers)
  row.names(A)[outliers] #"tcga-bh-a0hp"
  A <- A[-outliers,]

  outliers <- which(A[,9] < (-1.5))
  length(outliers)
  row.names(A)[outliers] #"tcga-ar-a0tx"
  A <- A[-outliers,]

  outliers <- which(A[,10] < (-0.8))
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

  save(mydata,file=paste("temp/verotypes_k",k,".RData",sep=''))

  mydata <- na.omit(A) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  save(mydata,file=paste("temp/verotypesWithOutL_k",k,".RData",sep=''))


# #paartitioning
# # Determine number of clusters
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
#
# # K-Means Cluster Analysis
# fit <- kmeans(mydata,10) # 5 cluster solution
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster)
#
#
# affiliationList = mydata$fit.cluster
# names(affiliationList) = row.names(mydata)
#
# save(affiliationList,file=paste("data/affiliationList_k",k,".RData",sep=''))
#
# source('~/CLIGEN_tgit/R/6-SurvivalAnalysis/6-9v1-SurvivalAnalysisofSVNSmooth.R', echo=TRUE)
# survivalAnalysis()
#
#
#
