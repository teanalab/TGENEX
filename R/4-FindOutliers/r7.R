#' version 1 (May/30/2018)
library("fpc")
library("NbClust")

#loadls("NbClust fpc")

rm(list = ls())
k<-7

# LoadFiles --------
load(paste("tempServer/bcCPk",k,".Rd",sep=''), envir = globalenv())

#' Projection of the factor matrices
# projection -------
  A <- bcCP$A  #patients x k

  hists <- function(A,k)
  {
    for(i in seq(1:k)){
      hist(A[,i], breaks = 100, density = T)
     # invisible(readline(prompt="Press [enter] to continue"))
    }
    par(ask=F)
  }

  #hists(A,k)


  #remove outliers
  case1 <- which(A[,3] < (-1))
  length(case1)
  row.names(A)[case1] #"tcga-e2-a10c"
  A[case1,]

  case2 <- which(A[,4] < (-1))
  length(case2)
  row.names(A)[case2] #"tcga-an-a0xw" "tcga-bh-a0bz" "tcga-bh-a0hp"
  A[case2,]


  case3 <- which(A[,5] > 2)
  length(case3)
  row.names(A)[case3] #"tcga-ar-a0tx"
  A[case3,]

  case4 <- which(A[,6] < (-2.5))
  length(case4)
  row.names(A)[case4] #"tcga-bh-a0dz"
  A[case4,]

  case5 <- which(A[,7] > (2))
  length(case5)
  row.names(A)[case5] #"tcga-bh-a0dz"
  A[case5,]



  A1 <- A[-unique(c(case1,case2,case3,case4,case5)),]

  # Prepare Data
  mydata <- na.omit(A1) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  k <- 7
  save(mydata,file=paste("temp/verotypes_k",k,".RData",sep=''))

  mydata <- na.omit(A) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  save(mydata,file=paste("temp/verotypesWithOutL_k",k,".RData",sep=''))

#
# #partitioning
#
# # Determine number of clusters ------
# # wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# # for (i in 2:15) wss[i] <- sum(kmeans(mydata,
# #                                      centers=i)$withinss)
# # plot(1:15, wss, type="b", xlab="Number of Clusters",
# #      ylab="Within groups sum of squares")
# #
# # pk1 <- pamk(mydata,krange=1:15,criterion="asw",critout=TRUE)
# # plot(pk1$crit)
# # pk1$nc #5
# # pk2 <- pamk(mydata,krange=1:15,criterion="multiasw",ns=2,critout=TRUE)
# # plot(pk2$crit)
# # pk2$nc #2
# # # "multiasw" is better for larger data sets, use larger ns then.
# # pk3 <- pamk(mydata,krange=1:15,criterion="ch",critout=TRUE)
# # plot(pk3$crit)
# # pk3$nc #2
# #
# # pka <- kmeansruns(mydata,krange=1:15,critout=TRUE,runs=2,criterion="asw")
# # plot(pka$crit)
# # pka$bestk #3
# #
# # pkc <- kmeansruns(mydata,krange=1:15,critout=TRUE,runs=2,criterion="ch")
# # plot(pkc$crit)
# # pkc$bestk #2
#
#
# save(mydata,file=paste("temp/verotypes_k",k,".RData",sep=''))
#
#
# # K-Means Cluster Analysis -----------
# set.seed(1234)
# fit <- kmeans(mydata,9, nstart = 25) # 5 cluster solution
#
#
#
# # get cluster means
# fit$centers
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
#
#
# #clusplot(mydata, fit$cluster, color=T, shade=T)
# plotcluster(mydata, fit$cluster)
#
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster)
#
#
#
# affiliationList = mydata$fit.cluster
# names(affiliationList) = row.names(mydata)
#
# save(affiliationList,file=paste("data/affiliationList_k",k,".RData",sep=''))
#
#
# source('~/CLIGEN_tgit/R/6-SurvivalAnalysis/6-9v1-SurvivalAnalysisofSVNSmooth.R', echo=TRUE)
# survivalAnalysis()
#
#
#
