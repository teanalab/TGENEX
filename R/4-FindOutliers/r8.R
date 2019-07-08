#' version 1 (May/30/2018)

rm(list = ls())
k<-8

# LoadFiles --------
load(paste("tempServer/bcCPk",k,".Rd",sep=''), envir = globalenv())

#' Projection of the factor matrices
# projection -------
  A <- bcCP$A  #patients x k

  hists <- function(A,k)
  {
    for(i in seq(1:k)){
      hist(A[,i], breaks = 100, density = T)
      #invisible(readline(prompt="Press [enter] to continue"))
    }
  }

  #hists(A,k)


  #remove outliers
  case1 <- which(A[,3] < (-1))
  length(case1)
  row.names(A)[case1] #"tcga-e2-a10c"
  A[case1,]

  case2 <- which(A[,4] < (-1))
  length(case2)
  row.names(A)[case2] #"tcga-a8-a09g" "tcga-bh-a0dz"
  A[case2,]


  case3 <- which(A[,5] > 1)
  length(case3)
  row.names(A)[case3] #"tcga-a1-a0so" "tcga-a2-a0t0" "tcga-a8-a07c"
  A[case3,]

  case4 <- which(A[,6] < (-2))
  length(case4)
  row.names(A)[case4] # "tcga-an-a0xw" "tcga-bh-a0bz"
  A[case4,]

  case5 <- which(A[,7] > (1))
  length(case5)
  row.names(A)[case5] #"tcga-a8-a09w" "tcga-bh-a0hp"
  A[case5,]

  case6 <- which(A[,8] > (1))
  length(case6)
  row.names(A)[case6] #"tcga-ar-a0tx"
  A[case6,]


  A1 <- A[-unique(c(case1,case2,case3,case4,case5,case6)),]


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
# fit <- kmeans(mydata,9) # 5 cluster solution
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




