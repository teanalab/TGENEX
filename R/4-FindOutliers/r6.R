#' version 1 (May/30/2018)

rm(list = ls())
k<-6
# LoadFiles --------
load("tempServer/bcCPk6.Rd", envir = globalenv())

#' Projection of the factor matrices
# projection -------
  A <- bcCP$A  #patients x k

  #remove outliers
  case1 <- which(A[,1] < (-5))
  length(case1)
  row.names(A)[case1] #"tcga-e2-a10c"
  A[case1,]

  case2 <- which(A[,3] < (-3))
  length(case2)
  row.names(A)[case2] #"tcga-an-a0xw" "tcga-bh-a0bz" "tcga-bh-a0hp"
  A[case2,]

  case3 <- which(A[,5] > 2)
  length(case3)
  row.names(A)[case3] #"tcga-ar-a0tx"
  A[case3,]

  case4 <- which(A[,4] < (-2))
  length(case4)
  row.names(A)[case4] #"tcga-bh-a0dz"
  A[case4,]

  A1 <- A[-unique(c(case1,case2,case3,case4)),]

  # Prepare Data
  mydata <- na.omit(A1) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  save(mydata,file=paste("temp/verotypes_k",k,".RData",sep=''))

  mydata <- na.omit(A) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  save(mydata,file=paste("temp/verotypesWithOutL_k",k,".RData",sep=''))

#
# #paartitioning
# # Determine number of clusters
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")
#
# # K-Means Cluster Analysis
# fit <- kmeans(mydata, 11) # 5 cluster solution
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster)
#
#
# affiliationList = mydata$fit.cluster
# names(affiliationList) = row.names(mydata)
#
# save(affiliationList,file="data/affiliationList_k6.RData")





