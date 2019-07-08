#' version 1 (May/30/2018)

rm(list = ls())
# LoadFiles --------
  load("tempServer/bcCPk4.Rd", envir = globalenv())

#' Projection of the factor matrices
# projection -------

  A <- bcCP$A  #patients x k
  hist(A[,1], breaks = 100, density = T)
  hist(A[,2], breaks = 100, density = T)
  hist(A[,3], breaks = 100, density = T)
  hist(A[,4], breaks = 100, density = T)


  #remove outlier
  case1 <- which(A[,4] < (-4))
  length(case1)
  row.names(A)[case1] #"tcga-ar-a0tx"
  A[case1,]
  A1 <- A[-case1,]

  # Prepare Data
  mydata <- na.omit(A1) # listwise deletion of missing
  mydata <- scale(mydata) # standardize variables

  k <- 4
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
# fit <- kmeans(mydata, 5) # 5 cluster solution
# # get cluster means
# aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# mydata <- data.frame(mydata, fit$cluster)
#
#
#   affiliationList = mydata$fit.cluster
#   names(affiliationList) = row.names(mydata)
#
#
#   save(affiliationList,file="data/affiliationList_k4.RData")





