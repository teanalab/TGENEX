library("NbClust")

rm(list = ls())

#prefix = "dataExample/Filtered/verotypesWithOutL_k"
#prefix = "dataExample/Filtered/verotypes_k"
prefix = "temp/patientFactors_r"

rRange = c(3:78)
clustersK <- c(2:15)


# K-Means Cluster Analysis -----------
for(r in rRange)
{
  for(k in clustersK )
  {
    load(file=paste(prefix, r,".RData",sep=''))
    set.seed(1234)
    fit <- kmeans(mydata, centers = k , nstart = 25) 
    # get cluster means
    #fit$centers
    fwnss <- fit$withinss
    #aggregate(mydata,by=list(fit$cluster),FUN=mean)
    #clusplot(mydata, fit$cluster, color=T, shade=T)
    #plotcluster(mydata, fit$cluster)
    affiliationList = fit$cluster
    
    # remove outliers - do not use this
    # singlePoints <- which(fwnss == 0)
    # while (length(singlePoints) > 0 ) {
    #   toRemove <- which( affiliationList %in% singlePoints )
    #   mydata <- mydata[-toRemove,]
    #   set.seed(1234)
    #   fit <- kmeans(mydata, centers = k , nstart = 25) # 5 cluster solution
    #   fwnss <- fit$withinss
    #   affiliationList = fit$cluster
    #   singlePoints <- which(fwnss == 0)
    # } 
    
    save(affiliationList,file=paste("temp/affiListKmeans_r",r,"_k",k,".RData",sep=''))
    #save(affiliationList,file=paste("temp/affiliationListWithOutL_k",r,"_k",k,".RData",sep=''))
    #save(affiliationList,file=paste("temp/affiliationList_k",r,"_k",k,".RData",sep=''))
    rm(list = ls()[-which(ls() %in% c("prefix","k", "r", "mydata", "rRange", "clustersK")  )])
  }
  rm(list = ls()[-which(ls() %in% c("prefix","k","r","rRange", "clustersK")  )])
}

