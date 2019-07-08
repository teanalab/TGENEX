library("purrr")
library("cluster")
library("fpc")
library("NbClust")

rm(list = ls())

#prefix = "dataExample/Filtered/verotypesWithOutL_k"
#prefix = "dataExample/Filtered/verotypes_k"
prefix = "temp/patientFactors_r"

rRange = c(10,20,30,40)
clustersK <- c(2:15)

# K-Means Cluster Analysis -----------
for(r in rRange)
{
  #withoutOutliers
  #load(file=paste("temp/verotypes_k",k,".RData",sep=''))
  #with ourliers
  #load(file=paste("temp/verotypesWithOutL_k",k,".RData",sep=''))
  load(file=paste(prefix, r,".RData",sep=''))
  set.seed(1234)
  

  # Dissimilarity matrix
  d <- dist(mydata, method = "euclidean")
  # # Hierarchical clustering using Complete Linkage
  # Ward's method
  hc5 <- hclust(d, method = "ward.D2" )


  for(k in clustersK )
  {
    # Cut tree into  groups
    sub_grp <- cutree(hc5, k)
    # Number of members in each cluster
    #table(sub_grp)
    affiliationList = sub_grp
    #save(affiliationList,file=paste("temp/affiliationHC_k",k,"_i",i,".RData",sep=''))
    #save(affiliationList,file=paste("temp/affiliationHCwithOL_r",r,"_k",k,".RData",sep=''))
    save(affiliationList,file=paste("temp/affiListHierClu_r",r,"_k",k,".RData",sep=''))
    rm(list = ls()[-which(ls() %in% c("prefix","k","r","rRange", "clustersK", "hc5" )  )])
  }
  rm(list = ls()[-which(ls() %in% c("prefix","k","r","rRange", "clustersK")  )])
}


  # # Hierarchical clustering using Complete Linkage
  # hc1 <- hclust(d, method = "complete" )
  # # Plot the obtained dendrogram
  # plot(hc1, cex = 0.6, hang = -1)
  # # Compute with agnes
  # hc2 <- agnes(mydata, method = "complete")
  # # Agglomerative coefficient
  # hc2$ac
  # # methods to assess
  # m <- c( "average", "single", "complete", "ward")
  # names(m) <- c( "average", "single", "complete", "ward")
  # # function to compute coefficient
  # ac <- function(x) {
  #   agnes(mydata, method = x)$ac
  # }
  # map_dbl(m, ac)
  # hc3 <- agnes(mydata, method = "ward")
  # pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")
  # Ward's method
  