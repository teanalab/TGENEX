library("skmeans")
rm(list = ls())

inputFolder = "dataExample/Filtered" 
# inputFolder = "temp"

# Spherical Kmeans Cluster Analysis -----------
for(r in c(3:10))
{
  load(file=paste(inputFolder,"/verotypesWithOutL_k",r,".RData",sep=''))
  #load(file=paste(inputFolder, "/verotypes_k",r,".RData",sep=''))
  for(i in c(2:15))
  {
    set.seed(1234)
    fit <- skmeans(mydata, i) 
    # plotcluster(mydata, fit$cluster)
    # append cluster assignment
    affiliationList = fit$cluster
    save(affiliationList,file=paste("temp/affiliationListWithOutL_sk",r,"_i",i,".RData",sep=''))
    #save(affiliationList,file=paste("temp/affiliationListSp_r",r,"_i",i,".RData",sep=''))
    rm(list = ls()[-which(ls() %in% c("inputFolder", "i", "r", "mydata")  )])
  }
  rm(list = ls()[-which(ls() %in% c("inputFolder", "i","r")  )])
}


