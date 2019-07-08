library(fpc)
library(NbClust)
library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(tidyr)
library(broom) #https://cran.r-project.org/web/packages/broom/vignettes/kmeans.html


rm(list = ls())
source("R/5-getCoxpVal.R")
source("R/8-writeRankedTables.R")

#kmeans r=6 k=8
#skemeans without OL k=13 r=9
r = 6
k = 8
threshold_groups = 10


inputFolder = "dataExample/Filtered" 
# inputFolder = "temp"

load(file=paste(inputFolder,"/verotypesWithOutL_k",r,".RData",sep=''))
#load(file=paste(inputFolder, "/verotypes_k",r,".RData",sep=''))


# K-Means Cluster Analysis -----------
set.seed(1234)
fit <- kmeans(mydata, centers = k , nstart = 25)
#aims to partition the points into k groups such that the sum of squares from points to the assigned cluster centres is minimized


# # Spherical Kmeans Cluster Analysis -----------
# set.seed(1234)
# fit <- skmeans(mydata, k) 



# Survival analysis to validate same clustering to table 
affiliationList = fit$cluster
survivalAnalysis(r, k, threshold_groups=threshold_groups)

#plotcluster(mydata, fit$cluster)
patientCenters <- fit$centers

# #summary per cluster level
# tidy(fit)
# 
# # plots
# kclusts <- tibble(k_ = 1:9) %>%
#   mutate(
#     kclust = map(k_, ~kmeans(mydata, .x)),
#     tidied = map(kclust, tidy),
#     glanced = map(kclust, glance),
#     augmented = map(kclust, augment, mydata)
#   )
# 
# clusters <- kclusts %>%
#   unnest(tidied)
# 
# assignments <- kclusts %>% 
#   unnest(augmented)
# 
# clusterings <- kclusts %>%
#   unnest(glanced, .drop = TRUE)
# 
# 
# p1 <- ggplot(assignments, aes(Comp.2, Comp.3)) +
#   geom_point(aes(color = .cluster)) + 
#   facet_wrap(~ k_)
# p1  #KmeansR6.pdf = Export -> save as pdf 7.49 x 6.29 inches



# # Determine number of clusters from: https://www.statmethods.net/advstats/cluster.html
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata, 
#                                      centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
#      ylab="Within groups sum of squares")




#testing distance function
fit$withinss
inGroup1 <- names(affiliationList[which(affiliationList==1)])
ss <- function(x) sum( scale(x, scale = FALSE)^2 )
ss(mydata[inGroup1,])
a_centroid = patientCenters[1,]
allDis <- c()
for ( i in inGroup1 ){
  a_point = mydata[i,]
  distance2centroid <- sum((a_centroid - a_point)^2 ) 
  allDis <- c(allDis, distance2centroid)
}
sum(allDis)

many_points = mydata[inGroup1,]
distance2centroid <- apply(many_points,1, function(x ) sum((a_centroid - x)^2 ) )
sum(distance2centroid)



# Get gene and clinical distance to centroids -----

load(paste("dataExample/Filtered/bcCPk",r,".Rd", sep=''), envir = globalenv())
patientsF <- bcCP$A  #p x k
genesF <- bcCP$B #genes x k
clinicalF <- bcCP$C #clinical x k.

# # traditional = non zero components
# ListComponents = projection(patientsF,genesF,clinicalF, k = r) 
# writeTopGenes(ListComponents, TopN = 200, r)
# writeTopClinical(ListComponents, TopN = 10, r)


# DISTANCE TO CENTROIDS
patientsF <- scale(patientsF) # standardize variables
genesF <- scale(genesF)
clinicalF <- scale(clinicalF)

distancesGenes = data.frame(matrix(nrow = length(row.names(genesF)), ncol = k) )
row.names(distancesGenes) = row.names(genesF)
names(distancesGenes) = as.character(paste0("Centroid",c(1:k))) 

for(i in seq(1,k))
  distancesGenes[,i] <- apply(genesF, 1, function(x) sum((patientCenters[i,] - x)^2 ) )

distancesClinical = data.frame(matrix(nrow = length(row.names(clinicalF)), ncol = k) )
row.names(distancesClinical) = row.names(clinicalF)
names(distancesClinical) = as.character(paste0("Centroid",c(1:k))) 

for(i in seq(1,k))
  distancesClinical[,i] <- apply(clinicalF, 1, function(x) sum((patientCenters[i,] - x)^2 ) )


writeTopGenesDC(distancesGenes, TopN = 200, k)
writeTopClinicalDC(distancesClinical, TopN = 10, k)
