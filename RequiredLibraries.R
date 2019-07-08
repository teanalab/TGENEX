rm(list = ls())
#load my libraries
libs<-c("Packages.R")
libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
sapply(libs, function(u) {source(u)})


#' Loading my util functions
#'
#' @return
#' @export
#'
#' @examples
loadMyLib <- function(){
  #load my libraries
  libs<-c("Packages.R")
  libs<-paste("https://gist.githubusercontent.com/datad/39b9401a53e7f4b44e9bea4d584ac3e8/raw/", libs,sep='')
  sapply(libs, function(u) {source(u)})
  #rm(list=ls()[-which(ls() %in% c("patients", "mutation", "clinicalLong"))])
  rm(list=ls()[which(ls() %in% c("libs"))])

  save.image("data/loadls.RData")
}

loadlib("reporttools")
loadlib("gridExtra",FALSE) #to save tables as pdfs with function grid.table


loadlib("gsubfn",F) #to return multiple variables with list[,] need 0.7-0 or later
# myfun <- function() list(a = 1, b = 2)
# list[a, b] <- myfun()


# NMF ----
loadls("plyr survival lmtest", F)
loadls("stats registry", F)
loadls("methods utils pkgmaker registry rngtools cluster", F)
loadlib("NMF", F) #For non-negative matrix factorization
install.extras('NMF')
#The paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
#The vignette: https://cran.r-project.org/web/packages/NMF/index.html



#plot survival analysis
loadls("plyr missForest survAUC prodlim survminer", F)
require(purrr)
require(cluster)
require(survival)


#k-means
loadls("fpc NbClust", F)
require("fpc")
require("NbClust")


load("data/mutationSmoothPositive.RData", envir = globalenv())


#RTCGA mutation
loadls("RTCGA.mutations dplyr reshape2",F)

#RTCGA clinical and survival
loadls("RTCGA.clinical dplyr reshape2",F)
require(RTCGA.clinical)
require(RTCGA.clinical)
require(dplyr)
