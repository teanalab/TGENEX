#visualization of t-sne
library("tsne")
library("Rtsne")
library("rgl")


#load rubik factorization
patientsFactors <- read.table(file="MATLAB-RUBIK/rubikOutput/R_10/u2.csv",
                       sep = ",", quote = "'",
                       header = FALSE, stringsAsFactors = FALSE,
                       skipNul = FALSE, blank.lines.skip = FALSE)

#to reduce dimensionality, tsne runs out of memory
library(tsne)
tsne_res = tsne(patientsFactors, k=3, perplexity=50, epoch=50)
tsne_resD = as.data.frame(tsne_res)
names(tsne_res) <- c("V1",'V2','V3')

#plot the three dimensions
library(plotly)
p <- plot_ly(tsne_resD, x = ~V1, y = ~V2, z = ~V3 ) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'V1'),
                      yaxis = list(title = 'V2'),
                      zaxis = list(title = 'V3')))
p
  
  
  #setting up plotly credentials. run this one time only
  Sys.setenv("plotly_username"="datad")
  Sys.setenv("plotly_api_key"="HAbBEEkwutxoLYgL0TqM")
  
  # Create a shareable link to your chart
  #chart_link = api_create(p, filename="r10_allGenes")
  #chart_link
  #https://plot.ly/~datad/7/#/