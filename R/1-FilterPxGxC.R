
#' Title
#'
#' @param boolMutation 
#' @param boolClinical 
#' @param PlotFrequencies 
#'
#' @return
#' @export
#'
#' @examples rm(list=ls())
#' load("data/boolClinical.RData")
#' load("data/mutationSmooth.RData")
#' boolMutation=mutationSmooth
#' source('R/0-ConfigurationVariables.R')
#' list[boolMutation, boolClinical] <- filterGenesMutation(boolMutation, boolClinical, PlotFrequencies=T) #plot on temp
filterGenesMutation <- function(boolMutation, boolClinical, PlotFrequencies=T){
  
  sumGenes4patient <- unlist(apply(boolMutation,1,sum)) #sum genes across 430 patients
  sumClinVars4patient <- unlist(apply(boolClinical,1,sum)) #sum clinical variables across 430 patients
  sumPatients4gene <- unlist(apply(boolMutation,2,sum)) #sum patients across 20596 and find min
  sumPatients4clinVar <- unlist(apply(boolClinical,2,sum)) #sum patients accross 106 clinical variables
  
  minGenes = min(sumGenes4patient)  #all patients have atleast this number of genes
  minPatients4g = min(sumPatients4gene) #all genes are mutated in atleast this number of patients 
  minClinical = min(sumClinVars4patient) #all patients have atleast this number of clinical variables
  minPatients4c = min(sumPatients4clinVar) #all clinical variables are present in at least this number of patents

  
  #Plots before filtering ----
  if(PlotFrequencies)
  {
    require(gridExtra)
    
    fileName <- "temp/frequencies_beforeFiltering.pdf"
    pdf(fileName,  height=8, width=8, compress = FALSE)  
    
    ## genes across 430 patients
    dataP <- sumGenes4patient
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of mutations across ",sum(Atable$Freq)," patients")
    ylab = "Frequency (number of patients)" 
    xlab = "number of mutated genes"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " patient has ", Atable[dim(Atable)[1],1], " mutated genes.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
      
    }
    
    
    ## patients across 20596 genes
    dataP <- sumPatients4gene 
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of patients across ",sum(Atable$Freq)," genes")
    ylab = "Frequency (number of genes)" 
    xlab = "number of patients"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " gene is mutated in ", Atable[dim(Atable)[1],1], " patients.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    ## Clinical vars across 430 patients
    dataP <- sumClinVars4patient
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of clinical vars. across ",sum(Atable$Freq)," patients")
    ylab = "Frequency (number of patients)" 
    xlab = "number of clinical variables"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " patients have ", Atable[dim(Atable)[1],1], " clinical variables.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    ## patients accross 106 clinical variables
    dataP <- sumPatients4clinVar
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of patients across ", sum(Atable$Freq), " clinical variables")
    ylab = "Frequency (number of clinical v.)" 
    xlab = "number of patients"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " clinical variable is present in ", Atable[dim(Atable)[1],1], " patients.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    dev.off()
  }
  
  #repeat until convergence
  while (minGenes <= minNumOfMutationsPerPatient || minPatients4g <= minNumOfPatientsPerGene ||
         minPatients4c <= minNumOfPatientsPerClinical)
  {
    #Number of mutations per patient
    sumGenes4patient <- unlist(apply(boolMutation,1,sum)) #sum genes across 430 patients
    #Patients to remove
    PtoDel <- which(sumGenes4patient<=minNumOfMutationsPerPatient)
    if( length(PtoDel) > 0  ){
      boolMutation <- boolMutation[-PtoDel,]
      boolClinical <- boolClinical[-PtoDel,]
    }
    
    #Number of patients per clinical
    sumPatients4clinVar <- unlist(apply(boolClinical,2,sum)) #sum patients accross 106 clinical variables
    #Patients to remove
    CtoDel <- which(sumPatients4clinVar<=minNumOfMutationsPerPatient)
    if( length(CtoDel) > 0  ){
      boolClinical <- boolClinical[,-CtoDel]
    }
    

    #Number of patients per gene
    sumPatients4gene <- unlist(apply(boolMutation,2,sum)) #sum patients across 20596 and find min
    #Genes to remove
    ToDel <- which(sumPatients4gene<=minNumOfPatientsPerGene)
    if ( length(ToDel) > 0 ){
      boolMutation <- boolMutation[,-ToDel]
    }
    
    sumClinVars4patient <- unlist(apply(boolClinical,1,sum)) #sum clinical variables across 430 patients
    
    minGenes = min(sumGenes4patient)  #all patients have atleast this number of genes
    minPatients4g = min(sumPatients4gene) #all genes are mutated in atleast this number of patients 
    minClinical = min(sumClinVars4patient) #all patients have atleast this number of clinical variables
    minPatients4c = min(sumPatients4clinVar) #all clinical variables are present in at least this number of patents
  }
  
  #Plots after filtering ----
  if(PlotFrequencies)
  {
    require(gridExtra)
    
    fileName <- "temp/frequencies_afterFiltering.pdf"
    pdf(fileName,  height=8, width=8, compress = FALSE)  
    
    ## genes across 430 patients
    dataP <- sumGenes4patient
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of mutations across ",sum(Atable$Freq)," patients")
    ylab = "Frequency (number of patients)" 
    xlab = "number of mutated genes"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " patient has ", Atable[dim(Atable)[1],1], " mutated genes.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
      
    }
    
    
    ## patients across 20596 genes
    dataP <- sumPatients4gene 
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of patients across ",sum(Atable$Freq)," genes")
    ylab = "Frequency (number of genes)" 
    xlab = "number of patients"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " gene is mutated in ", Atable[dim(Atable)[1],1], " patients.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    ## Clinical vars across 430 patients
    dataP <- sumClinVars4patient
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of clinical vars. across ",sum(Atable$Freq)," patients")
    ylab = "Frequency (number of patients)" 
    xlab = "number of clinical variables"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " patients have ", Atable[dim(Atable)[1],1], " clinical variables.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    ## patients accross 106 clinical variables
    dataP <- sumPatients4clinVar
    Atable <- as.data.frame(table(dataP))
    main = paste0("Frequency of patients across ", sum(Atable$Freq), " clinical variables")
    ylab = "Frequency (number of clinical v.)" 
    xlab = "number of patients"
    exampleInter = paste0("For example, ", Atable[dim(Atable)[1],2], " clinical variable is present in ", Atable[dim(Atable)[1],1], " patients.")
    if (dim(Atable)[1]<=7){
      plot.new()
      title(main=main)
      grid.table(Atable)
      legend("top",exampleInter)#,cex=0.5)
    } else {
      hist(dataP, main = main ,
           ylab = ylab ,
           xlab = xlab,
           breaks = 100)
      legend("topright",exampleInter)#,cex=0.5)
    }
    
    dev.off()
  }
  list(boolMutation, boolClinical)
  
}
