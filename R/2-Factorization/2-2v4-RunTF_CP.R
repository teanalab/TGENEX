#version 4 (September-4-2018)
rm(list = ls())

#Get all needed data
LoadMyData <- function()
{
  load("data/binaClinicalSmooth.Rd", envir = globalenv())
  load("data/matrixMutationSmooth.Rd", envir = globalenv())
  load("data/MatrixSmoothOrderPxGxC.RData", envir = globalenv())
  load("data/loadls.RData", envir = globalenv())

  NBSgenes <- read.table(file="/Users/diamac/GitLab/NBS_cligen/nbs_release_v02_wc/output/ova_genes_raw.csv",
                         sep = "\t", quote = "'", header = FALSE, stringsAsFactors = FALSE)

  NBSgenes<-unlist(NBSgenes)
  NBSgenesKeys <- as.numeric(unlist(names(binaMutation)))
  labelsGenes =  NBSgenes[NBSgenesKeys]
  labelsPA  = row.names(binaMutation)
  labelsClinical = names(binaClinical)

  loadls("ThreeWay") #For tensors
}
LoadMyData()

runCP <- function(MatrixOrderPxGxC)
{


#takes for ever
#Do this to prepare Mac before running it
#http://osxdaily.com/2012/10/11/mac-running-slow-reasons/
#Keep the activity Monitor open
#Close Every other program
#Turn off wifi
#Close RStudio, open with Atom, copy, close Atom, and run from console
#Open terminal
#cd Google\ Drive/current/CLINGEN/9source/CLIGEN
#R

#  Boolcp_ANOVA <- CP(MatrixOrderPxGxC, labelsPA, labelsGenes, labelsClinical)
#  507
#  12003
#  89
#  1



#careful! if somethings goes wrong I will
#lose all the time
#   Boolcp_3comp <- CP(MatrixOrderPxGxC, labelsPA, labelsGenes, labelsClinical)
#   503
#   11996
#   70
#   0
#   0
#   0
#   0
#   0
#   10
#   0
#   1e-6
#   1
#   10000
#   0
#   0
#   0
#   0
#   0
#   0
#   0
#   0
#
#
#   503 #number of A-mode entities
#   11996 #number of B-mode entities
#   70 #number of C-mode entities
#   0 #no ANOVA
#   0 #no Centered
#   0 #no normalized
#   0 #no PCA
#   0 #no scree test for deciding
#   10 #number of components
#   0 #unconstrained
#   1e-6 #convergence criterion
#   4 #additional runs
#   10000 #maximum number of iterations
#   time: 9:45pm
#   until: 9:00am
#   gave a fit of  8.5 %
#   Run no. 1
# Candecomp/Parafac function value at Start is  329091.729720201
# Candecomp/Parafac function value is 301792.413478017 after 11 iterations
# Fit percentage is 8.33386584514857 %
# Procedure used 8260.61 seconds
# Run no. 2
# Candecomp/Parafac function value at Start is  329240.124704498
# f= 303544.014325045 after 50 iters; diff.= 173.051829638542
# Candecomp/Parafac function value is 302817.105390185 after 60 iterations
# Fit percentage is 8.02262692033375 %
# Procedure used 4402.9 seconds
# Run no. 3
# Candecomp/Parafac function value at Start is  329240.456246504
# Candecomp/Parafac function value is 301230.82982118 after 16 iterations
# Fit percentage is 8.50444071889547 %
# Procedure used 1179.88 seconds
# Run no. 4
# Candecomp/Parafac function value at Start is  329240.186917669
# Candecomp/Parafac function value is 301248.740330289 after 42 iterations
# Fit percentage is 8.499000598278 %
# Procedure used 3303.74 seconds
# Run no. 5
# Candecomp/Parafac function value at Start is  329239.815783814
# Candecomp/Parafac function value is 302222.32023851 after 40 iterations
# Fit percentage is 8.20328638383207 %
# Procedure used 3079.3 seconds
#
# Fit (%) values from all runs:
# Start n.1 Start n.2 Start n.3 Start n.4 Start n.5
#      8.33      8.02      8.50      8.50      8.20
#
# Candecomp/Parafac analysis with  10  components, gave a fit of  8.5 %
# Simple check on degeneracy: inspect matrix of triple congruences
#         Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7 Comp.8 Comp.9 Comp.10
# Comp.1   1e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00 -1e-04   0e+00
# Comp.2   0e+00  1e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00 -3e-04   0e+00
# Comp.3   0e+00  0e+00  1e+00  0e+00  0e+00  0e+00  0e+00  0e+00 -7e-04   0e+00
# Comp.4   0e+00  0e+00  0e+00  1e+00  0e+00  0e+00  0e+00  0e+00 -4e-04   0e+00
# Comp.5   0e+00  0e+00  0e+00  0e+00  1e+00  0e+00  0e+00  0e+00 -5e-04   0e+00
# Comp.6   0e+00  0e+00  0e+00  0e+00  0e+00  1e+00  0e+00  0e+00 -1e-04   0e+00
# Comp.7   0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  1e+00  0e+00 -2e-04   0e+00
# Comp.8   0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  1e+00 -3e-04   0e+00
# Comp.9  -1e-04 -3e-04 -7e-04 -4e-04 -5e-04 -1e-04 -2e-04 -3e-04  1e+00  -2e-04
# Comp.10  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00  0e+00 -2e-04   1e+00
#   1 #scale solution in patients
#   #It is sometimes useful to SCALE solution, e.g., scale two matrices so that
#   #they have unit sum of squares compensating this scale in remaining matrix.
#   #If you want to scale components, specify '1':
#   0 #no permute or reflec
#   1 #stability check
#   0 #random sample
#
#
#   You can now manually PERMUTE and REFLECT columns of solution
#   If you want to reflect/permute columns, specify '1':
#   1:
#   Read 0 items
#
#   If you want to carry out a STABILITY CHECK on current or different solution, specify '1':
#   1:
#   Read 0 items
#
#   If you want to carry out a BOOTSTRAP procedure for
#   computing confidence intervals for the current solution, specify '1':
#   1:
#   Read 0 items
#
#   If you want to carry out a FITPARTITIONING on current solution, specify '1':
#   1:
#   Read 0 items
#
#   Press 'return' to conclude the analysis
#   1:
#   Read 0 items
#   The component matrices are normalized such that A and B have unit sums of squares
#   To see component matrices and fit and other results, digits: $Xprep, $A, $B, $C, $fit, $fitA, $fitB, $fitC, ...

}

save.image(file="temp/31-Boolcp_3comp_V2.RData")
