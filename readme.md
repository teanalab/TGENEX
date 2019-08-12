---
# TGENEX: Tensor based integration of clinical data, gene mutation, and GEne EXpression
---

# Step 1 - Obtain the raw data

The input for the whole pipeline is the mutation file and clinical records. These files can be automatically downloaded using RTCGA or manually downloaded from firebrowse and cBiolite.

### Option 1 - Using RTCGA

The binary mutation table is automatically constructed but you need to manually select the clinical variables that you want to use in your analysis and use as input. We included the list of clinical variables that are relevant to Breast Cancer and Ovarian cancer as part of the '1-createBoolClinical.R' script.

See 'Example 1 - BRCA RTCGA.R' for an example on reading the files using RTCGA.


#### Clinical variables

List of biomarkers for BRCA that we used is:

1.    Demographics and clinical history:
•    Age
•    Sex
•    Race
•    BMI
•    Prior cancer history
•    Family history of cancer
•    Diagnosis of diabetes
•    Menstrual status

2. Histopathology results:
•    Tumor size
•    Histology designations
•    Tumor grade
•    Cancer stage
•    Lymph node stage
•    Metastasis stage

3. Test results for molecular markers or their immunohistochemistry
surrogates:
•    EGFR protein
•    Cytokeratin 5/6 (CK 5/6)
•    ER/PR
•    HER1
•    HER2
•    TP53
•    CA-125
•    Prostate-specific antigen (PSA)
•    KRAS
•    ERBB2
•    UGT-1A
•    EML4
•    ALK
•    BRCA

4. Clinical events related to cancer treatment, survival and outcomes:
•    Surgery
•    Chemotherapy (including prescribed medication)
•    Radiation therapy
•    Cancer recurrence
•    Patient death


-    Based on the list of biomarkers, the following clinical variables where selected:

Tumor.Coded
Node.Coded
Metastasis.Coded                                
estrogen_receptor_status
progesterone_receptor_status
HER2.Final.Status                                 
Histology             
Stage                     
age_at_initial_pathologic_diagnosis  
patient.ethnicity                          
patient.gender

### Option 2 - Downloading the files manually

1. Download the raw files in the folder 'rawData'.

2. Write a file named 'clinical_useful.csv' with the list of clinical variables that you want to include in your analysis. Locate the file in the folder 'rawData'. We provide an example of the clinical variables that are relevant to Breast Cancer.
This csv file should have the following column names:
"column", "USEFUL", "type". "column" refers to the name of the clinical variable in the raw file, "USEFUL" has a 'Y' or 'N' value indicating if the clinical variable is useful for the analysis or not, and "type" can have the values 'nominal' or 'continuos' indicating the type of useful variables.

3. See the file 'Example 2'. Here we use a gdca_Level1 clinical data (obtained from firebrowse) and a  CBioportal clinical data file, both located in the 'rawData' folder. The gdca_Level1 file is named "BRCA.clin.merged.txt" and the CBioportal file is named "data_clinical.txt". Here we loaded 505 patients and 6066 genes




-    Clinical variables:

Based on the list of biomarkers, the following clinical variables where selected:
Tumor.Coded
Node.Coded
Metastasis.Coded                                
estrogen_receptor_status
progesterone_receptor_status
HER2.Final.Status                                 
Histology             
Stage                     
ageOfDiagnosis  
patient.ethnicity                          
patient.gender


Then, these variables were dicomotized resulting in the 32 following boolean variables:
[1] "patient.ethnicityhispanic.or.latino"                          
 [2] "patient.genderfemale"                                         
 [3] "Tumor.CodedT_other"                                           
 [4] "Node.CodedPositive"                                           
 [5] "Metastasis.CodedPositive"                                     
 [6] "patient.breast_carcinoma_estrogen_receptor_statusnegative"    
 [7] "patient.breast_carcinoma_estrogen_receptor_statuspositive"    
 [8] "patient.breast_carcinoma_progesterone_receptor_statusnegative"
 [9] "patient.breast_carcinoma_progesterone_receptor_statuspositive"
[10] "HER2.Final.StatusEquivocal"                                   
[11] "HER2.Final.StatusNegative"                                    
[12] "HER2.Final.StatusPositive"                                    
[13] "X.clinical.patients..cols..infiltrating.ductal.carcinoma"     
[14] "X.clinical.patients..cols..infiltrating.lobular.carcinoma"    
[15] "X.clinical.patients..cols..mixed.histology..please.specify."  
[16] "X.clinical.patients..cols..other..specify"                    
[17] "X.clinical.patients..cols..No_Conversion"                     
[18] "X.clinical.patients..cols..Stage.I"                           
[19] "X.clinical.patients..cols..Stage.IIA"                         
[20] "X.clinical.patients..cols..Stage.IIB"                         
[21] "X.clinical.patients..cols..Stage.IIIA"                        
[22] "X.clinical.patients..cols..Stage.IIIC"                        
[23] "NA[26,36]"                                                    
[24] "NA(36,43]"                                                    
[25] "NA(43,49]"                                                    
[26] "NA(49,54]"                                                    
[27] "NA(54,59]"                                                    
[28] "NA(59,64]"                                                    
[29] "NA(64,70]"                                                    
[30] "NA(70,76]"                                                    
[31] "NA(76,82]"                                                    
[32] "NA(82,89]"  

#### Pruning genes


- Genes

We first obtained a histogram for clinical variables and genes in descending order to set threshold (frequency) to prune more aggressively the genes.

To produce these histograms run:

```{r}
#source('~/CLIGEN_tgit/R/1-TensorDefinition/1-1v3-descriptive_statitics_mutation.R')

"current number of genes:"
# [1] 11971
# [1] "X"
# [1] 1
# [1] "How many patients have less than or exactly to X mutations?"
# [1] 0
# [1] "How many genes appear mutated less than or exactly X times?"
# [1] 5905


2
8940

3
10352

4
11060

5
11432

6
11646

7
11756

8
11822

9
11867

10
11897


# out of 11971 genes


# X    | How many genes appear mutated less than or exactly X times?    | Genes left for analysis
0    0    11971
1    5905    6066
2    8940    3031
3    10352    1619
4    11060    911
5    11432    539
6    11646    325
7    11756    215
8    11822    149
9    11867    104
10    11897    74

```

[Fig3](../../temp/1-1v3-Fig3.pdf)

[Fig4](../../temp/1-1v3-Fig4.pdf)


Based on these figures, we conclude that the threshold is [[X = 5]]


After deciding the threshold we change this variable in the "0-ConfigurationVariables.R" file.



# Documentation

Prerequisites:

1. Folder data and temp need to be present.
2. Tested on:

    R version 3.5.0 (2018-04-23) -- "Joy in Playing"
    Copyright (C) 2018 The R Foundation for Statistical Computing
    Platform: x86_64-apple-darwin15.6.0 (64-bit)


# Step 2 - Construct the tensor


The script "2-BoolMatrices.R" generates two csv files: "data/10-binaC.csv" and "data/10-binaM.csv". These files are the binary matrices of clinical variables and mutation. '10-binaM.csv' has 2714 rows (genes) and 456 columns (patients). '10-binaC.csv' has 38 rows (clinical variables) and 456 columns (patients).
The input for this script is "boolMutation.RData"
and "data/boolClinical.RData". For details in the generation of these files, check Step 1.

### Smoothing the tensor

To smooth the tensor we used an SVD-based approach. First, we smooth the Mutation Matrix using SVD. The complete script and details can be found in "2-Smooth_MutationM_SVD.R". Then, we construct the tensor as the vector multiplication between the smoothed mutation matrix and the clinical matrix. For more details, check out the script "2-Construct_tensor.R".



# Step 3 - Tensor Factorization

Use the script "3-RunTF.R" to perform CADECOM/PARAFAC tensor factorization with a particular number-of-components parameter (e.g. we run for r= 2,3,4,5,6,7,8,9 and 10 components).
Change the value of R in line 30.


[//]: # (takes for ever
Do this to prepare Mac before running it
http://osxdaily.com/2012/10/11/mac-running-slow-reasons/
Keep the activity Monitor open
Close Every other program
Turn off wifi
Close RStudio, open with Atom, copy, close Atom, and run from console
Open terminal
cd Google\ Drive/current/CLINGEN/9source/CLIGEN
R)


# Step 4 - Obtaining subtyping

To obtain the subtypes we cluster the factors using using the components as features.
For the nine different factorizations (Components 2-10), we cluster using k-means and hierarchical clustering, with and without outliers.

We found the outliers of the nine factorizations (r = 2, 3, 4, 5, 6, 7, 8, 9, 10) and stored the subtypes. For all details see the scripts in the folder "R/4-FindOutliers". These scripts store the subtypes in the files: "verotypesWithOutL_kK.RData" and
"verotypes_kK.RData".


### k-means

To run k-means use the script
'4-kmeans.R' which will store the affiliation lists when including outliers and excluding outliers in the files "affiliationList_kK_iI.RData"
and
"affiliationListWithOutL_kK_iI.RData"
respectively.


## Hierarchical clustering:

To run hierarchical clustering use the script
'4-hclust.R'which will store the affiliation lists when including outliers and excluding outliers in the files
"affiliationHC_kK_iI.RData"
and
"affiliationHCwithOL_kK_iI.RData".


# Step 5,6,7 - Survival Analysis:

Run Survival Analysis on all the clustering and obtain a table with the Cox log rank test p-values using the script '5-getPvalsSurvAnal.R'.

To change the limit in the number of groups, change the value of the variable 'threshold_groups'.


To plot the Kaplan-meier curves of any clustering k, i run the script "6-survivalAnalysis.R" and comment out depending of the type of clustering you want to load.


'7-clinicalSubtypes.R' generates the subtypes using pam50.


# Step 8 - Obtain genes and clinical variables per cluster:

We find the distance between genes and clinical variables to each cluster centroid and ranked them by sum of squares. 



Diana Diaz.
PhD candidate.
Computer Science Department.
Wayne State University.
http://www.cs.wayne.edu/dmd
