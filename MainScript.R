## ---------------------------
##
## This script was created with R version 4.2.1 (2022-06-23)
##
##
## Title: Patient-specific identification of genome-wide DNA-methylation differences between intracranial and extracranial melanoma metastases
##
## Purpose of script: Script to rebuild analysis of manuscript
##
## Author: Theresa Kraft
##
## Date created: 08-15-2022
## Date last modified: 09-13-2022
##
## R-code and data usage license: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
##
## Copyright (c) Theresa Kraft, 2022
## Email: theresa.kraft@tu-dresden.de
##
## ---------------------------
##
## Notes: - Please set your directories before you start. 
##        - All required R-functions can be found in the R-script HelperFunctions.R.
##
## ---------------------------


################ 
## DIRECTORIES 
################ 

# workingDirectory
workingDirectory <- "~YOUR DIRECTORY~"

# This directory is where data for figures are collected and where lateron, figures are created separately 
figureDirectory <- "~YOUR DIRECTORY~/Figures/"

# This directory is where the methylation data is stored
# please download data from Zenodo 10.5281/zenodo.7074484
dataDirectory <- "~YOUR DIRECTORY~/Data/"

# This directory is where all annotation files are stored
annotationDirectory <- "~YOUR DIRECTORY~/Annotations/"

# This directory is where you can find the ARHMM.jar file 
hmmDirectory <- "~YOUR DIRECTORY~/Programs/"


setwd(workingDirectory)


################ 
## LIBRARIES 
################ 

library(stringr)
library(pvclust)
library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
library(pbapply)
library(plyr)
library(clusterProfiler)

################ 
## START 
################ 


# read in methylation data file 
methylationSamples <- read.csv(file=paste0(dataDirectory,
              "Table_methylationSamples.csv"),
              sep = ",", check.names = F) 

################ 
## FIGURE 1
# stability dendrogram
################ 

# create data for dendrogram 
methylationValues <- methylationSamples[,7:43]


# use pvclust with 1000 bootstraping repetitions for stability analysis
# parallel = T or define parallel = as.integer(numCores) if parallel computing is desired
# took approx. 25hours in a parallel mode with 50 cores

t1 <- Sys.time()
valuesManhD2.pv <- pvclust(methylationValues, method.hclust = "ward.D2", method.dist = "manhattan", parallel = as.integer(50), nboot = 1000)

print(difftime(Sys.time(), t1) )
timePassed <- difftime(Sys.time(), t1)

outputFileFig1 <- paste0(figureDirectory, "Figure1-Stability-Dendrogram/stabilityHclust_wardD2_manhattan.RData" ) 
save(valuesManhD2.pv, file = outputFileFig1)


################ 
## Metastases pairs 
## creation 
################ 

# create metastases pairs: Intracranial methylation value minus extacranial methylation value 
methylationPairs <- data.frame("ID" = methylationSamples$ID, "Chr" = methylationSamples$Chr, "Pos" = methylationSamples$Pos, "Gene" = methylationSamples$Gene, 
                               "RegFeature" = methylationSamples$RegFeature, "RegFeaturePos" = methylationSamples$RegFeaturePos, 
                               "P42-BLym-1" = methylationSamples$P42_B - methylationSamples$`P42_Lym-1`, "P42-BLym-2" = methylationSamples$P42_B  - methylationSamples$`P42_Lym-2`,
                               "P39-BLun" = methylationSamples$P39_B  - methylationSamples$P39_Lun, "P28-BSki" = methylationSamples$P28_B - methylationSamples$P28_Ski,
                               "P18-BLun-1" = methylationSamples$P18_B - methylationSamples$`P18_Lun-1`, "P18-BLun-2" = methylationSamples$P18_B - methylationSamples$`P18_Lun-2`,
                               "P17-BLym-1" = methylationSamples$P17_B - methylationSamples$`P17_Lym-1`, "P17-BLym-2" = methylationSamples$P17_B - methylationSamples$`P17_Lym-2`,
                               "P06-BLym-1" = methylationSamples$`P06_B-1` - methylationSamples$P06_Lym, "P06-BLym-2" = methylationSamples$`P06_B-2` - methylationSamples$P06_Lym, 
                               "P04-BSki-1" = methylationSamples$P04_B - methylationSamples$`P04_Ski-1`, "P04-BSki-2" = methylationSamples$P04_B - methylationSamples$`P04_Ski-2`,
                               "P03-BLun" = methylationSamples$P03_B - methylationSamples$P03_Lun, "P02-BLym" = methylationSamples$P02_B - methylationSamples$P02_Lym,    
                               "P09-BLiv" = methylationSamples$P09_B - methylationSamples$P09_Liv, "P08-BSof-1" = methylationSamples$P08_B - methylationSamples$`P08_Sof-1`,    
                               "P08-BSof-2" = methylationSamples$P08_B - methylationSamples$`P08_Sof-2`, "P08-BSof-3" = methylationSamples$P08_B - methylationSamples$`P08_Sof-3`, 
                               "P67-BLun-1" = methylationSamples$`P67_B-1` - methylationSamples$`P67_Lun-1`,   "P67-BLun-2" = methylationSamples$`P67_B-1` - methylationSamples$`P67_Lun-2`,
                               "P67-BLun-3" = methylationSamples$`P67_B-2` -  methylationSamples$`P67_Lun-1`,   "P67-BLun-4" = methylationSamples$`P67_B-2` - methylationSamples$`P67_Lun-2`,
                               "P64-BLun" =   methylationSamples$P64_B - methylationSamples$P64_Lun, "P16-BLun" = methylationSamples$P16_B - methylationSamples$P16_Lun,check.names = F     )

################ 
## FIGURE 2 
# autocorrelation
################ 

# calculate autocorrelation for 1000 random CpG orders
set.seed(1)
rows <- sample(nrow(methylationPairs))
permutatedDF <- methylationPairs[rows, ]

allAcfPermutated <- create_methyl_acf_df(permutatedDF)

allAcfPermutated <- allAcfPermutated[,c(1,5)]

for(i in 2:1000){
  set.seed(i)
  rows <- sample(nrow(methylationPairs))
  permutatedDF <- methylationPairs[rows, ]
  
  allAcfPermutated2 <- create_methyl_acf_df(permutatedDF)
  
  allAcfPermutated <- cbind(allAcfPermutated, allAcfPermutated2$quant_50 )
}

allAcfPermutatedCalc <- allAcfPermutated[,-1]
neededQuantiles = c(0 ,0.1 ,0.25 ,0.5 ,0.75, 0.9, 1 )
allMethAcfsQuantiles <- NULL

q = t(sapply(1:nrow(allAcfPermutatedCalc),FUN = function(lag) quantile(allAcfPermutatedCalc[lag,],neededQuantiles)))
colnames(q) = gsub("\\%","",paste0("quant_",colnames(q)))
allMethAcfsQuantiles = data.frame(lag = 0:100, as.data.frame(q))

# calculate acf for biological data 
allMethAcfQuantilesChrCondensed <- create_methyl_acf_df(methylationPairs)

autocorrelationData <- list(allMethAcfQuantilesChrCondensed,allMethAcfsQuantiles)

outputFileFig2 <- paste0(figureDirectory, "Figure2-Autocorrelation/autocorrelationData.RData" ) 
save(autocorrelationData, file = outputFileFig2)



################ 
## HMM training 
################ 


# create training data in the desired format for ARHMM
hmmTrainingData <- data.frame("ID" = methylationPairs$ID, "chr" = methylationPairs$Chr, 
                              "position" = methylationPairs$Pos, methylationPairs[,7:30],check.names = F)


outputFileHMMTrain <- paste0(workingDirectory, "DNA_Methylation_TrainingData.txt")
write.table( hmmTrainingData, file = outputFileHMMTrain, row.names = FALSE,
             col.names = TRUE, dec = ".", sep = "\t", quote = FALSE )

# bash script for HMM training using optimal parameters
# will create the StatePosteriorDecoding file that contains the most probable state for each CpG in each metastases pair
# ARHMM_Trainer.jar from Seifert et. al, 2014 (PMID: 24955771) 
bash <- paste0("java -jar ",hmmDirectory,"ARHMM_Trainer.jar -hmmOrder 1 -arOrder 0 -initialStateDist 0.1 0.8 0.1 -means -3 0 3 -sds 0.3 0.5 0.3 -scaleMeans 1e+06 1000 1e+06 -shapeSds 5e+05 10 5e+05 -scaleSds 1E-4 1E-4 1E-4 -statePosterior -modelBasics -model model1 -dataSet DNA_Methylation_TrainingData.txt")
system(bash)


methylationClassification <- read.delim(paste0(workingDirectory, 
                                               "model1_DNA_Methylation_TrainingData.txt_StatePosteriorDecoding.txt"))



################ 
## Create big table including HMM predictions  
################ 

# merge methylation classification and methylation values to one dataframe
methylationPairsRun <- merge(methylationClassification, methylationPairs, by.x = "Genes",  by.y = "ID" )

colnames(methylationPairsRun)[1] <- "ID"



################ 
## Annotation  
################ 


### add CpG island/shore/shelf annotation
IslandAnnotation <- as.data.frame(Islands.UCSC)
methylationPairs.annotated <- merge(methylationPairsRun, IslandAnnotation, by.x = "ID",  by.y = 0 )




### add pathways annotation 
# from Seifert et. al 2016 (PMID: 27716417) 
pathwayPath <- paste0(annotationDirectory, "geneset_with_pathways_merged.txt")
pathways <- read.table( file=pathwayPath, header=T, check.names = FALSE)


pathwayArray <- c()
for(i in 1:nrow(pathways)){
  gene <- as.character(pathways[i,7])
  path <- ""
  for(j in 8:32){
    if(as.character(pathways[i,j]) == "1"){
      path <- paste0(path, colnames(pathways)[j],";")
    }
  }
  pathwayArray <- rbind(pathwayArray, c(gene, path))
}
pathwayArray <- as.data.frame(pathwayArray)

annodf <- pbapply(methylationPairs.annotated, FUN =  function(x) extract_paths(x), MARGIN = 1)

methylationPairs.annotated <- data.frame(methylationPairs.annotated, "pathways" = annodf)




### Enhancer / Gene body / Promoter annotation
enhancerAnnotPath <- paste0(annotationDirectory, "enhancerAnnotation.txt") 
enhancerAnnot <- read.table( file=enhancerAnnotPath, header=TRUE, sep = "\t", check.names = FALSE)

methylationPairs.annotated <- merge(methylationPairs.annotated, enhancerAnnot, by.x = "ID", by.y = "ID")


methylationPairs.annotated$GenomicRegion <- ""
# define all CpGs that are associated with genes to a gene body
methylationPairs.annotated$GenomicRegion[str_detect(methylationPairs.annotated$Gene, ".")] <- "body"

# might overwrite body in case of Promoter / Enhancer annotation
methylationPairs.annotated$GenomicRegion[str_detect(methylationPairs.annotated$Enhancer, ".") & 
                                           !(str_detect(methylationPairs.annotated$RegFeature, "^Promoter"))] <- "enhancer"
methylationPairs.annotated$GenomicRegion[!str_detect(methylationPairs.annotated$Enhancer, ".") & 
                                           (str_detect(methylationPairs.annotated$RegFeature, "^Promoter"))] <- "promoter"
methylationPairs.annotated$GenomicRegion[str_detect(methylationPairs.annotated$Enhancer, ".") & 
                                           (str_detect(methylationPairs.annotated$RegFeature, "^Promoter"))] <- "enhancer;promoter"

methylationPairs.annotated$Enhancer <- NULL

# change colnames, metastases pairs *_p = methylation prediction
colnames(methylationPairs.annotated) <-  c("ID","P42.BLym.1_p","P42.BLym.2_p","P39.BLun_p","P28.BSki_p","P18.BLun.1_p","P18.BLun.2_p","P17.BLym.1_p",
                                           "P17.BLym.2_p","P06.BLym.1_p","P06.BLym.2_p","P04.BSki.1_p","P04.BSki.2_p","P03.BLun_p","P02.BLym_p","P09.BLiv_p",
                                           "P08.BSof.1_p","P08.BSof.2_p","P08.BSof.3_p","P67.BLun.1_p","P67.BLun.2_p","P67.BLun.3_p","P67.BLun.4_p",
                                           "P64.BLun_p","P16.BLun_p","Chr","Pos","Gene","RegFeature","RegFeaturePos","P42.BLym.1","P42.BLym.2",
                                           "P39.BLun","P28.BSki","P18.BLun.1","P18.BLun.2","P17.BLym.1","P17.BLym.2","P06.BLym.1",
                                           "P06.BLym.2","P04.BSki.1","P04.BSki.2","P03.BLun","P02.BLym","P09.BLiv","P08.BSof.1",
                                           "P08.BSof.2","P08.BSof.3","P67.BLun.1","P67.BLun.2","P67.BLun.3","P67.BLun.4",
                                           "P64.BLun","P16.BLun","Islands_Name","Relation_to_Island","pathways",
                                           "GenomicRegion")

################ 
## FIGURE 3 
# genomic regions
################ 

###
## Figure 3A 

# mapping of samples names "predictions" = c("values")

sampleNameDF <- data.frame("P03.BLun_p" = c("P03.BLun"), "P16.BLun_p" = c("P16.BLun"), "P18.BLun.1_p" = c("P18.BLun.1"), 
                           "P18.BLun.2_p" = c("P18.BLun.2"), "P39.BLun_p" = c("P39.BLun"), "P08.BSof.1_p" = c("P08.BSof.1"), 
                           "P08.BSof.2_p" = c("P08.BSof.2"), "P08.BSof.3_p" = c("P08.BSof.3"), "P42.BLym.1_p"  = c("P42.BLym.1"),
                           "P42.BLym.2_p" = c("P42.BLym.2"), "P04.BSki.1_p" = c("P04.BSki.1"), "P28.BSki_p" = c("P28.BSki"), 
                           "P17.BLym.1_p" = c("P17.BLym.1"),  "P17.BLym.2_p" = c("P17.BLym.2"), "P06.BLym.1_p" = c("P06.BLym.1"), 
                           "P06.BLym.2_p" = c("P06.BLym.2"), "P04.BSki.2_p" = c("P04.BSki.2"), "P02.BLym_p" = c("P02.BLym"), 
                           "P09.BLiv_p" = c("P09.BLiv"), "P67.BLun.1_p" = c("P67.BLun.1"), "P67.BLun.2_p" = c("P67.BLun.2"), 
                           "P67.BLun.3_p" = c("P67.BLun.3"), "P67.BLun.4_p" = c("P67.BLun.4"), "P64.BLun_p" = c("P64.BLun"),
                           stringsAsFactors = FALSE)



###
## calculate number of CpGs +/- of all metastases pairs 

PatientID <- c()
Methylation <- c()
value <- c()
overviewDF <- NULL
for(i in 1:length(sampleNameDF)){
  sampleValues <- as.character(sampleNameDF[i])
  samplePredictions <- as.character(colnames(sampleNameDF)[i])
  
  
  highMethyl <- sum(methylationPairs.annotated[,samplePredictions] == "+")
  lowMethyl <- sum(methylationPairs.annotated[,samplePredictions] == "-")
  
  PatientID <- c(PatientID, sampleValues)
  Methylation <- c(Methylation, "Decreased Methylation")
  value <- c(value, lowMethyl)
  
  PatientID <- c(PatientID, sampleValues )
  Methylation <- c(Methylation, "Increased Methylation")
  value <- c(value, highMethyl)
}


ranking <- c()
for(i in 1:length(PatientID)){
  if(str_detect(PatientID[i], "Lun")){
    ranking <- c(ranking,5)
  }
  if(str_detect(PatientID[i], "Lym")){
    ranking <- c(ranking,4)
  }
  if(str_detect(PatientID[i], "Ski")){
    ranking <- c(ranking,3)
  }
  if(str_detect(PatientID[i], "Liv")){
    ranking <- c(ranking,1)
  }
  if(str_detect(PatientID[i], "Sof")){
    ranking <- c(ranking,2)
  }
}

barplotDF <- data.frame(PatientID, Methylation, value, ranking)

for(i in 1:length(barplotDF[,1])){
  if(as.character(barplotDF[i,2]) == "Decreased Methylation"){
    barplotDF[i,3] <- - barplotDF[i,3]
  }
}
barplotOverviewDF <-  barplotDF


outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/barplotOverviewDF.RData" ) 
save(barplotOverviewDF, file = outputFileFig3)




###
## Figure 3B

promAnnotDF <- methylationPairs.annotated[str_detect(methylationPairs.annotated$GenomicRegion, "promoter"),]
enhancerDF <- methylationPairs.annotated[str_detect(methylationPairs.annotated$GenomicRegion, "enhancer"),]
geneAnnotDF <- methylationPairs.annotated[str_detect(methylationPairs.annotated$GenomicRegion, "body"),]

listofDFs <- list(geneAnnotDF,promAnnotDF,enhancerDF)
listofAnnotations <- list("Gene", "Promoter", "Enhancer")

# creates a dataframe ready to plot with ggplot
enhancerGenePlotDf <- create_annotation_analysis_df(listOfAnnotations = listofAnnotations, listOfDFs = listofDFs)
enhancerGenePlotDf$association = factor(enhancerGenePlotDf$association, levels = c("Gene", "Promoter", "Enhancer"), ordered = TRUE)

outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/enhancerGenePlotDf.RData" ) 
save(enhancerGenePlotDf, file = outputFileFig3)



###
## Figure 3C

pparDF <-  methylationPairs.annotated[str_detect(methylationPairs.annotated$pathway, "ppar"),]
mapkDF <-  methylationPairs.annotated[str_detect(methylationPairs.annotated$pathway, "mapk"),]
erbbDF <-  methylationPairs.annotated[str_detect(methylationPairs.annotated$pathway, "erbb"),]
cytokineDF <-  methylationPairs.annotated[str_detect(methylationPairs.annotated$pathway, "cytokine"),]
cell_cycleDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "cell_cycle"),]
p53DF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "p53"),]
mtorDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "mtor"),]
pi3k_aktDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "pi3k_akt"),]
apoptosisDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "apoptosis"),]
wntDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "wnt"),]
tgfbDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "tgfb"),]
vegfDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "vegf"),]
focal_adhDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "focal_adh"),]
ECMDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "ECM"),]
adh_junDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "adh_jun"),]
jak_statDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "jak_stat"),]
notchDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "notch"),]
hedgehogDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "hedgehog"),]
replicationDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "replication"),]
BERDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "BER"),]
NERDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "NER"),]
HRDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "HR"),]
NHEJDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "NHEJ"),]
mismatchDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "mismatch"),]
telomereDF <-  methylationPairs.annotated[ str_detect(methylationPairs.annotated$pathway, "telomere"),]

listofDFs <- list(pparDF, mapkDF , erbbDF , cytokineDF , cell_cycleDF , p53DF , mtorDF , pi3k_aktDF , apoptosisDF , wntDF , tgfbDF , vegfDF , focal_adhDF ,
                  ECMDF , adh_junDF , jak_statDF , notchDF , hedgehogDF , replicationDF , BERDF , NERDF , HRDF , NHEJDF , mismatchDF , telomereDF)
listofAnnotations <- list("ppar","mapk","erbb","cytokine","cell_cycle","p53","mtor","pi3k_akt","apoptosis","wnt","tgfb",	
                          "vegf","focal_adh","ECM","adh_jun","jak_stat","notch","hedgehog","replication","BER","NER","HR","NHEJ","mismatch","telomere")

pathPlotDf <- create_annotation_analysis_df(listOfAnnotations = listofAnnotations, listOfDFs = listofDFs)




pathPlotDfAllLook <- pathPlotDf[pathPlotDf$significanceNeg != "" | pathPlotDf$significancePos != "",]
table(pathPlotDfAllLook$association)

# most frequently altered are: cytokine (14), ECM (10), mapk (11), pi3k_akt (8)

pathPlotDfAllFreq <- pathPlotDf[pathPlotDf$association == "cytokine" | pathPlotDf$association == "ECM" | 
                                  pathPlotDf$association == "mapk" |
                                  pathPlotDf$association == "pi3k_akt",]
pathPlotDfAllFreq$association <- as.character(pathPlotDfAllFreq$association)
pathPlotDfAllFreq$association[pathPlotDfAllFreq$association == "mapk"] <-"MAPK" 
pathPlotDfAllFreq$association[pathPlotDfAllFreq$association == "pi3k_akt"] <-"Pi3K/Akt" 

outputFileFig3 <- paste0(figureDirectory, "Figure3-Genomic-Regions/pathPlotDfAllFreq.RData" ) 
save(pathPlotDfAllFreq, file = outputFileFig3)



################ 
## Average multiple metastases pairs of one patient
################ 

# calculate average run from all patients with multiple pairs: P18, P67, P06, P17, P42, P04, P08
P18 <- create_mean_from_runDataCols(methylationPairs.annotated$P18.BLun.1_p, methylationPairs.annotated$P18.BLun.2_p)
P67 <- create_mean_from_runDataCols(methylationPairs.annotated$P67.BLun.1_p, methylationPairs.annotated$P67.BLun.2_p, 
                                    methylationPairs.annotated$P67.BLun.3_p, methylationPairs.annotated$P67.BLun.4_p)
P06 <- create_mean_from_runDataCols(methylationPairs.annotated$P06.BLym.1_p, methylationPairs.annotated$P06.BLym.2_p)
P17 <- create_mean_from_runDataCols(methylationPairs.annotated$P17.BLym.1_p, methylationPairs.annotated$P17.BLym.2_p)
P42 <- create_mean_from_runDataCols(methylationPairs.annotated$P42.BLym.1_p, methylationPairs.annotated$P42.BLym.2_p)
P04 <- create_mean_from_runDataCols(methylationPairs.annotated$P04.BSki.1_p, methylationPairs.annotated$P04.BSki.2_p)
P08 <- create_mean_from_runDataCols(methylationPairs.annotated$P08.BSof.1_p, methylationPairs.annotated$P08.BSof.2_p, methylationPairs.annotated$P08.BSof.3_p)

# create a runDataframe containing one HMM prediction for each of the 14 patients
runPatientData <- data.frame("ID" = methylationPairs.annotated$ID,"P02-BLym" = methylationPairs.annotated$P02.BLym_p, "P03-BLun" = methylationPairs.annotated$P03.BLun_p, 
                             "P04-BSki" = P04, "P06-BLym" = P06, "P08-BSof" = P08, 
                             "P09-BLiv" = methylationPairs.annotated$P09.BLiv_p, "P16-BLun" = methylationPairs.annotated$P16.BLun_p, 
                             "P17-BLym" = P17, "P18-BLun" = P18, "P28-BSki" = methylationPairs.annotated$P28.BSki_p, 
                             "P39-BLun" = methylationPairs.annotated$P39.BLun_p, "P42-BLym" = P42, 
                             "P64-BLun" = methylationPairs.annotated$P64.BLun_p, "P67-BLun" = P67)


methylationPatientCount <- data.frame( "ID" = runPatientData$ID, "countNeg" = apply(runPatientData, MARGIN = 1, function(x) sum(x == "-")), "countPos" = apply(runPatientData, MARGIN = 1, function(x) sum(x == "+")) )
methylationPairs.annotated_count <- merge(methylationPatientCount, methylationPairs.annotated, by.x = "ID",  by.y = "ID" )


################ 
## FIGURE 4 - Meta-Information
################ 

# all meta information available for the patients (see manuscript Table 1)
metaInfoDF <- data.frame("Patient" = c("P02-Blym","P03-Blun","P04-Bski","P06-Blym","P08-Bsof","P09-Bliv","P16-Blun","P17-Blym",
                                       "P18-Blun","P28-Bski","P39-Blun","P42-Blym","P64-Blun","P67-Blun"), 
                         "TherapyTissue" = c("None", "None", "None", "Intracranial only",  "Both", "None",  "None", "None", "None", 
                                             "Both","Both", "Intracranial only",  "None", "Unknown"),
                         "therapyOption" = c("No \n treatment", "No \n treatment", "No \n treatment", "Chemo-& \n Immuno- \n therapy", 
                                             "Chemo-& \n Immuno- \n therapy", "No \n treatment",  "No \n treatment", "No \n treatment", "No \n treatment", 
                                             "Chemo-& \n Immuno- \n therapy","Immuno- \n therapy", "Chemo- \n therapy",  "No \n treatment", "Unknown"),
                         "timeBetween" = c( -1,3,5,74,0,-2,4,8,3,-1,6,21,7,0),
                         "survivalBM" = c( 10,11,0,130,4,6,4,1,14,6,61,18,4,0), 
                         "numPosCpG" = c(sum(methylationPairs.annotated$P02.BLym_p == "+"),sum(methylationPairs.annotated$P03.BLun_p == "+"),
                                         (sum(methylationPairs.annotated$P04.BSki.1_p == "+") + sum(methylationPairs.annotated$P04.BSki.1_p == "+") ) / 2 ,
                                         (sum(methylationPairs.annotated$P06.BLym.1_p == "+") + sum(methylationPairs.annotated$P06.BLym.2_p == "+") ) / 2 ,
                                         (sum(methylationPairs.annotated$P08.BSof.1_p == "+") + sum(methylationPairs.annotated$P08.BSof.2_p == "+") + sum(methylationPairs.annotated$P08.BSof.3_p == "+") ) / 3 ,
                                         sum(methylationPairs.annotated$P09.BLiv_p == "+"),sum(methylationPairs.annotated$P16.BLun_p == "+"),
                                         (sum(methylationPairs.annotated$P17.BLym.1_p == "+") + sum(methylationPairs.annotated$P17.BLym.2_p == "+") ) / 2 ,
                                         (sum(methylationPairs.annotated$P18.BLun.1_p == "+") + sum(methylationPairs.annotated$P18.BLun.2_p == "+") ) / 2 ,
                                         sum(methylationPairs.annotated$P28.BSki_p == "+"),sum(methylationPairs.annotated$P39.BLun_p == "+"),
                                         (sum(methylationPairs.annotated$P42.BLym.1_p == "+") + sum(methylationPairs.annotated$P42.BLym.2_p == "+") ) / 2 ,
                                         sum(methylationPairs.annotated$P64.BLun_p == "+"),
                                         (sum(methylationPairs.annotated$P67.BLun.1_p == "+") + sum(methylationPairs.annotated$P67.BLun.2_p == "+") + 
                                            sum(methylationPairs.annotated$P67.BLun.3_p == "+") + sum(methylationPairs.annotated$P67.BLun.4_p == "+")) / 4 ),
                         "numNegCpG" = c(sum(methylationPairs.annotated$P02.BLym_p == "-"),sum(methylationPairs.annotated$P03.BLun_p == "-"),
                                         (sum(methylationPairs.annotated$P04.BSki.1_p == "-") + sum(methylationPairs.annotated$P04.BSki.1_p == "-") ) / 2 ,
                                         (sum(methylationPairs.annotated$P06.BLym.1_p == "-") + sum(methylationPairs.annotated$P06.BLym.2_p == "-") ) / 2 ,
                                         (sum(methylationPairs.annotated$P08.BSof.1_p == "-") + sum(methylationPairs.annotated$P08.BSof.2_p == "-") + sum(methylationPairs.annotated$P08.BSof.3_p == "-") ) / 3 ,
                                         sum(methylationPairs.annotated$P09.BLiv_p == "-"),sum(methylationPairs.annotated$P16.BLun_p == "-"),
                                         (sum(methylationPairs.annotated$P17.BLym.1_p == "-") + sum(methylationPairs.annotated$P17.BLym.2_p == "-") ) / 2 ,
                                         (sum(methylationPairs.annotated$P18.BLun.1_p == "-") + sum(methylationPairs.annotated$P18.BLun.2_p == "-") ) / 2 ,
                                         sum(methylationPairs.annotated$P28.BSki_p == "-"),sum(methylationPairs.annotated$P39.BLun_p == "-"),
                                         (sum(methylationPairs.annotated$P42.BLym.1_p == "-") + sum(methylationPairs.annotated$P42.BLym.2_p == "-") ) / 2 ,
                                         sum(methylationPairs.annotated$P64.BLun_p == "-"),
                                         (sum(methylationPairs.annotated$P67.BLun.1_p == "-") + sum(methylationPairs.annotated$P67.BLun.2_p == "-") + 
                                          sum(methylationPairs.annotated$P67.BLun.3_p == "-") + sum(methylationPairs.annotated$P67.BLun.4_p == "-")) / 4 ))

metaInfoDF$numDifCpG <- metaInfoDF$numPosCpG + metaInfoDF$numNegCpG
metaInfoDF$negScore <-  (( metaInfoDF$numPosCpG -  metaInfoDF$numNegCpG ) / metaInfoDF$numDifCpG)


outputFileFig4 <- paste0(figureDirectory, "Figure4-Meta-Information/metaInfoDF.RData" ) 
save(metaInfoDF, file = outputFileFig4)


################ 
## FIGURE 5 - genome-wide gene candidates 
################ 


## Figure 5A 
#

genesPatientDF <- methylationPairs.annotated_count[,c("ID","Gene","countNeg","countPos","pathways","Chr","RegFeature")]
numCand <- get_number_of_candidates_on_cutoff(genesPatientDF)

numCand <- numCand[numCand$Patients != 0,]
outputFileFig5 <- paste0(figureDirectory, "Figure5-genome-wide-candidates/numCand.RData" ) 
save(numCand, file = outputFileFig5)



## Figure 5B 
#

cancerGeneAnnotation <- read.table( file= paste0(annotationDirectory, "annotationTable_new.txt"), header = T)
                                    
allGenes <-  unique(str_split(paste(methylationPairs.annotated_count$Gene, collapse = ";"), ";")[[1]])

# create a datatable containing annotation to each gene 
allGeneAnnotation <- data.frame("Gene" = c("a"), "pathway" = c("a"), "cancerPath" = c("a"), stringsAsFactors = F)
for(i in 1:length(allGenes)){
  gene <- allGenes[i]
  paths <- pathways[pathways$Gene == gene,8:32]
  if(length(paths[,1]) > 0 ){
    paths <- colnames(paths[which(paths == 1)])
    paths <- paste(paths, collapse = ";")
  } else {
    paths <- ""
  }
  
  cancers <- cancerGeneAnnotation[cancerGeneAnnotation$Gene == gene,8:16]
  
  if(length(cancers[,1]) > 0 ){
    cancers <- colnames(cancers[which(cancers == 1)])
    cancers <- paste(cancers, collapse = ";")
  } else{
    cancers <- ""
  }
  
  
  allGeneAnnotation <- rbind(allGeneAnnotation, c("Gene" = gene, "pathway" = paths, "cancerPath" = cancers))
}

allGeneAnnotation <- allGeneAnnotation[-1,]

allGeneAnnotation <- allGeneAnnotation[allGeneAnnotation$Gene != "",] 

# create candidate gene set with > 50% of all patients affected
candGenesDF <- methylationPairs.annotated_count[methylationPairs.annotated_count$countNeg >= 8|methylationPairs.annotated_count$countPos >= 8, ]

candGenes <- data.frame("Gene" = unique(str_split(paste(candGenesDF$Gene, collapse = ";"), ";")[[1]])) 

candGenes <- data.frame("Gene" = candGenes[candGenes$Gene != "",])

candGeneAnnotation <- merge(x=candGenes,y=allGeneAnnotation,by= 'Gene' ,all.x=TRUE)
candGeneAnnotation <- candGeneAnnotation[complete.cases(candGeneAnnotation),]

# create pathways annotation
compositionPathDF <- NULL

# get all pathways
pathAnnots <- colnames(pathways[8:32])
numCandGenes <- length(candGeneAnnotation$pathway)
numAllGenes <- length(allGeneAnnotation$pathway)

# calculate enrichment using Fisher's exact test 
for(p in 1:length(pathAnnots)){
  path <- pathAnnots[p]
  perc <- round((length(candGeneAnnotation$pathway[str_detect(candGeneAnnotation$pathway, path)]) / (length(allGeneAnnotation$pathway[str_detect(allGeneAnnotation$pathway, path)])) * 100),2)
  
  upperleft <- length(candGeneAnnotation$pathway[str_detect(candGeneAnnotation$pathway, path)])
  lowerleft <- numCandGenes - upperleft
  upperright <- length(allGeneAnnotation$pathway[str_detect(allGeneAnnotation$pathway, path)]) - upperleft
  lowerright <- numAllGenes - length(allGeneAnnotation$pathway[str_detect(allGeneAnnotation$pathway, ".")]) - lowerleft
  
  qval <- fisher.test(matrix(c(upperleft,lowerleft,upperright,lowerright), nrow = 2), alternative = "greater")$p.value
  
  compositionPathDF <- rbind(compositionPathDF, c(path, perc, qval))
}

compositionPathDF <- as.data.frame(compositionPathDF)
colnames(compositionPathDF) <- c("pathway", "percentage" , "qvalFisher" ) 
compositionPathDF$qvalFisher <- p.adjust(as.numeric(as.character(compositionPathDF$qvalFisher)), method = "fdr")



outputFileFig5 <- paste0(figureDirectory, "Figure5-genome-wide-candidates/compositionPathDF.RData" ) 
save(compositionPathDF, file = outputFileFig5)





## Figure 5C
#

# get all cancer gene annotations
cancAnnots <- colnames(cancerGeneAnnotation[8:16])

compositionCancDF <- NULL

for(c in 1:length(cancAnnots)){
  canc <- cancAnnots[c]
  perc <- round((length(candGeneAnnotation$cancerPath[str_detect(candGeneAnnotation$cancerPath, canc)]) / (length(allGeneAnnotation$cancerPath[str_detect(allGeneAnnotation$cancerPath, canc)])) * 100),2)
  
  
  upperleft <- length(candGeneAnnotation$cancerPath[str_detect(candGeneAnnotation$cancerPath, canc)])
  lowerleft <- numCandGenes - upperleft
  upperright <- length(allGeneAnnotation$cancerPath[str_detect(allGeneAnnotation$cancerPath, canc)]) - upperleft
  lowerright <- numAllGenes - length(allGeneAnnotation$cancerPath[str_detect(allGeneAnnotation$cancerPath, ".")]) - lowerleft
  
  qval <- fisher.test(matrix(c(upperleft,lowerleft,upperright,lowerright), nrow = 2), alternative = "greater")$p.value
  
  
  compositionCancDF <- rbind(compositionCancDF, c(canc, perc, qval))
}
compositionCancDF <- as.data.frame(compositionCancDF)
colnames(compositionCancDF) <- c("cancerPath", "percentage" , "qvalFisher" ) 

compositionCancDF$percentage <- as.numeric(as.character(compositionCancDF$percentage))
compositionCancDF$qvalFisher <- as.numeric(as.character(compositionCancDF$qvalFisher))
compositionCancDF$qvalFisher <- p.adjust(compositionCancDF$qvalFisher, method = "fdr")


outputFileFig5 <- paste0(figureDirectory, "Figure5-genome-wide-candidates/compositionCancDF.RData" ) 
save(compositionCancDF, file = outputFileFig5)



## Figure 5D
#



candGene.df <- bitr(candGenes$Gene, fromType = "SYMBOL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)

allGene.df <- bitr(allGenes, fromType = "SYMBOL",
                   toType = c("ENTREZID"),
                   OrgDb = org.Hs.eg.db)


ego <- enrichGO(gene          = candGene.df$ENTREZID,
                universe      = names(allGene.df$ENTREZID),
                OrgDb         = org.Hs.eg.db,
                ont           = "all",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                readable      = TRUE)



enrichdf <- data.frame("GOType" = c("a"), "Description" = c("a"), "qval" = c(0), "perc" = c(0),stringsAsFactors = F)
egoATable <- ego@result
egoATable <- egoATable[order(egoATable$qvalue),]
for (i in 1:20){
  go <- egoATable$ONTOLOGY[i]
  descr <- egoATable$Description[i]
  qval <- egoATable$qvalue[i]
  perc <- round((as.numeric(str_split(as.character(egoATable$GeneRatio[i]), "/")[[1]][1]) / 
                   as.numeric(str_split(as.character(egoATable$GeneRatio[i]), "/")[[1]][2]) ) * 100,2)
  enrichdf <- rbind(enrichdf, c(go, descr, qval, perc))
  
}
enrichdf$qval <- as.numeric(enrichdf$qval)
enrichdf <- enrichdf[-1,]
enrichdf2 <- enrichdf
enrichdf2$qval <- - log10(enrichdf2$qval)


outputFileFig5 <- paste0(figureDirectory, "Figure5-genome-wide-candidates/enrichdf2.RData" ) 
save(enrichdf2, file = outputFileFig5)


################ 
## FIGURE 6 - FreqGenes heatmap
################ 


mostFreqCpGs <- candGenesDF[str_detect(candGenesDF$pathways, "mapk") | 
                               str_detect(candGenesDF$pathways, "pi3k_akt") | 
                               str_detect(candGenesDF$pathways, "ECM") | 
                               str_detect(candGenesDF$pathways, "cytokine") ,c(1,33:56)]
rownames(mostFreqCpGs) <- mostFreqCpGs[,1]
mostFreqCpGs <- mostFreqCpGs[,-1]


mostFreqCpGsAnnot <- candGenesDF[candGenesDF$ID %in% rownames(mostFreqCpGs),c("ID","GenomicRegion","pathways")]

annotationGeneName <- NULL
for(i in 1:length(rownames(mostFreqCpGs))){
  genes <- ""
  if(rownames(mostFreqCpGs)[i] %in% candGenesDF$ID){
    genes <- as.character(candGenesDF$Gene[candGenesDF$ID == rownames(mostFreqCpGs)[i]])
  }
  if(genes != ""){
    gene <- paste0(unique(str_split(genes, pattern = ";")[[1]]), collapse = ";")
    gene <- paste0(gene, "_" , rownames(mostFreqCpGs)[i])
    annotationGeneName <- c(annotationGeneName, gene)
  } else {
    annotationGeneName <- c(annotationGeneName,  rownames(mostFreqCpGs)[i])
  }
  
}

rownames(mostFreqCpGs) <- annotationGeneName

rownames(mostFreqCpGsAnnot) <- annotationGeneName
mostFreqCpGsAnnot <-  mostFreqCpGsAnnot[,-1, drop = F]



annotationRowsN <- mostFreqCpGsAnnot
annotationRowsN$GenomicRegion[annotationRowsN$GenomicRegion == "enhancer;promoter"] <- "promoter"
annotationRowsN$cytokine <- as.character(str_detect(annotationRowsN$pathways, "cytokine"))
annotationRowsN$Pi3kAkt <- as.character(str_detect(annotationRowsN$pathways, "pi3k_akt"))
annotationRowsN$ECM <- as.character(str_detect(annotationRowsN$pathways, "ECM"))
annotationRowsN$MAPK <- as.character(str_detect(annotationRowsN$pathways, "mapk"))
annotationRowsN <- annotationRowsN[,-2]

colnames(annotationRowsN) <- c("Region","cytokine","Pi3kAkt","ECM","MAPK")

mostFreqAnnotationRowList <- list(mostFreqCpGs, annotationRowsN)

outputFileFig6 <- paste0(figureDirectory, "Figure6-FreqGenes Heatmap/mostFreqAnnotationRowList.RData" ) 
save(mostFreqAnnotationRowList, file = outputFileFig6)



################ 
## FIGURE 7 - Joint analysis of gene body methylation and gene expression
################ 

literatureCpGs <- mostFreqCpGs

literatureCpGs$Gene <- unlist(lapply(rownames(literatureCpGs), FUN = function(x){ str_split(x,"_")[[1]][1] } ) )
literatureCpGs$Region <- annotationRowsN$Region

# read in expression data
expressionData <- read.csv(file=paste0(dataDirectory,
                                       "Table_availableExpressionData.tsv"),
                           sep = "\t", check.names = F) 


candidateGenes <- expressionData$`Associated Gene Name`
candidateGenesDecr <- candidateGenes[candidateGenes != "MECOM" & candidateGenes != "BMP4"]


# mean expresion daata for multiple samples of the same metastases
expressionDataMeaned <- data.frame("Ensembl.Gene.ID" = expressionData$`Ensembl Gene ID`,"Associated.Gene.Name" = expressionData$`Associated Gene Name`,
                                   "Description" = expressionData$Description,"Chr" = expressionData$Chr,"Start" = expressionData$Start,"End" = expressionData$End,
                                   "Strand" = expressionData$Strand,"Length" = expressionData$Length,"P18_B" = expressionData$P18_B,
                                   "P18_Lun" = (expressionData$`P18_Lun-1` + expressionData$`P18_Lun-2`) / 2,"P16_B" = expressionData$P16_B,
                                   "P16_Lun" = expressionData$P16_Lun,"P3_B" = expressionData$P03_B,"P3_Lu" = expressionData$P03_Lu,
                                   "P39_B" = expressionData$P39_B,"P39_Lu" = expressionData$P39_Lu,"P4_B" = expressionData$P04_B,
                                   "P4_Ski" = (expressionData$`P04_Ski-1` + expressionData$`P04_Ski-2` ) /2,"P42_B" = expressionData$P42_B ,
                                   "P42_Lym" = (expressionData$`P42_Lym-1` + expressionData$`P42_Lym-2`) / 2,"P6_B" = (expressionData$`P06_B-1` + expressionData$`P06_B-2` ) / 2,
                                   "P6_Lym" = expressionData$P06_Lym,"P8_B" = expressionData$P08_B,"P8_Sof" = (expressionData$`P08_Sof-1` + expressionData$`P08_Sof-2` + expressionData$`P08_Sof-3` ) /3)

# create paired expression data (just like methylation pairs)
expressionDataPaired <- data.frame("Ensembl.Gene.ID" = expressionData$`Ensembl Gene ID`,"Associated.Gene.Name" = expressionData$`Associated Gene Name`,
                                   "Description" = expressionData$Description,"Chr" = expressionData$Chr,"Start" = expressionData$Start,"End" = expressionData$End,
                                   "Strand" = expressionData$Strand,"Length" = expressionData$Length,
                                   "P18_BLun_1" = expressionData$P18_B - expressionData$`P18_Lun-1`, "P18_BLun_2" = expressionData$P18_B - expressionData$`P18_Lun-2`,
                                   "P16_BLun" = expressionData$P16_B - expressionData$P16_Lun,"P03_BLun" = expressionData$P03_B - expressionData$P03_Lu,
                                   "P39_BLun" = expressionData$P39_B - expressionData$P39_Lu,
                                   "P04_BSki_1" = expressionData$P04_B - expressionData$`P04_Ski-1` , "P04_BSki_2" = expressionData$P04_B - expressionData$`P04_Ski-2` , 
                                   "P42_BLym_1" = expressionData$P42_B - expressionData$`P42_Lym-1`,"P42_BLym_2" = expressionData$P42_B - expressionData$`P42_Lym-2`, 
                                   "P06_BLym_1" = expressionData$`P06_B-1` - expressionData$P06_Lym, "P06_BLym_2" = expressionData$`P06_B-2` - expressionData$P06_Lym,
                                   "P08_BSof_1" = expressionData$P08_B - expressionData$`P08_Sof-1`, "P08_BSof_2" = expressionData$P08_B - expressionData$`P08_Sof-2`, 
                                   "P08_BSof_3" = expressionData$P08_B - expressionData$`P08_Sof-3` ) 

allSharedPatientsM <- c("P03.BLun","P39.BLun","P42.BLym.1","P42.BLym.2","P04.BSki.1","P04.BSki.2","P06.BLym.1","P06.BLym.2","P08.BSof.1","P08.BSof.2","P08.BSof.3","P16.BLun",
                        "P18.BLun.1","P18.BLun.2")
allSharedPatientsE <- c("P03_BLun","P39_BLun","P42_BLym_1","P42_BLym_2","P04_BSki_1","P04_BSki_2","P06_BLym_1", "P06_BLym_2","P08_BSof_1","P08_BSof_2","P08_BSof_3","P16_BLun",
                        "P18_BLun_1","P18_BLun_2")

# create data frame that connects gene expression and body methylation
exprMethylConnection <- NULL
for(g in 1:length(expressionDataMeaned$Associated.Gene.Name)){
  gene <- as.character(expressionDataMeaned$Associated.Gene.Name[g])
  for(p in 1:length(allSharedPatientsE)){
    expressionVal <- expressionDataPaired[expressionDataPaired$Associated.Gene.Name == gene, allSharedPatientsE[p]]
    
    methylationPeakPlus <- max(literatureCpGs[str_detect(literatureCpGs$Gene,gene) & ! str_detect(literatureCpGs$Region,"^promoter"), allSharedPatientsM[p]] )
    methylationPeakMinus <- min(literatureCpGs[str_detect(literatureCpGs$Gene,gene) & ! str_detect(literatureCpGs$Region,"^promoter"), allSharedPatientsM[p]] )
    
    methylationPeak <- ifelse(methylationPeakPlus < abs(methylationPeakMinus), methylationPeakPlus, methylationPeakMinus)
    
    
    if(methylationPeak != -Inf){
      exprMethylConnection <- rbind(exprMethylConnection, c("Gene" = gene, "Patient" = allSharedPatientsE[p], "value" = expressionVal, 
                                                            "data" = "expression"))
      exprMethylConnection <- rbind(exprMethylConnection, c("Gene" = gene, "Patient" = allSharedPatientsE[p], "value" = methylationPeak, 
                                                            "data" = "methylation"))
    }
  }
  
}

exprMethylConnection <- as.data.frame(exprMethylConnection)
colnames(exprMethylConnection) <- c("Gene","Patient","value","data")

# group genes according to their differences in median
medianDiff <- NULL
for(g in 1:length(unique(exprMethylConnection$Gene))){
  gene <- unique(as.character(exprMethylConnection$Gene))[g]
  
  medianDiff <- rbind(medianDiff, c("Gene" = gene, "medianDiff" =  
                                      median(as.numeric(as.character(exprMethylConnection$value[exprMethylConnection$Gene == gene & exprMethylConnection$data == "methylation"])))-
                                      median(as.numeric(as.character(exprMethylConnection$value[exprMethylConnection$Gene == gene & exprMethylConnection$data == "expression"]))) ) ) 
}

exprMethylConnection <- merge(exprMethylConnection, medianDiff, by = "Gene")
exprMethylConnection$medianDiff <- as.numeric(as.character(exprMethylConnection$medianDiff))

exprMethylConnection$value <- as.numeric(as.character(exprMethylConnection$value))

exprMethylConnection$data <- factor(exprMethylConnection$data, levels = c("methylation", "expression"))


outputFileFig7 <- paste0(figureDirectory, "Figure7-expression/exprMethylConnection.RData" ) 
save(exprMethylConnection, file = outputFileFig7)



# test all genes for an expression of metastases pairs unequal to zero
ExprTest <- c()
for(e in 1:length(unique(exprMethylConnection$Gene))){
  gene <- as.character(unique(exprMethylConnection$Gene)[e])
  diffZeroPval <- t.test(exprMethylConnection$value[exprMethylConnection$Gene == gene & 
                                                      exprMethylConnection$data == "expression"],
                         mu = 0)$p.value
  
  ExprTest <- rbind(ExprTest,c(gene, mean(exprMethylConnection$value[exprMethylConnection$Gene == gene & 
                                                                       exprMethylConnection$data == "expression"]), diffZeroPval ))
}
ExprTest <- as.data.frame(ExprTest)
colnames(ExprTest) <- c("gene", "mean expression","pvalue")
ExprTest$pvalue <- as.numeric(as.character(ExprTest$pvalue))
ExprTest$adj.pvalue <- p.adjust(ExprTest$pvalue, method ="fdr")


# expression intra vs extra for top ranking genes TNXB JAK3 MECOM

topGenes <- c("TNXB","JAK3","MECOM")

# create dataframe to direct intracranial and extracranial expression comparison of the three candidate genes
ExpressionComparisonPlot <- NULL
for(g in 1:length(topGenes)){
  gene <- topGenes[g]
  
  geneExpressionDataMeaned <- expressionDataMeaned[expressionDataMeaned$Associated.Gene.Name == gene,9:24]
  
  for(n in 1:length(colnames(geneExpressionDataMeaned))){
    ExpressionComparisonPlot <- rbind(ExpressionComparisonPlot, c("gene" = gene,"value" = geneExpressionDataMeaned[[n]], 
                                                                  "patient" = colnames(geneExpressionDataMeaned)[n],
                                                                  "region" =  ifelse(str_detect(colnames(geneExpressionDataMeaned)[n],"_B"), "intracranial", "extracranial")))
  }
}

ExpressionComparisonPlot <- as.data.frame(ExpressionComparisonPlot)

ExpressionComparisonPlot$value <- as.numeric(as.character(ExpressionComparisonPlot$value))
ExpressionComparisonPlot$region <- factor(ExpressionComparisonPlot$region, levels = c("intracranial", "extracranial"))



outputFileFig7 <- paste0(figureDirectory, "Figure7-expression/ExpressionComparisonPlot.RData" ) 
save(ExpressionComparisonPlot, file = outputFileFig7)

