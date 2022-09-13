## ---------------------------
##
## This script was created with R version 4.2.1 (2022-06-23)
##
##
## Title: Patient-specific identification of genome-wide DNA-methylation differences between intracranial and extracranial melanoma metastases
##
## Purpose of script: Helper functions for script MainScript.R
##
## Author: Theresa Kraft
##
## Date Created: 08-15-2022
## Date last modified: 09-13-2022
##
## R-code and data usage license: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
##
## Copyright (c) Theresa Kraft, 2022
## Email: theresa.kraft@tu-dresden.de
##
## ---------------------------
##
## Notes: This is a collections of all helper functions used in the MainScript.R
##   
##
## ---------------------------



#############################################################################################################
##
### function calculates autocorrelation of any dataframe in a specific order
### measurementDF: Dataframe containing measurement values in specific order
##
############################################################################################################

create_methyl_acf_df <- function(measurementDF){
  rnaDataM = as.matrix(measurementDF[,7:30])
  
  chrIdx = sapply(unique(measurementDF$Chr), function(chr) which(measurementDF$Chr == chr))
  
  neededChr = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
  allMethAcfs = lapply (neededChr, function(chr) {
    sapply(1:ncol(rnaDataM), function(sampleIdx) {
      acfRes = acf(rnaDataM[chrIdx[[chr]],sampleIdx], lag.max = 100,plot = F)
      acfRes$acf[,1,1]
    })
  })
  names(allMethAcfs) = neededChr
  
  
  neededQuantiles = c(0 ,0.1 ,0.25 ,0.5 ,0.75, 0.9, 1 )
  
  allMethAcfsQuantiles = NULL
  for(chr in neededChr) {
    q = t(sapply(1:nrow(allMethAcfs[[chr]]),function(lag) quantile(allMethAcfs[[chr]][lag,],neededQuantiles)))
    colnames(q) = gsub("\\%","",paste0("quant_",colnames(q)))
    allMethAcfsQuantiles = rbind(allMethAcfsQuantiles, data.frame(chr=chr,lag = 1:nrow(allMethAcfs[[chr]]), q))
  }
  
  
  allMethAcfQuantilesChrCondensed = NULL
  for(lag in sort(unique(allMethAcfsQuantiles$lag))) {
    allMethAcfQuantilesChrCondensed =
      rbind(allMethAcfQuantilesChrCondensed,
            data.frame(lag=lag,t(colMeans(allMethAcfsQuantiles[ allMethAcfsQuantiles$lag==lag,-c(1,2)]))))
  }
  
  
  
  return(allMethAcfQuantilesChrCondensed)
}






#############################################################################################################
##
### function that takes x = column of data and extracts all pathways from annotated genes
### x: column of dataframe to get all paths from
##
############################################################################################################

extract_paths <- function(x){
  genes <- unique(strsplit(as.character(x[["Gene"]]), split = ";")[[1]])
  path <- c()
  if(length(genes) > 0){
    for(g in 1:length(genes)){
      gene <- genes[g]
      if(length(pathwayArray[pathwayArray$V1 == gene,2]) != 0){
        path <- c(path, strsplit(as.character(pathwayArray[pathwayArray[,1] == gene,2]), split = ";")[[1]])
      }
    }
  }
  path <- sort(unique(path))
  path <- paste(path, collapse = ";")
  return(path)
} 




#############################################################################################################
##
### function that creates a dataframe to make a barplot with 
###   listOfAnnotations: list of how the annotations should be called (e.g. "Gene", "Promoter", "Others")
###   listOfDFs: list of dataframes which contain CpGs that are in corresponding annotations than listOfAnnotations 
###       (e.g. geneDF - containing gene CpGs only, promoterDF - containing promoter CpGs only, otherDF - 
##        containing CpGs in no gene no promoter)    
##
############################################################################################################
create_annotation_analysis_df <- function(listOfAnnotations, listOfDFs){
  
  if(length(listOfAnnotations) != length(listOfDFs)){
    stop("Length of lists are not equal! Returning without calculations.")
  }
  ##
  # calculates the percentages of CpGs in + and - state dependent on their annotationDF
  ##
  
  # reset all values
  PatientID <- c()
  association <- c() 
  valuePositive <- c() 
  valueNegative <- c() 
  
  # calculate percentage
  for(i in 1:length(sampleNameDF)){
    sampleValues <- as.character(sampleNameDF[i])
    samplePredictions <- as.character(colnames(sampleNameDF)[i])
    
    for(l in 1:length(listOfAnnotations)){
      PatientID <- c(PatientID, sampleValues )
      
      association <- c(association, listOfAnnotations[[l]])
      
      singleDF <- listOfDFs[[l]]
      highMethyl<- nrow(singleDF[singleDF[,samplePredictions] == "+",])
      lowMethyl <- nrow(singleDF[singleDF[,samplePredictions] == "-",])
      valuePositive <- c(valuePositive, highMethyl / nrow(singleDF) ) 
      valueNegative <- c(valueNegative, lowMethyl/ nrow(singleDF)) 
    }
    
  }
  
  # form barplot dataframe in right format
  barplotDFPos <- data.frame(PatientID, association, value = valuePositive)
  barplotDFNeg <- data.frame(PatientID, association, value =  -valueNegative )
  barplotDF <- rbind(barplotDFPos,barplotDFNeg)
  
  ##
  # rank the sample names according to their tissue type 
  ##
  ranking <- c()
  for(i in 1:length(barplotDF[,1])){
    if(str_detect(barplotDF[i,1], "Lun")){
      ranking <- c(ranking,5)
    }
    if(str_detect(barplotDF[i,1], "Lym")){
      ranking <- c(ranking,4)
    }
    if(str_detect(barplotDF[i,1], "Ski")){
      ranking <- c(ranking,3)
    }
    if(str_detect(barplotDF[i,1], "Liv")){
      ranking <- c(ranking,1)
    }
    if(str_detect(barplotDF[i,1], "Sof")){
      ranking <- c(ranking,2)
    }
  }
  barplotDF <- cbind(barplotDF,ranking)
  
  ##
  # add significance according to fisher test to each row 
  ##
  
  fisherTestDF <- as.data.frame(create_fisher_test_dataframe(listOfAnnotations, listOfDFs))
  fisherTestDF$valuePV <- p.adjust(as.numeric(as.character(fisherTestDF$valuePV)), method = "fdr")
  
  significancePos <- c() 
  significanceNeg <- c() 
  pval <- c()
  or <- c()
  for(i in 1:length(barplotDF[,1])){
    association <- as.character(barplotDF[i,2])
    
    if(barplotDF[i,3] > 0){
      value <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                      fisherTestDF[,4] == association & 
                                                      fisherTestDF[,5] == "Increased",][1,2]))
      pv <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                   fisherTestDF[,4] == association & 
                                                   fisherTestDF[,5] == "Increased",][1,3]))
      
      # significance check
      if(value > 2 && pv < (0.05 / (24 * length(listOfAnnotations)))){
        significancePos <- c(significancePos,"X")
      } else {
        significancePos <- c(significancePos, "")
      }
      significanceNeg <- c(significanceNeg,"")
      
    } else {
      value <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                      fisherTestDF[,4] == association & 
                                                      fisherTestDF[,5] == "Decreased",][1,2]))
      pv <- as.numeric(as.character(fisherTestDF[fisherTestDF[,1] == as.character(barplotDF[i,1]) & 
                                                   fisherTestDF[,4] == association & 
                                                   fisherTestDF[,5] == "Decreased",][1,3]))
      
      # significance check
      if(pv < 0.05 ){
        significanceNeg <- c(significanceNeg,"X")
      } else {
        significanceNeg <- c(significanceNeg, "")
      }
      significancePos <- c(significancePos,"")
      
    }
    
    
    pval <- c(pval, pv )
    or <- c(or, value)
    
  }
  
  barplotDF <- cbind(barplotDF,significancePos,significanceNeg, pval, or)
  
  
  return(barplotDF)
  
}



############################################################################################################
##
### function for creation of a fisherTest dataframe containing the Odds Ratio of each sample and positive/negative methylation
###   listOfAnnotations: list of how the annotations should be called (e.g. "Gene", "Promoter", "Others")
###   listOfDFs: list of dataframes which contain CpGs that are in corresponding annotations than listOfAnnotations 
###       (e.g. geneDF - containing gene CpGs only, promoterDF - containing promoter CpGs only, otherDF - containing CpGs in no gene no promoter)    
##
############################################################################################################
create_fisher_test_dataframe <- function(listOfAnnotations, listOfDFs){
  fisherTestAll <- c(samples = as.character(sampleNameDF[1,]))
  
  for(l in 1:length(listOfDFs)){
    assocDdf <- listOfDFs[[l]]
    
    fisherTestAll <- cbind(fisherTestAll, create_fisher_all_samples(assocDdf, T,F), 
                           create_fisher_all_samples(assocDdf, T,T), 
                           create_fisher_all_samples(assocDdf, F,F), 
                           create_fisher_all_samples(assocDdf, F,T))
  }
  
  fishertestOR <- NULL
  
  for(i in 1:length(fisherTestAll[,1])){
    a <- 1
    for(j in seq(2,length(fisherTestAll[1,]),by = 4)){
      fishertestOR <- rbind(fishertestOR, cbind(samples = as.character(fisherTestAll[i,1]), 
                                                valueOR = fisherTestAll[i,j], valuePV = fisherTestAll[i,j + 1],
                                                Association = listOfAnnotations[[a]] , Methylation = "Increased"))
      fishertestOR <- rbind(fishertestOR, cbind(samples = as.character(fisherTestAll[i,1]), 
                                                valueOR = fisherTestAll[i,j + 2], valuePV = fisherTestAll[i,j+3],
                                                Association = listOfAnnotations[[a]] , Methylation = "Decreased"))
      a <- a + 1
    }
  }
  
  return(fishertestOR)
}

############################################################################################################
##
### function to perform a fisher test for all samples (needs function fisher_test_by_sample) 
## associationDF: dataframe containing CpGs associated with annotation category
## positive: boolean value to check if we are looking at increased methylation (otherwise: decreased methylation)
## pvalue: boolean value to check if we are looking at pvalue (otherwise: estimate)
##
############################################################################################################
create_fisher_all_samples <- function(associationDF, positive, pvalue){
  
  if(pvalue){
    return (c(fisher_test_by_sample(2,associationDF,positive)$p.value,fisher_test_by_sample(3,associationDF,positive)$p.value,
              fisher_test_by_sample(4,associationDF,positive)$p.value,fisher_test_by_sample(5,associationDF,positive)$p.value,
              fisher_test_by_sample(6,associationDF,positive)$p.value,fisher_test_by_sample(7,associationDF,positive)$p.value,
              fisher_test_by_sample(8,associationDF,positive)$p.value,fisher_test_by_sample(9,associationDF,positive)$p.value,
              fisher_test_by_sample(10,associationDF,positive)$p.value,fisher_test_by_sample(11,associationDF,positive)$p.value,
              fisher_test_by_sample(12,associationDF,positive)$p.value,fisher_test_by_sample(13,associationDF,positive)$p.value,                                              
              fisher_test_by_sample(14,associationDF,positive)$p.value,fisher_test_by_sample(15,associationDF,positive)$p.value,
              fisher_test_by_sample(16,associationDF,positive)$p.value,fisher_test_by_sample(17,associationDF,positive)$p.value,
              fisher_test_by_sample(18,associationDF,positive)$p.value,fisher_test_by_sample(19,associationDF,positive)$p.value,
              fisher_test_by_sample(20,associationDF,positive)$p.value,fisher_test_by_sample(21,associationDF,positive)$p.value,
              fisher_test_by_sample(22,associationDF,positive)$p.value,fisher_test_by_sample(23,associationDF,positive)$p.value, 
              fisher_test_by_sample(24,associationDF,positive)$p.value,fisher_test_by_sample(25,associationDF,positive)$p.value))
  }
  if(!pvalue){
    return (c(fisher_test_by_sample(2,associationDF,positive)$estimate,fisher_test_by_sample(3,associationDF,positive)$estimate,
              fisher_test_by_sample(4,associationDF,positive)$estimate,fisher_test_by_sample(5,associationDF,positive)$estimate,
              fisher_test_by_sample(6,associationDF,positive)$estimate,fisher_test_by_sample(7,associationDF,positive)$estimate,
              fisher_test_by_sample(8,associationDF,positive)$estimate,fisher_test_by_sample(9,associationDF,positive)$estimate,
              fisher_test_by_sample(10,associationDF,positive)$estimate,fisher_test_by_sample(11,associationDF,positive)$estimate,
              fisher_test_by_sample(12,associationDF,positive)$estimate,fisher_test_by_sample(13,associationDF,positive)$estimate,                                              
              fisher_test_by_sample(14,associationDF,positive)$estimate,fisher_test_by_sample(15,associationDF,positive)$estimate,
              fisher_test_by_sample(16,associationDF,positive)$estimate,fisher_test_by_sample(17,associationDF,positive)$estimate,
              fisher_test_by_sample(18,associationDF,positive)$estimate,fisher_test_by_sample(19,associationDF,positive)$estimate,
              fisher_test_by_sample(20,associationDF,positive)$estimate,fisher_test_by_sample(21,associationDF,positive)$estimate,
              fisher_test_by_sample(22,associationDF,positive)$estimate,fisher_test_by_sample(23,associationDF,positive)$estimate, 
              fisher_test_by_sample(24,associationDF,positive)$estimate,fisher_test_by_sample(25,associationDF,positive)$estimate))
  }
  
}


############################################################################################################
##
### function to perform a fisher test on enrichment of association in sample 
### sampleNr: index number of column with state annotation in annotationDF and associationDF 
### associationDF: dataframe which only contains the CpGs associated with the annotation category of interest
### positive: boolean for positive or negative methylation analysis 
##
############################################################################################################
fisher_test_by_sample <- function(sampleNr, associationDF, positive){
  
  # number of all, all that are positive and all that are negative in this sample 
  numAll <- sum(methylationPairs.annotated[,sampleNr] == "+") + sum(methylationPairs.annotated[,sampleNr] == "-") + sum(methylationPairs.annotated[,sampleNr] == "=")  
  numAllPos <- sum(methylationPairs.annotated[,sampleNr]== "+")  
  numAllNeg <- sum(methylationPairs.annotated[,sampleNr]== "-")
  
  # number of all associated and all positive and negative associated CpGs
  numAssoc <- sum(associationDF[,sampleNr] == "+") +
    sum(associationDF[,sampleNr] == "-") + 
    sum(associationDF[,sampleNr] == "=")
  numAssocPos <- sum(associationDF[,sampleNr] == "+")
  numAssocNeg <- sum(associationDF[,sampleNr] == "-")
  
  ##
  # Negative methylation
  ##
  if(! positive){
    # associated with region; in - state
    obenLinks <- numAssocNeg
    
    # not associated with region: in - state
    obenRechts <- numAllNeg - numAssocNeg
    
    # associated with region; not in - state
    untenLinks <- numAssoc - numAssocNeg
    
    # not associated with region; not in - state
    untenRechts <- (numAll - numAllNeg) - (numAssoc - numAssocNeg)
  }
  
  ##
  # Positive methylation
  ##
  if(positive){
    # associated with region; in + state
    obenLinks <- numAssocPos
    
    # not associated with region: in - state
    obenRechts <- numAllPos - numAssocPos
    
    # associated with region; not in - state
    untenLinks <- numAssoc - numAssocPos
    
    # not associated with region; not in - state
    untenRechts <- (numAll - numAllPos) - (numAssoc - numAssocPos)
  }
  
  fisherTestMatrix <- matrix(c(obenLinks,untenLinks,obenRechts, untenRechts), ncol = 2 , nrow = 2)
  return(fisher.test(fisherTestMatrix, alternative = "greater"))
  
  
}


############################################################################################################
##
### function that counts the differences in methylation predictions between two metastases pairs
### pat1: Name of metastases pair 1
### pat2: Name of metastases pair 2
##
############################################################################################################

countDiffsBetweenPairs <- function(pat1, pat2){
  col1 <- methylationPairs.annotated[,names(sampleNameDF[which(sampleNameDF == pat1)])]
  col2 <- methylationPairs.annotated[,names(sampleNameDF[which(sampleNameDF == pat2)])]
  comparison <- paste0(pat1, " vs ", pat2)
  
  return(c(comparison, 
           sum((col1 == "+" & col2 == "+")),
           sum((col1 == "-" & col2 == "-")),
           sum((col1 == "=" & col2 == "=")),
           sum((col1 == "+" & col2 == "=")),
           sum((col1 == "=" & col2 == "+")),
           sum((col1 == "-" & col2 == "=")),
           sum((col1 == "=" & col2 == "-")),
           sum((col1 == "+" & col2 == "-")),
           sum((col1 == "-" & col2 == "+"))))
}


############################################################################################################
##
### function that creates a mean run from multiple metastases pairs
### firstRow: first metastases pair of one patient
### secondRow: second metastases pair of one patient
### thirdRow (optional): third metastases pair of one patient
### fourthRow (optional): fourth metastases pair of one patient
##
############################################################################################################



create_mean_from_runDataCols <- function(firstRow, secondRow, thirdRow = NULL, fourthRow = NULL){
  patientRun <- c()
  firstRow <- convert_runData_to_numeric(firstRow)
  secondRow <- convert_runData_to_numeric(secondRow)
  allRows <- data.frame(firstRow,secondRow)
  
  if(!is.null(thirdRow)){
    thirdRow <- convert_runData_to_numeric(thirdRow)
    allRows <- data.frame(firstRow,secondRow, thirdRow)
  }
  if(!is.null(fourthRow)){
    fourthRow <- convert_runData_to_numeric(fourthRow)
    allRows <- data.frame(firstRow,secondRow, thirdRow, fourthRow)
  } 
  
  
  patientRNr <- pbapply::pbapply(allRows, FUN = mean, MARGIN = 1)
  
  patientR <- patientRNr
  patientR[patientRNr > 0 ] <- "+"
  patientR[patientRNr < 0 ] <- "-"
  patientR[patientRNr == 0 ] <- "="
  
  return(patientR)
}

convert_runData_to_numeric <- function(data){
  data <- as.character(data)
  data[data == "+"] <- 1
  data[data == "-"] <- -1
  data[data == "="] <- 0
  data <- as.numeric(data)
  return(data)
}



############################################################################################################
##
### function calculates number of candidate genes dependent on specific patient cutoff
### genesPatientDF: dataframe that includes countNeg and countPos 
##
############################################################################################################


get_number_of_candidates_on_cutoff <- function(genesPatientDF){
  numCand <- data.frame("Patients" = c(0), "Number" = c(0), "Direction" = c("0"), "Level" = c("0"), stringsAsFactors = F)
  for (i in 0:14){
    numPos <- genesPatientDF$Gene[genesPatientDF$countPos >= i]
    numNeg <- genesPatientDF$Gene[genesPatientDF$countNeg >= i]
    genesP <- unique(strsplit(paste(as.character(numPos), collapse = ";"), split = ";")[[1]])
    genesN <- unique(strsplit(paste(as.character(numNeg), collapse = ";"), split = ";")[[1]])
    genesP <- genesP[genesP != ""]
    genesN <- genesN[genesN != ""]
    
    numCand <- rbind(numCand, c("Patients" = i, "Number" = length(genesP), "Direction" = "Increased", "Level" = "Gene"))
    numCand <- rbind(numCand, c("Patients" = i, "Number" = length(genesN), "Direction" = "Decreased", "Level" = "Gene"))
  }  
  numCand <- numCand[-1,]
  numCand$Number <- as.numeric(as.character(numCand$Number))
  numCand$Patients <- as.numeric(as.character(numCand$Patients))
  
  return(numCand)
}

