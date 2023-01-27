#Manuscript Title: Identifying genes related to Retinitis Pigmentosa in Drosophila Melanogaster using gene expression data
#Supplementary Material
#R-code for the implemented algorithm
#code built using RStudio Version 4.1.2

#necessary Libraries to implement parallel programming
library(parallel)
library(foreach)
library(doMC)

#Register the required number of cores
registerDoMC(64)


inputFileGeneExpression = "dgrp_expression_Female.txt";
inputFileEyeSizes = "Rh1G69D.txt";
outputFileGenes = "out.txt";

# Gene expression file's expected format is as in original source. 
# Eye size file's expected format is as in the following example:
# Header1 Header2<newline>
# Line1Name Eye size<newline>
# and so on.
# First line is the header line.
# Content of all other lines follows the original source.
# Column separator is either the space, the comma or the tab character.
possibleColumnSeparators="[ ,\t]"

# Read the expression data file, first.
print("Reading expression data...")
expressionRawCon<-file(inputFileGeneExpression)
expressionRaw <- readLines(expressionRawCon)
close(expressionRawCon)

# split the lines into fields by the tab character
print("Preprocessing expression data...")
expressionFieldsLines <- Reduce(c,
                                lapply(expressionRaw, 
                                       function(x) strsplit(x,
                                                            possibleColumnSeparators)))

# Read the phenotype file
print("Reading eye size data...")
con1<-file(inputFileEyeSizes)
eyeSizeRaw <- readLines(con1)
close(con1)

# split each line into fields by the tab character
print("Preprocessing eye size data...")
eyeSizesAndLines <- Reduce(c,
                           lapply(eyeSizeRaw[2:length(eyeSizeRaw)],
                                  function(x) strsplit(x,possibleColumnSeparators)))
strainNames <- Reduce(c,lapply( eyeSizesAndLines, function(x) x[1]))

# get only the eye sizes and put them in ''eyeSizes''
eyeSizes <- Reduce(c, lapply(eyeSizesAndLines, function(x) as.numeric(x[2])))
maxEyeSizes <- max(eyeSizes)
minEyeSizes <- min(eyeSizes)

#select the replicates whether they match the selected extreme lines above
splitFunction<-function(x){
  numStr <- strsplit(strsplit(x,":")[[1]][1], "_")[[1]][2]; 
  sprintf("%03d",as.numeric(numStr))
}

getExpressionFieldsIndexVector <- function(indexVector) unlist(lapply(expressionFieldsLines[[1]], function(x) paste("RAL",splitFunction(x), sep = "") %in% strainNames[indexVector]))
print("Done reading and preprocessing input data.")

#settings
print(sprintf("Min eye size: %f", minEyeSizes))
print(sprintf("Max eye size: %f", maxEyeSizes))

quantilePercents = c(0.209, 0.872)
quantileNumbers = quantile(c(minEyeSizes,maxEyeSizes), quantilePercents)

extremeTopThreshold = quantileNumbers[[2]]

extremeBottomThreshold = quantileNumbers[[1]]

print(sprintf("Top threshold of eye size: %f", extremeTopThreshold))
print(sprintf("Bottom threshold of eye size: %f", extremeBottomThreshold))

# Assumption: all the lines in the phenotype file are sorted according to alphebetical order of the line names(a-z, 0-1000). Same assumption is for the column headers(replicate names) in the expression file. (21 is in front of 210)
# get the index vectors for the selected lines that satisfy the extreme condition

smallEyeSizesIndexVector <- sapply(eyeSizes, function(x) x < extremeBottomThreshold)
largeEyeSizesIndexVector <- sapply(eyeSizes, function(x) x > extremeTopThreshold)

eyeSizesBooleanIndexVector <- (smallEyeSizesIndexVector| largeEyeSizesIndexVector)
eyeSizeIndex <- c(which(smallEyeSizesIndexVector),which(largeEyeSizesIndexVector))

numberOfSelectedLines <- length(which(eyeSizesBooleanIndexVector))

strainList <- strainNames[which(eyeSizesBooleanIndexVector)] 

print(sprintf("Number of selected lines: %d", numberOfSelectedLines))
smallEyeSizes <- eyeSizes[((1:length(eyeSizes))[smallEyeSizesIndexVector])]
largeEyeSizes <- eyeSizes[((1:length(eyeSizes))[largeEyeSizesIndexVector])]

largeEyeSizesLineName <- strainNames[((1:length(eyeSizes))[largeEyeSizesIndexVector])]
largeEyeSizesRaw <- eyeSizeRaw[((1:length(eyeSizes))[largeEyeSizesIndexVector])+1]

smallEyeSizesLineName <- strainNames[((1:length(eyeSizes))[smallEyeSizesIndexVector])]
smallEyeSizesRaw <- eyeSizeRaw[((1:length(eyeSizes))[smallEyeSizesIndexVector])+1]

expressionFieldsBottomIndexVector<-getExpressionFieldsIndexVector(smallEyeSizesIndexVector)
expressionFieldsTopIndexVector<-getExpressionFieldsIndexVector(largeEyeSizesIndexVector)

rangeOfBottom <- 1:(length(which(smallEyeSizesIndexVector)))
rangeOfTop <- (1+length(which(largeEyeSizesIndexVector))):(numberOfSelectedLines)

expressionFieldsIndexVector <- c(which(expressionFieldsBottomIndexVector), 
                                 which(expressionFieldsTopIndexVector))

# get all the replicates name just to double check. Also,get their column numbers.
globalSelectedReplicateNames <- expressionFieldsLines[[1]][expressionFieldsIndexVector]
acceptedExpressionLinePosition <- expressionFieldsIndexVector

# get all the rows with only the selected columns(replicates) as well as the header
expressionFieldsExtremesWithHeader <- sapply(expressionFieldsLines, 
                                             function(x) x[expressionFieldsIndexVector])
expressionMatrixReplRowGeneCol <<- matrix(
  as.numeric(
    expressionFieldsExtremesWithHeader[,2:length(expressionFieldsLines)]), 
  nrow = numberOfSelectedLines * 2, 
  byrow = FALSE)


# get all the gene names
geneNamesTemp <- sapply(expressionFieldsLines, function(x) x[[1]])
geneNames <- geneNamesTemp[2:length(geneNamesTemp)]
numberOfGeneNames <- length(geneNames)

algorithm1 <- function (X, Y, C, Ess, exprLevelsMatrix ) {
  quantilePercents = c(X,Y)
  quantileNumbers = quantile(c(minEyeSizes,maxEyeSizes), quantilePercents)
  
  extremeTopThreshold = quantileNumbers[[2]]
  extremeBottomThreshold = quantileNumbers[[1]]
  
  smallEyeSizesIndexVector <- sapply(Ess, function(x) x < extremeBottomThreshold)
  largeEyeSizesIndexVector <- sapply(Ess, function(x) x > extremeTopThreshold)
  
  eyeSizesBooleanIndexVector <- (smallEyeSizesIndexVector| largeEyeSizesIndexVector)
  
  eyeSizeIndex <- c(which(smallEyeSizesIndexVector),which(largeEyeSizesIndexVector))
  selSizes <- Ess[eyeSizeIndex]
  
  strainList <- strainNames[which(eyeSizesBooleanIndexVector)] 
  
  replCombs <- algorithm2 (strainList)
  result <- algorithm3 ( selSizes, exprLevelsMatrix, replCombs, C )
  
  print("Replicate combination with maximum score:")
  print(result[[1]])
  unSortedGenes <- result[[2]]
  names(unSortedGenes) <- geneNames
  sortedGenes <- unSortedGenes[order(-abs(unSortedGenes))]
  genesDataFrame <- as.data.frame(sortedGenes)
  write.table(genesDataFrame, 
              row.names=TRUE, 
              col.names=FALSE,
              file=outputFileGenes, 
              quote=FALSE,
              sep=",")
}

algorithm2 <- function(strainArr) {
  variableN <- length(strainArr)
  replicateFlags <- sapply(seq(0,2^variableN - 1,1),
                           function(x){ as.integer(intToBits(x)[1:variableN])})
  
  replicateCombinations <-
    sapply(1:2^variableN, 
           function (y) sapply(seq(1,variableN,1),
                               function (x) {
                                 replicateFlags[,y][x]+2 * (x-1)+1}))
}

algorithm3<-function(selSizes, exprLevelsMatrix, replCombs, threshold){
  maxCount <- 0
  f<-function(x){
    temp <- as.vector(cor(selSizes,
                          exprLevelsMatrix[replCombs[,x],][,1:numberOfGeneNames]))
    temp1 <- (temp < (-threshold)); 
    temp2 <- (temp > (threshold)); 
    temp <- temp1 | temp2
    tableRes <- table(temp);
    tempTableTrueCount <- 0 
    if (length(tableRes) > 1) {
      tempTableTrueCount <- tableRes[[2]]
    }
    if (tempTableTrueCount > maxCount) {
      maxCount <<- tempTableTrueCount
      bestCombination <<- x;
      correlationResultsOfTheBestCombination<<-temp
    }
    temp<-NA
    
    c(maxCount,bestCombination)
    
  }
  
  #Custom Defined Reduction Function
  mymin<-function(a,b)
  {
    if(a[1][1]<b[1][1]) b else a
  }
  
  #parallel section using dynamic scheduling
  res<-foreach(x=1:2^numberOfSelectedLines,.combine = mymin) %dopar%
    {
      f(x)
    }
  
  
  CorrVals <- cor(selSizes,
                  exprLevelsMatrix[replCombs[,res[2]],][,1:numberOfGeneNames])
  list(globalSelectedReplicateNames[replCombs[,res[2]]], 
       CorrVals)
}

silent_output = algorithm1(quantilePercents[1], quantilePercents[2], 0.6, eyeSizes, expressionMatrixReplRowGeneCol )
silent_output = NA


