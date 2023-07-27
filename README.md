# Identifying Genes Related to Retinitis Pigmentosa in Drosophila Melanogaster using Gene Expression Data

## Introduction
This repository contains the R code implementation for identifying genes related to Retinitis Pigmentosa (RP) in Drosophila Melanogaster using gene expression data. The algorithm aims to discover genes that are correlated with eye size, which could potentially be associated with RP in the fruit fly.

## Supplementary Material
The repository includes the following supplementary material:
- R-code for the implemented algorithm
- DGRP Expression Data: The gene expression data file, which is in a tabular format with replicates as columns and genes as rows.
- Eye Size Data: The eye size data file, containing the average eye sizes for each strain in the Drosophila Melanogaster DGRP dataset.

## Environment and Libraries
The code is built using RStudio Version 4.1.2. The following R libraries are necessary for implementing parallel programming:
- `parallel`
- `foreach`
- `doMC`

## Execution Time
The original sequential execution time of the algorithm breakdown is as follows:
- Algorithm1: 16.4 minutes
- Algorithm2: 7.21 seconds
- Algorithm3: 16.28 minutes
Overall code execution time: 20.73 minutes

After implementing parallel processing using the `doMC` library, the execution time significantly improved:
- Parallel (32 cores): 3.47 minutes (208.688 seconds)
- Parallel (64 cores): 1.657 minutes (99.436 seconds)
The parallel implementation reduced the code execution time to more than half, enhancing the efficiency of the analysis.

## Input Data
- **Gene Expression Data**: The gene expression file is expected to be in a tabular format, where the first row represents the header with replicate names and subsequent rows contain expression values for each gene.
- **Eye Size Data**: The eye size file is expected to be in a tabular format with the first column representing the strain names and the second column containing the average eye size for each strain.

## Output
The algorithm identifies genes correlated with eye size and outputs the results to the file specified as `out.txt`. The output file contains a list of genes sorted by their correlation values in descending order.

## Usage
1. Ensure RStudio Version 4.1.2 is installed along with the required libraries.
2. Place the gene expression data file (`dgrp_expression_Female.txt`) and the eye size data file (`Rh1G69D.txt`) in the same directory as the R code.
3. Execute the R code in RStudio or any R environment.

Please note that the code uses parallel processing to speed up execution time, so it is recommended to have a sufficient number of CPU cores available for optimal performance.
