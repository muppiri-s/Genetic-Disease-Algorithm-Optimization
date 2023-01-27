# Genetic-Disease-Algorithm-Optimization

The paper uses two input text data files: the first dataset which consists of the average eye sizes for each strain from the DGRP dataset. The second is the gene expression table. The current implementation of the sequential algorithm to identify the genes which is correlated with the eye size accurately gives the correlation. 
Overall code: 20.73 minutes

Break down:
Algorithm1: 16.4 mins
Algorithm2: 7.21 secs
Algorithm3: 16.28 mins 
Overall code: 20.73 mins (Combined time of algorithm 1 2, and 3)


After using parallel library in R. 
Code type	Execution Time
Parallel (32 cores)	3.47 minutes  (208.688) seconds
Parallel (64 cored)	1.657 minutes (99.436 seconds)

Reduced the code execution time to more than half. 
