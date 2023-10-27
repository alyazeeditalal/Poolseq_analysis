# Poolseq pipeline

This is a bash pipeline to analyse poolseq data. The analysis produces genome wide population genetics data, with a main focus on statistical analyses such as Tajima's D and Fst. 

The pipeline can also be used to produce genome wide bulk segregate analysis results (BSA) from an F2 segregating population. 

The input data could be poolseq genomic data from parents for population genetics analysis and a poolseq data from segregating populaiton for bulk segregate analysis. 

The population genetic analysis will be conducted for all samples and pairwise Fst will be computed between all samples. The population genetic data will be produced at 
1- the SNP level 
2- a determined window size 
3- Gene-wise    
