# Poolseq pipeline

This is a bash pipeline to analyse poolseq data. The analysis produces genome wide population genetics data including Tajima's D, Fst and allele frequency difference that is useful to identify QTLs from pooled segregating individuals.

The pipeline can also be used to produce genome wide bulk segregate analysis results (BSA) from an F2 segregating population. 

The input data could be poolseq genomic data for population genetics analysis or a poolseq data from segregating populaiton for bulk segregate analysis. 

The population genetic analysis will be conducted for all samples and pairwise Fst will be computed between all samples. 

1. Fst at the SNP level 
2. a determined window size 
3. Gene-wise    

A snp table will be produced from freebayes variants that can be used for BSA analysis. 

_______________________________________________________
The pipeline steps can be run together or seperately by running each step of the script on its own. 

To run this pipeline consider creating a conda environemnt from the yaml file 

```
conda env create --file poolseq.yml -n env_name
````

Start by saving the fastq files in the 1_data/ directory 

to run the entire pipeline use the run_all.bash

```
./run_all.bash
```

__________________________________________________


The pipeline steps include 