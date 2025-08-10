#!/bin/bash

snps_Data='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.table.tsv'

# Create output directories
mkdir -p ../8_BSA_analysis/results/snps/
mkdir -p ../8_BSA_analysis/plots/snps/

# Install R libraries
Rscript scripts/install.libraries.R

# Run BSA analysis
Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=FUM_Alive_test --highBulk_size=50 --lowBulk=FUM_dead_test --lowBulk_size=50 --windowSize=1e6 --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/