#!/bin/bash

vcf_input='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.vcf'
out_dir='../7_population_genetics/'
log_dir='../log/grenedalf/'

# Create directories
mkdir -p $out_dir
mkdir -p $log_dir

# Generate sync file
grenedalf sync --threads 8 --vcf-path $vcf_input --out-dir $out_dir --file-prefix population --log-file $log_dir/population.sync.log 

# Run CMH test
cmh-test.pl --input ../7_population_genetics/population.sync --output ../7_population_genetics/population.cmh --min-count 12 --min-coverage 50 --max-coverage 200 --population 1-2