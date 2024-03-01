#!/bin/bash

#----------------------------------------------------------------
# This script is used to carry out cmh test between pooled samples uing popoolation
#----------------------------------------------------------------
#location of config files and input files

#In vcf file 
vcf_input='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.vcf'

#out directory
out_dir='../7_population_genetics/'

#log file
log_grenedalf='../log/grenedalf/'

#location of config file
source config.sh

#----------------------------------------------------------------
#Generating a sync file 
echo "${YELLOW}
..........Generating a sync file............" 

grenedalf sync \
--threads $t \
--vcf-path $vcf_input \
--out-dir $out_dir \
--file-prefix population \
--log-file $log_grenedalf/population.sync.log 

#----------------------------------------------------------------
#cmh popoolation 
echo "${RED}
.........running CMH test ........." 

perl ~/programs/popoolation2_1201/cmh-test.pl --input ../7_population_genetics/populationsync.sync \
 --output ../7_population_genetics/17F2_population.cmh --min-count 12 --min-coverage 50 --max-coverage 200 --population $population
