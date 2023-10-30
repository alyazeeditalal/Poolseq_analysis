#!/bin/bash
#----------------------------------------
#This script will use gatk tools to convert the SNPs to a table
#----------------------------------------

#location of config file
source config.sh

#----------------------------------------
#Convert the snps to table
table=$(basename "$vcf_to_table" .vcf)
if [ -s "../5_freebayes_variants/${table}.table.tsv" ]; then 
   echo ${table}.table.tsv "file has been already generated"
else 
   echo "processing $table"  
   
gatk --java-options "-Xmx8G" VariantsToTable \
-V ${vcf_to_table} \
-F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
-O ../5_freebayes_variants/${table}.table.tsv 2> $log_vcf_table/${table}.log
fi

#----------------------------------------
#Convert the indels to table
table=$(basename "$indels_to_table" .vcf)
if [ -s "../5_freebayes_variants/${table}.table.tsv" ]; then 
   echo ${table}.table.tsv "file has been already generated$"
else 
   echo "processing $table"
   
gatk --java-options "-Xmx8G" VariantsToTable \
-V ${indels_to_table} \
-F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF GQ -GF PL \
-O ../5_freebayes_variants/${table}.table.tsv 2> $log_vcf_table/${table}.log
fi

#FN----------------------------------------