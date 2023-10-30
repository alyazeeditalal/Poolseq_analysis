#!/bin/bash
#----------------------------------------
#In this script I will annotate the filtered variants using snpeff
#----------------------------------------

#location of config file
source config.sh

#----------------------------------------

#Building the database for a new genome if the database doesn't exist
#snpEff build -gff3 $genome_DB 2> $log_snpeff/database.log

#----------------------------------------

#Annotating SNPs
if [ -s $annot_vcf ]; then 
   echo $annot_vcf "file has been already generated"
else 
echo "generating" $annot_vcf

snpEff $genome_DB -stats ../5_freebayes_variants -v $in_vcf > $annot_vcf 2> $log_snpeff/annotation.log
fi 

