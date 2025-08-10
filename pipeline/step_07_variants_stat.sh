#!/bin/bash
#----------------------------------------
#This script is intended to filter varaints called using varscan
#----------------------------------------

#------------------------
#Generating a stat file using bcftools
#------------------------

for vcf in ../*_variants/*.vcf
do 

echo "processing $vcf file" 

# Use basename to extract the base name without the .vcf extension
vcf_name=$(basename "$vcf" .vcf)


#Generate a stat for the vcf file 
if [ -s "../0_QC/vcf_stat/${vcf_name}.stats" ]; then 
   echo "../0_QC/vcf_stat//${vcf_name}.stats file has been already generated"
else 
   echo "Generating $vcf_name.stats" 
   bcftools stats $vcf > ../0_QC/vcf_stat/${vcf_name}.stats

fi 

done 

#FN----------------------------------------