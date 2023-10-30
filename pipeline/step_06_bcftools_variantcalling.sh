#!/bin/bash
#----------------------------------------
#This script is intended to call varaints from the bam files 
#----------------------------------------

#location of config file
source config.sh

#----------------------------------------
#create a bam list to run samtools mpileup 
echo "creating bams list for samtools mpileup" 
ls ../3_processing/*sorted-md-rg.bam > pool.bam.txt

#adjust the name of the vcf head to remove the path from each name 
ls ../3_processing/*sorted-md-rg.bam | sed 's/\..\/3_processing\///; s/-sorted-md-rg\.bam$//' > samples_name.txt

#----------------------------------------
#bcftools variant calling  

if [ -s $vcf_bcftools ]; then 
   echo $vcf_bcftools "file has been already generated"
else 
echo "generating" $vcf_bcftools

bcftools mpileup --threads 8 -B -Ou -f "$loc_genome" --bam-list pool.bam.txt --samples-file samples_name.txt | bcftools call --threads 8 -mv -Ov -o $vcf_bcftools 2> $log_bcftools/bcftools.log
fi 

#----------------------------------------
filtered_vcf=$(basename "$vcf_bcftools" .vcf)

if [ -s "../6_bcftools_variants/${filtered_vcf}.filtered.vcf" ]; then 
   echo "${filtered_vcf}.filtered.vcf file has been already generated$"
else 
   echo ".....Generating...." ${filtered_vcf}.filtered.vcf
   bcftools view -e '%QUAL<=20 && FMT/GT="./." && INFO/DP<100' -m2 -M2 -v snps $vcf_bcftools -o ../6_bcftools_variants/${filtered_vcf}.filtered.vcf
fi

echo "All Done!" 
#FI----------------------------------------