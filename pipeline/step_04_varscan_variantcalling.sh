#!/bin/bash
#----------------------------------------
#This script is intended to call varaints from the bam files using varscan
#----------------------------------------
#location of config file
source config.sh

#----------------------------------------
# I will start by creating an mpileup file from the bam files 
#----------------------------------------

#index reference genome using samtools 
samtools faidx ${loc_genome}

#create a bam list to run samtools mpileup 
echo "creating bams list for samtools mpileup" 
ls ../3_processing/*sorted-md-rg.bam > pool.bam.fofn

#Run samtools mpileup 
if [ -s $mpileup ]; then 
   echo $mpileup "$file has been already generated"
else 
   echo "$.....Generating...." $mpileup
   samtools mpileup -B -f $loc_genome -b pool.bam.fofn --output $mpileup 2> "$log_varscan/mpileup.log"
fi

#----------------------------------------
#varscan variant calling  

#adjust the name of the vcf head to remove the path from each name 
ls ../3_processing/*sorted-md-rg.bam | sed 's/\..\/3_processing\///; s/-sorted-md-rg\.bam$//' > Varscan_samples_name

if [ -s $vcf_varscan ]; then 
   echo $vcf_varscan "file has been already generated"
else 
echo "generating" $vcf_varscan

$varscan mpileup2snp $mpileup \
--min-var-freq 0.01 \
--min-coverage 10 \
--vcf-sample-list Varscan_samples_name \
--output-vcf > $vcf_varscanc
fi 

#Filtering
#filtered missing genotypes
#Keeping only snps 
#removing snps with GQ<99 and DP<100
#Keep biallelic SNPs

filtered_vcf=$(basename "$vcf_varscan" .vcf)

if [ -s "../4_varscan_variants/${filtered_vcf}.filtered.vcf" ]; then 
   echo "${filtered_vcf}.filtered.vcf file has been already generated"
else 
   echo ".....Generating...." ${filtered_vcf}.filtered.vcf
   bcftools view -e '%QUAL<=20 && GT="./." && GQ<99 && DP<100' -m2 -M2 -v snps $vcf_varscan -o ../4_varscan_varinats/${filtered_vcf}.filtered.vcf
fi

echo "All Done!" 
#FI----------------------------------------
