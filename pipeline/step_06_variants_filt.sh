#!/bin/bash
#----------------------------------------
#In this script I will filter the variants for freebayes
#----------------------------------------

#location of config file
source config.sh

#output files containing oneperline effect
#Determine a vcf file path depending on the desired variant calling method
annot_vcf='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.vcf'

snps_oneperline_file='../5_freebayes_variants/snps_oneperline.effect.txt'

indels_oneperline_file='../5_freebayes_variants/indels_oneperline.effect.txt'

#----------------------------------------

filtered_vcf=$(basename "$annot_vcf" .vcf)

if [ -s "../5_freebayes_variants/${filtered_vcf}.snps.vcf" ] && [ -s "../5_freebayes_variants/${filtered_vcf}.indels.vcf" ]; then 
   echo "${filtered_vcf}.snps.vcf file has been already generated$"
else 
   echo ".....Filtering...." ${filtered_vcf}.snps.vcf

   bcftools view -e '%QUAL<=20 && FMT/GT="./." && GQ<99 && FMT/DP<100' -m2 -M2 -v snps $annot_vcf -o ../5_freebayes_variants/${filtered_vcf}.snps.vcf
   
   bcftools query -l ../5_freebayes_variants/${filtered_vcf}.snps.vcf > ../${filtered_vcf}.snps.sample.list.txt 

   echo ".....Filtering...." ${filtered_vcf}.indel.vcf
   
   bcftools view -e '%QUAL<=20 && FMT/GT="./." && GQ<99 && FMT/DP<100' -m2 -M2 -v indels $annot_vcf -o ../5_freebayes_variants/${filtered_vcf}.indels.vcf
fi

#----------------------------------------

#Generate a one per line effect 
#Running the oneperline script
if [ -s "$snps_oneperline_file" ] && [ -s "$indels_oneperline_file" ]; then 
   echo $snps_oneperline_file "file has been already generated"
   echo $indels_oneperline_file "file has been already generated"
else 
   echo "generating one per line files"

cat ../5_freebayes_variants/${filtered_vcf}.snps.vcf | $oneperline |\
SnpSift extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].IMPACT" "EFF[*].GENE" "EFF[*].CODON" "EFF[*].AA" \
> $snps_oneperline_file

cat ../5_freebayes_variants/${filtered_vcf}.indels.vcf | $oneperline |\
SnpSift extractFields - CHROM POS REF ALT "ANN[*].EFFECT" "EFF[*].IMPACT" "EFF[*].GENE" "EFF[*].CODON" "EFF[*].AA" \
> $indels_oneperline_file
fi

#FIN---------------------------------------