#!/bin/bash
#-----------------------------------------------------
#This script creates required directory for the pipeline
#-----------------------------------------------------

#create pipeline directories 
mkdir -p ../0_QC/fasta_qc #fastqc
mkdir -p ../0_QC/bam_dup #duplication logs 
mkdir -p ../0_QC/vcf_stat #vcf stat folder
mkdir -p ../0_QC/bams_QC #bam qc
mkdir ../2_mapping #directory for mapped files   
mkdir ../3_processing #directory for processed files
mkdir ../4_varscan_variants #directory for varscan vcf
mkdir ../5_freebayes_variants #directory for freebayes vcf 
mkdir ../6_bcftools_variants #directory for bcftools vcf 
mkdir ../7_population_genetics #directry for the population genetics
mkdir -p ../8_BSA_analysis/results/snps #directory for BSA analysis
mkdir -p ../8_BSA_analysis/plots/snps #directory for BSA analysis

#create a directory for log files
mkdir -p ../log/fasta_qc #fastqc logs 
mkdir -p ../log/bwa_align_log #alignemnts logs
mkdir -p ../log/align_process #alignemnts files processing
mkdir -p ../log/bamqc  #bam qc
mkdir -p ../log/varscan  #varscan_vcf
mkdir -p ../log/freebayes #freebayes
mkdir -p ../log/bcftools #bcftools
mkdir -p ../log/snpeff #snpeff logs 
mkdir -p ../log/vcf_table #vcf to table 
mkdir -p ../log/grenedalf  #grendalf logs


