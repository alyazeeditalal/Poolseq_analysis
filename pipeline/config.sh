#!/bin/bash
#-----------------------------------------------------
#This script is a config file that contains all the parameters for the pipeline
#-----------------------------------------------------

#prepare required data 

#make a directory and call it 1_data and move all the fastqc files there  

# Download the latest version of the genome. 
# In a newly created directory called genome
#creat a genome directory 
#mkdir ../genome

#Manually add a genome to the directory
#wget https://vectorbase.org/common/downloads/Current_Release/AfunestusFUMOZ/fasta/data/VectorBase-65_AfunestusFUMOZ_Genome.fasta


# Download the latest version of the genome annotation 
# in a newly created directory called annotation
#creat annotation directory 
#mkdir ../annotation

#Manually add annotation file 
#wget https://vectorbase.org/common/downloads/Current_Release/AfunestusFUMOZ/gff/data/VectorBase-65_AfunestusFUMOZ.gff

#Extract genes only GFF file 
#awk '{if($3=="protein_coding_gene") print $0}' sorted.VectorBase-54_AfunestusFUMOZ.gff > gene.only.gff

#Make a directory and call it input 
#Add input file called pools.txt 
#In this file add the name of the pools and the size of the pools in a comma seperated format 

#-----------------------------------------------------
####fastqc####

#output directory for the fastQC files 
#out_fastqc='../0_QC/fasta_qc'

#-----------------------------------------------------
####indexing the genome####

#location of the genome
loc_genome='../genome/VectorBase-65_AfunestusFUMOZ_Genome.fasta'

#-----------------------------------------------------
####alignment to genome####

#location of the input reads directory
input_reads='../1_data/'

#location of output reads 
map_file='../2_mapping'
#-----------------------------------------------------
####picard processing####

#location directory of bam files
input_bam='../2_mapping'

#location of output bams 
picard_bam='../3_processing'

output_dup='../0_QC/bam_dup' #location of duplicate files 

log_align_process='../log/align_process' #Alignment files processing log directory  
#-----------------------------------------------------
####bam_qc####

#location of output bams 
bams_QC='../0_QC/bams_QC'

log_bamqc='../log/bamqc'

#input is the picard bam
#-----------------------------------------------------
####freebayes variant calling####

## location of output vcf file
vcf_freebayes='../5_freebayes_variants/pooled.continuous.freebayes.parallel.vcf'

## location of temporary directory 
tmp="../../../lstm_scratch"

#log files 
log_freebayes='../log/freebayes'

#-----------------------------------------------------
####variants annotation using snpEff####

#If the database is not added to the snpEff database then it has to be added manually. 
#To do so find the snpEff config file 

#find $PWD -type f -name snpEff.config   #This will identify the config file 
#Add this line for the funestus genome for instance 
#A.funestus.genome : Anopheles funestus genome from vector base 
#Download the fasta and gff file 
#mkdir -p data/A.funestus
#cd /path/to/snpEff/data/A.funestus
#wget https://vectorbase.org/common/downloads/Current_Release/AfunestusFUMOZ/gff/data/VectorBase-65_AfunestusFUMOZ.gff 
#wget https://vectorbase.org/common/downloads/Current_Release/AfunestusFUMOZ/fasta/data/VectorBase-65_AfunestusFUMOZ_Genome.fasta
#mv VectorBase-65_AfunestusFUMOZ_Genome.fasta sequences.fa
#mv VectorBase-65_AfunestusFUMOZ.gff genes.gff

genome_DB='A.funestus'
log_snpeff='../log/snpeff'

#input vcf
in_vcf='../5_freebayes_variants/pooled.continuous.freebayes.parallel.vcf'

#annotated vcf
annot_vcf='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.vcf'
#-----------------------------------------------------
####variants stats####

#assign a stat ouptut
stat_output='../0_QC/vcf_stat' 

#-----------------------------------------------------
####VCF to table####

#Determine a vcf file path depending on the desired variant calling method
vcf_to_table="../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.vcf"

indels_to_table='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.indels.vcf'

log_vcf_table='../log/vcf_table'

#-----------------------------------------------------
####Fst, diversity and allele frequency difference####

#In vcf file 
vcf_input='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.vcf'
indel_input='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.indels.vcf'

#gene regions gff 
genes_only='../annotation/gene.only.gff'

#out directory
out_dir='../7_population_genetics'

#log file
log_grenedalf='../log/grenedalf/'

#pool sizes
pool_sizes='../input/pools.txt'

#threads 
t=8

#Width of each window along the chromosome,
window_sliding_width='50000'

#Stride between windows along the chromosome
window_sliding_stride='25000'

#filter minimum coverage and minimum count
min_cov='10'
min_count='2'

#-----------------------------------------------------
#####variants filteration#####

#find $PWD -type f -name vcfEffOnePerLine.pl to find the oneperline script
#location of theone perline script   
oneperline='/home/alyazeedit/miniconda3/envs/poolseq/share/snpsift-5.1-0/scripts/vcfEffOnePerLine.pl'

#output files containing oneperline effect
snps_oneperline_file='../5_freebayes_variants/snps_oneperline.effect.txt'

indels_oneperline_file='../5_freebayes_variants/indels_oneperline.effect.txt'
####Varscan variant calling####

#location of varscan 
varscan='java -jar -Xmx100g ~/programs/varscan-master/VarScan.v2.4.4.jar'

#mpileup file 
mpileup='../4_varscan_variants/all_samples.mpileup'

#the vcf file 
vcf_varscan='../4_varscan_variants/all_samples.varscan.vcf'

#samtools log 
log_varscan='../log/varscan'

#-----------------------------------------------------
####bcftool variant calling####

## location of output vcf file
vcf_bcftools='../6_bcftools_variants/all_samples.bcftools.vcf'
log_bcftools='../log/bcftools'

#-----------------------------------------------------
#color parameters for the pipeline messages 

YELLOW=`tput setaf 3 bold`

RED=`tput setaf 1 bold`

BLUE=`tput setaf 4 bold`

GREEN=`tput setaf 10 bold`

rest=`tput init`