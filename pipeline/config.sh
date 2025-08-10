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
#wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz

# Download the latest version of the genome annotation 
# in a newly created directory called annotation
#creat annotation directory 
#mkdir ../annotation

#Manually add annotation file 
#wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/gff3/anopheles_gambiae/Anopheles_gambiae.AgamP4.57.gff3.gz

#Extract genes only GFF file 
#awk '{if($3=="gene") print $0}' Anopheles_gambiae.AgamP4.57.gff3 > gene.only.gff

#Make a directory and call it input 
#Add input file called pools.txt 
#In this file add the name of the pools and the size of the pools in a comma seperated format 

#-----------------------------------------------------
####fastqc####

#Reads extension
#input_fastqc='../1_data/*.fq.gz'

#-----------------------------------------------------
####indexing the genome####

#location of the genome
loc_genome='../genome/Anopheles_gambiae.AgamP4.dna.toplevel.fa'

#-----------------------------------------------------
####alignment to genome####

#Reads extension  
fa_ex_1='_1.fq.gz'
fa_ex_2='_2.fq.gz'

#-----------------------------------------------------
####freebayes variant calling####

## location of directory to save temporary files 
tmp="/tmp"

#-----------------------------------------------------
####variants annotation using snpEff####

#find the name of the species database in the snpeff database 
#java -jar snpEff.jar databases | grep -i gambiae

#If the database is not added to the snpEff database then it has to be added manually following the steps outlined below. 
#To do so find the snpEff config file 

#find $PWD -type f -name snpEff.config   #This will identify the config file 
#Add this line for the funestus genome for instance 
#A.funestus.genome : Anopheles funestus genome from vector base 
#Download the fasta and gff file 
#mkdir -p data/A.funestus
#cd /path/to/snpEff/data/A.funestus
#wget gff
#wget fasta
#mv fasta sequences.fa
#mv gff genes.gff

genome_DB='Anopheles_gambiae'
log_snpeff='../log/snpeff'

#-----------------------------------------------------
#####variants filteration#####

#find $PWD -type f -name vcfEffOnePerLine.pl to find the oneperline script
#location of theone perline script   
oneperline='/Users/talal.alyazeedi/miniforge3/envs/poolseq/share/snpeff-5.2-1/scripts/vcfEffOnePerLine.pl'
#-----------------------------------------------------
####Fst, diversity and allele frequency difference####

#gene regions gff 
genes_only='../annotation/gene.only.gff'

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
#color parameters for the pipeline messages 

YELLOW=`tput setaf 3 bold`

RED=`tput setaf 1 bold`

BLUE=`tput setaf 4 bold`

GREEN=`tput setaf 10 bold`

rest=`tput init`

