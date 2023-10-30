#!/bin/bash
#----------------------------------------
#This script is intended to call varaints from the bam files and filter them
#----------------------------------------

#location of config file
source config.sh

#----------------------------------------
#indexing the genome to produce fai indexed genome 
if [ -s "${loc_genome}.fai" ]; then
   echo "${loc_genome}.fai has been generated" 
else 
   echo "Generating .fai index" 
    samtools faidx $loc_genome
fi 

#----------------------------------------
#run freebayes

#create a bam list to run samtools mpileup 
echo "creating bams list for samtools mpileup" 
ls ../3_processing/*sorted-md-rg.bam > ../input/pool.bam.txt

#Varaint calling analysis using freebayes requires a very long time especially
#with the large number of ployidy in the pooled samples 

if [ -s $vcf_freebayes ]; then 
   echo $vcf_freebayes "file has been already generated"
else 
echo "generating" $vcf_freebayes

TMPDIR=$tmp freebayes-parallel \
  <(fasta_generate_regions.py ${loc_genome}.fai 100000) 8 \
   --fasta-reference ${loc_genome} \
   --use-best-n-alleles 4 \
   --pooled-continuous \
   --bam-list ../input/pool.bam.txt > $vcf_freebayes 2> $log_freebayes/freebayes.log

fi


echo "freebayes variant calling is Done!" 
#----------------------------------------