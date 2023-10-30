#!/bin/bash

#----------------------------------------
#This script aligns reads using BWA to the reference genome and immediately coverted to a bam file
#----------------------------------------
#location of config file
source config.sh

#----------------------------------------
for file in $input_reads/*_R1.fastq.gz; do 
individual=$(basename $(echo $file) _R1.fastq.gz) #assign individaul name 
out_file="../2_mapping/${individual}.bam" #assign outfile 
log_file="../log/bwa_align_log/${individual}.log" #Assign a log file 

  if [ -f "$map_file/${individual}.bam" ]; then 
      echo "bam fies for $individual already exist skipping the command"
  else 
      echo "bam files for $individual do not exist...
            ....running BWA alignment"

     bwa-mem2 mem -M -t 16 \
	    $loc_genome $input_reads/${individual}_R1.fastq.gz \
	    $input_reads/${individual}_R2.fastq.gz 2> $log_file \
      | samtools view -buS - > $out_file
   
  fi 

done 

#FIN-----------------------------------------------------