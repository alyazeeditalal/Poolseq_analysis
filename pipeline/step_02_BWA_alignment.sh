#!/bin/bash

#----------------------------------------
#This script aligns reads using BWA to the reference genome and immediately coverted to a bam file
#----------------------------------------
#location of config file
source config.sh

#----------------------------------------
for file in ../1_data/*${fa_ex_1}; do 
individual=$(basename $(echo $file) ${fa_ex_1}) #assign individaul name 
out_file="../2_mapping/${individual}.bam" #assign outfile 
log_file="../log/bwa_align_log/${individual}.log" #Assign a log file 

  if [ -f "../2_mapping/${individual}.bam" ]; then 
      echo "bam fies for $individual already exist skipping the command"
  else 
      echo "bam files for $individual do not exist...
            ....running BWA alignment"

     bwa-mem2 mem -M -t 16 \
	    $loc_genome ../1_data/"${individual}${fa_ex_1}" \
	    ../1_data/"${individual}${fa_ex_2}" 2> $log_file \
      | samtools view -buS - > $out_file
   
  fi 

done 

#FIN-----------------------------------------------------