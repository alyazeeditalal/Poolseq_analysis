#!/bin/bash


#In this script I will merge the fasta fies of each pair if there are more than one file per pair 
########################################################################

#output directory
output_dir='../1_data' #Directory were the fasta files will be saved 
#input directory
input_dir='../../01.RawData' #Directory contaiing subdirectories with fasta files 

#merge fasta files that are of the same pair 
for dir in $input_dir/*; do 
individual=$(basename $(echo $dir))

 for fa in $individual; do 
  cat $input_dir/${individual}/*_1.fq.gz >  $output_dir/${individual}_1.fq.gz
  cat $input_dir/${individual}/*_2.fq.gz >  $output_dir/${individual}_2.fq.gz
 done
done 