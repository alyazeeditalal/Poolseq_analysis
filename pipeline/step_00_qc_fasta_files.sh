#!/bin/bash

#----------------------------------------------------------------
#This script is intending to produce QC report for all the fast file       
#----------------------------------------------------------------
#output directory for the fastQC files 
#out_fastqc='../0_QC/fasta_qc'

#----------------------------------------------------------------
#Run fast QC if files do not exist
#----------------------------------------------------------------

#Check if *fastqc* files already exist in the output directory
if [ -n "$(find ../0_QC/fasta_qc -maxdepth 1 -type f -name '*fastqc*' -print -quit)" ]; then
#if exist echo a message
    echo "FastQC files already exist. Skipping the command"
else
    echo "......fastqc files do not exist..... 
    ........Generating fastqc files........"
    parallel "fastqc {} -o ../0_QC/fasta_qc" ::: ../1_data/*.fastq.gz
fi

#FN----------------------------------------------------------------