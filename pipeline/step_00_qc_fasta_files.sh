#!/bin/bash

#----------------------------------------------------------------
#This script is intending to produce QC report for all the fast file       
#----------------------------------------------------------------

#----------------------------------------------------------------
#Run fast QC if files do not exist
#----------------------------------------------------------------

#Check if *fastqc* files already exist in the output directory
if ls ../0_QC/fasta_qc/*fastqc* 1> /dev/null 2>&1; then
#if exist echo a message
    echo "FastQC files already exist. Skipping the command"
else
    echo "......fastqc files do not exist..... 
    ........Generating fastqc files........"
    parallel "fastqc {} -o ../0_QC/fasta_qc" ::: ../1_data/*
fi

#FN----------------------------------------------------------------