#!/bin/bash
#-------------------------------------------------------
#This script is intended to download the genome and index it using BWA
#------------------------------------------------------

#Creat a genome directory and download a genome manualy

#location of config file
source config.sh
#-------------------------------------------------------
#check if BWA index files exist

if [ -f "${loc_genome}" ] && \
 [ -f "${loc_genome}.0123" ] && \
 [ -f "${loc_genome}.amb" ] && \
 [ -f "${loc_genome}.ann" ] && \
 [ -f "${loc_genome}.bwt.2bit.64" ] && \
 [ -f "${loc_genome}.pac" ]; then
    echo "BWA index files already exist. Skipping indexing."
else 
    echo "Indexing genome using BWA..."
    bwa-mem2 index $loc_genome #indexing the genome using BWA
fi

#FIN-----------------------------------------------------