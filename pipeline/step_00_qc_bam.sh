#!/bin/bash

#----------------------------------------
#This script is intending to produce QC report for all the bam file       
#----------------------------------------

#location of picard bams 
picard_bam='../3_processing'

#location of output bams 
bams_QC='../0_QC/bams_QC'

log_bamqc='../log/bamqc'
#----------------------------------------

for bam in $picard_bam/*-sorted-md-rg.bam; do
    individual=$(basename $(echo $bam) -sorted-md-rg.bam)

if [ -s "../0_QC/bams_QC/${individual}.flagstat" ] && [ -d "../0_QC/bams_QC/${individual}" ]; then 
     echo "Quality check for $individual bam is completed"

else       
    
 samtools flagstat -@8 $bam > $bams_QC/$individual.flagstat
     
 qualimap bamqc -nt 8 --java-mem-size=4G -bam $bam -c -outdir $bams_QC/${individual} 2> $log_bamqc/${individual}.log

fi
done