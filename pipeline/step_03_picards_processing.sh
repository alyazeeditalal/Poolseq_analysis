#!/bin/bash
#----------------------------------------
#This script is intended to process alignment bams using picard tools
#----------------------------------------
#location of config file
source config.sh

#----------------------------------------

for bam in $map_file/*.bam; do
    individual=$(basename $(echo $bam) .bam)

if [ -s "$picard_bam/${individual}-sorted-md-rg.bam" ]; then 
     echo "bam file for $individual has been 
     1-sorted 
     2-duplicated reads tagged 
     3- read groups were changed
     4- bams indexed"

else 
    flowcell=$(samtools view $bam |
           head -n 1 |
           sed 's:\t.*::' |
           sed 's/.*:.*:\(.*\):\(.*\):.*:.*:.*/\1.\2/')


    echo "bam files for $individual have not been sorted...
          ....Sorting using Samtools...."
      
    samtools sort $bam\
      --thread 8\
      -o $picard_bam/${individual}-sorted.bam 2> $log_align_process/${individual}.sorted.log
      #-T $temp/${individual}.temp 
 
    echo "Tagging duplicates for $individual bam...
            ....Tagging duplicated using picard...."

    gatk --java-options "-Xmx8G" MarkDuplicates \
    -I $picard_bam/${individual}-sorted.bam \
    -O $picard_bam/${individual}-sorted-md.bam \
    --REMOVE_DUPLICATES true \
    --QUIET true \
    -M $output_dup/${individual}-md-metrics.txt 2> $log_align_process/${individual}.sorted-md.log

 
    echo "changing read groups for $individual bam...
            ....Changing read group using picard...."

    gatk --java-options "-Xmx8G" AddOrReplaceReadGroups \
    -I $picard_bam/${individual}-sorted-md.bam \
    -O $picard_bam/${individual}-sorted-md-rg.bam \
    -RGID ${individual} \
    -RGLB ${individual} \
    -RGPL illumina \
    -RGPU "$flowcell".$individual \
    -RGSM ${individual} \
    --QUIET true 2> $log_align_process/${individual}.sorted-md-rg.log

    samtools index $picard_bam/${individual}-sorted-md-rg.bam


    rm $picard_bam/${individual}-sorted.bam #Remove unwanted files 
    rm $picard_bam/${individual}-sorted-md.bam  #remove unwanted files 

fi

done

echo 'Bams processing done!'
#FIN-----------------------------------------------------