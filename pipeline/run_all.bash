#!/bin/bash
#-----------------------------------------------------
#This script links all the scripts in one pipelie 
#-----------------------------------------------------
#Activate conda environment
#!/bin/bash
#-----------------------------------------------------
# Activate conda environment
if [[ $CONDA_DEFAULT_ENV != "poolseq" ]]; then
    echo "Activating poolseq environment..."
    conda activate poolseq
fi

#-----------------------------------------------------
#location of config file
source config.sh

#-----------------------------------------------------

#Making pipeline directores 
make_dir='step_00_creating_pipeline_dir.sh'

#location of the fasta file
fastqc='step_00_qc_fasta_files.sh'

#location of the BWA index file 
BWA_index='step_01_BWA_index_genome.sh' 

#location of the BWA alignemnt script
BWA_align='step_02_BWA_alignment.sh'

#location of the picard processing 
picard='step_03_picards_processing.sh'

#location of the bam qc 
bam_qc='step_00_qc_bam.sh'

#location of freebayes
freebayes='step_04_freebayes_variantcalling.sh'

#location of variants annotation 
variants_annotation='step_05_variants_annotation.sh'

#location of variants filteration 
variants_filt='step_06_variants_filt.sh'

#location of variants filteration
variants_stat='step_07_variants_stat.sh'

#location of vcf to table 
vcf_to_table='step_08_SNPs_to_table.sh' 

#location of grenedalf snps 
grenedalf_snps='step_09_grenedalf_snps.sh'

#-----------------------------------------------------

echo "${YELLOW}
     ......RUNNING THE POOLSEQ PIPELINE........ ${rest}
     " 

#make directories
echo "${GREEN}
..........Making pipeline directories ..........." 
bash $make_dir

#Run fastqc 
echo "${BLUE}
..........Running fastqc ..........." 
bash $fastqc 

#Run BWA indexing
echo "${GREEN}
..........Index genome If the genome is not indexed by BWA.........." 
bash $BWA_index

#Running BWA alignment 
echo "${RED}
..........Running BWA alignemnt.........." 
bash $BWA_align

#Running picard processing 
echo "${YELLOW}
..........Running picard processing.........." 
bash $picard

#Running bam quality check 
echo "${BLUE}
..........Running quality check of bam files.........." 
bash $bam_qc

#Running variant calling using freebayes
echo "${GREEN}
..........Running freebayes.........." 
bash $freebayes

#Running the variants annotation 
echo "${BLUE}
.........running snpEff variants annotation analyses ........." 
bash $variants_annotation

#Running the variants filteration  
echo "${RED}
.........running variants filteration ........." 
bash $variants_filt

#Obtaining variants stats
echo "${YELLOW}
..........Variants stats............" 
bash $variants_stat

#converts variants to a table 
echo "${BLUE}
...........convert variants to a table........." 
bash $vcf_to_table

#running population genetics analyses 
echo "${GREEN}
.........running grenedalf population genetics analyses on SNPs........." 
bash $grenedalf_snps

#running multiQC 
echo "${BLUE}
.........Multi QC........." 
multiqc ../0_QC --outdir ../0_QC/ 2> ../log/multiqc.log

echo "${YELLOW} 
.......ANALYSIS IS COMPLETE....."