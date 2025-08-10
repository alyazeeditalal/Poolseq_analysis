# Pooled Sequencing (Pool-seq) Analysis Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: Bash](https://img.shields.io/badge/Language-Bash-green.svg)](https://www.gnu.org/software/bash/)
[![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-blue.svg)](https://en.wikipedia.org/wiki/Bioinformatics)

A bash-based bioinformatics pipeline for analyzing pooled sequencing data with maximum flexibility for customisation and integration across computing environments.

The pipeline was used to produce results in the following papers 

1. Genetic mapping of resistance: A QTL and associated polymorphism conferring resistance to alpha-cypermethrin in Anopheles funestus
T AL-Yazeedi, G Djuifo, L Mugenzi, A Muhammad, J Hearn, CS Wondji
bioRxiv, 2025.04. 16.649184 [link](https://www.biorxiv.org/content/10.1101/2025.04.16.649184v1.abstract)

## Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Additional Analyses](#additional-analyses)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

## Overview

This pipeline processes paired-end FASTQ files from pooled DNA samples through complete variant calling and population analysis:

**Pipeline Steps:**
1. **Quality Control** - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. **Alignment** - [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) 
3. **BAM Processing** - [GATK/Picard](https://github.com/broadinstitute/gatk)
4. **Variant Calling** - [FreeBayes](https://github.com/freebayes/freebayes)
5. **Variant Annotation** - [SnpEff](https://github.com/pcingola/SnpEff)
6. **Population Genetics** - [grenedalf](https://github.com/lczech/grenedalf)
7. **BSA QTLsqr** - [qtlseqr](https://github.com/bmansfeld/QTLseqr)

**Additional Analyses:**
- Bulk Segregant Analysis (BSA) for QTL mapping
- Cochran-Mantel-Haenszel (CMH) test for population differentiation popoolation2 (https://sourceforge.net/p/popoolation2/wiki/Tutorial/)

## Quick Start

```bash
# Clone and setup
git clone https://github.com/your-username/poolseq-pipeline.git
cd poolseq-pipeline
mkdir 1_data genome annotation input

# Add data and run
# (Add FASTQ files to 1_data/, reference to genome/, annotation to annotation/, pools.txt to input/)
cd pipeline && bash run_all.bash
```

## Installation

### Prerequisites
- Linux or macOS
- Conda/Mamba package manager

### Create Conda Environment

```bash
# Create new environment
conda create -n poolseq python=3.12

# Activate environment
conda activate poolseq

# Install core bioinformatics tools
mamba install -c bioconda bwa-mem2 samtools bcftools freebayes gatk4 snpeff snpsift grenedalf fastqc multiqc qualimap parallel 

# Install R packages for BSA analysis (optional)
mamba install -c conda-forge r-base r-essentials r-tidyverse r-ggplot2 r-devtools r-docopt

# Install population genetics tools (optional)
mamba install -c bioconda popoolation2
```

### 2. Verify Installation

```bash
# Test key tools
bwa-mem2 version
samtools --version
freebayes --version
gatk --version
grenedalf --version
```

## Directory Setup

### 1. Download Pipeline Scripts

```bash
# Clone or download pipeline scripts to a directory called 'pipeline'
cd /path/to/your/project
mkdir pipeline
# Copy all .sh scripts to pipeline/ directory
```

### 2. Create Required Directories

```bash
# Navigate to the cloned repository
cd poolseq-pipeline

# Create data and reference directories (pipeline/ already exists)
mkdir 1_data genome annotation input

# Final directory structure:
# poolseq-pipeline/
# ├── 1_data/           # FASTQ files
# ├── genome/           # Reference genome
# ├── annotation/       # GFF annotation files
# ├── input/            # Pool configuration
# ├── pipeline/         # Pipeline scripts (from GitHub)
# └── README.md         # Documentation
```

## Data Preparation

### 1. Prepare FASTQ Files

> **💡 Tip**: File naming is configurable via `config.sh` parameters

```bash
cd 1_data

# Place your paired-end FASTQ files here
# Naming convention is defined in config.sh:
# - fa_ex_1='_1.fq.gz'  (forward reads extension)
# - fa_ex_2='_2.fq.gz'  (reverse reads extension)

# Files should be named: SAMPLE_NAME + extension
# Examples:
# Pool_Alive_test_1.fq.gz    # Forward reads for Pool_Alive_test sample
# Pool_Alive_test_2.fq.gz    # Reverse reads for Pool_Alive_test sample
# Pool_dead_test_1.fq.gz     # Forward reads for Pool_dead_test sample
# Pool_dead_test_2.fq.gz     # Reverse reads for Pool_dead_test sample

# Alternative naming (modify config.sh accordingly):
# Sample.R1.fastq.gz / Sample.R2.fastq.gz  -> fa_ex_1='.R1.fastq.gz', fa_ex_2='.R2.fastq.gz'
# Sample_read1.fq.gz / Sample_read2.fq.gz  -> fa_ex_1='_read1.fq.gz', fa_ex_2='_read2.fq.gz'

# Important: All samples must follow the same naming convention!
```

### 2. Download Reference Genome

```bash
cd genome

# Download Anopheles gambiae reference genome
wget http://ftp.ensemblgenomes.org/pub/metazoa/release-57/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz
gunzip Anopheles_gambiae.AgamP4.dna.toplevel.fa.gz

# For testing with smaller data, use single chromosome:
# wget "http://ftp.ensemblgenomes.org/pub/metazoa/release-57/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna_sm.chromosome.2L.fa.gz"
# gunzip Anopheles_gambiae.AgamP4.dna_sm.chromosome.2L.fa.gz
# mv Anopheles_gambiae.AgamP4.dna_sm.chromosome.2L.fa Anopheles_gambiae.AgamP4.dna.toplevel.fa
```

### 3. Download Annotation Files

```bash
cd annotation

# Download GFF annotation
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-57/gff3/anopheles_gambiae/Anopheles_gambiae.AgamP4.57.gff3.gz
gunzip Anopheles_gambiae.AgamP4.57.gff3.gz

# Extract genes only
awk '{if($3=="gene") print $0}' Anopheles_gambiae.AgamP4.57.gff3 > gene.only.gff

# For single chromosome testing:
# awk '{if($3=="gene" && $1=="2L") print $0}' Anopheles_gambiae.AgamP4.57.gff3 > gene.only.gff
```

### 4. Create Pool Configuration

```bash
cd input

# Create pools.txt file with sample names and pool sizes
# Format: sample_name,pool_size
cat > pools.txt << EOF
Pool_Alive_test,50
Pool_dead_test,50
EOF
```

## Customization

### File Naming Conventions

The pipeline uses `config.sh` to define file naming patterns:

```bash
# Edit pipeline/config.sh to match your file naming:

# For files named Sample_1.fq.gz / Sample_2.fq.gz (default):
fa_ex_1='_1.fq.gz'
fa_ex_2='_2.fq.gz'

# For files named Sample.R1.fastq.gz / Sample.R2.fastq.gz:
fa_ex_1='.R1.fastq.gz'  
fa_ex_2='.R2.fastq.gz'

# For files named Sample_read1.fq.gz / Sample_read2.fq.gz:
fa_ex_1='_read1.fq.gz'
fa_ex_2='_read2.fq.gz'
```

### Species and Reference Genome

```bash
# For different species, update these in config.sh:

# Reference genome location
loc_genome='../genome/YOUR_SPECIES.fa'

# SnpEff database name (check available databases)
genome_DB='YOUR_SPECIES_DATABASE'

# Annotation file
genes_only='../annotation/YOUR_SPECIES_genes.gff'
```

### Analysis Parameters

```bash
# Population genetics parameters in config.sh:
window_sliding_width='50000'    # Window size for sliding window analysis
window_sliding_stride='25000'   # Step size between windows
min_cov='10'                   # Minimum coverage filter
min_count='2'                  # Minimum allele count
```

## Running the Pipeline

### 1. Activate Environment and Run

```bash
conda activate poolseq
cd pipeline

# Run complete pipeline
bash run_all.bash
```

The pipeline executes these steps automatically:
1. **QC**: FastQC on raw reads
2. **Indexing**: BWA-MEM2 genome indexing  
3. **Alignment**: Read alignment to reference
4. **Processing**: BAM file processing (sort, mark duplicates, add read groups)
5. **QC**: BAM file quality control
6. **Variant Calling**: FreeBayes variant calling
7. **Annotation**: SnpEff variant annotation
8. **Filtering**: Variant filtering and separation
9. **Statistics**: Variant statistics
10. **Conversion**: VCF to table format
11. **Population Genetics**: Fst, diversity, and frequency analysis
12. **MultiQC**: Combined quality report

### 2. Run Individual Steps (Optional)

> **💡 Flexibility**: The bash design allows running individual pipeline components

```bash
# Run only quality control
bash step_00_qc_fasta_files.sh

# Run only alignment steps
bash step_01_BWA_index_genome.sh
bash step_02_BWA_alignment.sh

# Run only variant calling
bash step_04_freebayes_variantcalling.sh

# Run only population genetics
bash step_09_grenedalf_snps.sh
```

## Configuration

### Update config.sh

> **⚠️ Important**: Update the SnpSift oneperline script path for your environment

```bash
cd pipeline

# Find the correct path for SnpSift oneperline script
find $CONDA_PREFIX -name "vcfEffOnePerLine.pl" 2>/dev/null

# Edit config.sh and update the oneperline path:
# oneperline='/path/to/vcfEffOnePerLine.pl'
```

## Additional Analyses

### CMH Test (Population Differentiation)
Identifies significantly differentiated SNPs between populations using the Cochran-Mantel-Haenszel test.

**Requirements:** Completed main pipeline, grenedalf, popoolation2

**Example:**
```bash
cd pipeline
bash cmh_analysis.sh
```

**Output:** `7_population_genetics/population.cmh` containing SNPs with p-values

**Customization:**
```bash
# Edit cmh_analysis.sh for multiple population comparisons:
--population 1-2,1-3,2-3

# Adjust filtering:
--min-count 12 --min-coverage 50 --max-coverage 200
```

### BSA Analysis (Bulk Segregant Analysis)
Maps quantitative trait loci (QTL) by analyzing allele frequency differences between extreme phenotype pools.

**Requirements:** Completed main pipeline, R packages (tidyverse, ggplot2, docopt, QTLseqr)

**Example:**
```bash
cd pipeline
bash BSA_analysis.sh
```

**Outputs:**
- `8_BSA_analysis/results/snps/*.csv` - QTL analysis results
- `8_BSA_analysis/plots/snps/*.pdf` - Statistical plots (SNP distribution, Δ(SNP-index), G' statistic)

**Customization:**
```bash
# Edit BSA_analysis.sh for different samples:
--highBulk=Resistant_pool --highBulk_size=100
--lowBulk=Susceptible_pool --lowBulk_size=100  
--windowSize=1e6

# Multiple comparisons:
Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=Pool1 --lowBulk=Pool2 ...
Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=Pool1 --lowBulk=Pool3 ...
```

### Output File Structure

```
poolseq-pipeline/
├── 0_QC/                    # Quality control reports
│   ├── fasta_qc/           # FastQC reports
│   ├── bams_QC/            # BAM quality reports
│   └── multiqc_report.html # Combined QC report
├── 2_mapping/              # BAM alignment files
├── 3_processing/           # Processed BAM files
├── 5_freebayes_variants/   # Variant files (VCF, tables)
├── 7_population_genetics/  # Population genetics results
├── 8_BSA_analysis/         # BSA results (if run)
│   ├── results/snps/       # CSV result files
│   └── plots/snps/         # PDF plots
└── log/                    # Log files for all steps
```

###  Key Result Files

- **Variants**: `5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.vcf`
- **SNP Table**: `5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.table.tsv`
- **Population Stats**: `7_population_genetics/snps_windowed_*.csv`
- **QC Report**: `0_QC/multiqc_report.html`

## Citation

If you use this pipeline, please cite the tools and our papers:
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [GATK/Picard](https://github.com/broadinstitute/gatk)  
- [FreeBayes](https://github.com/freebayes/freebayes)
- [SnpEff](https://github.com/pcingola/SnpEff)
- [grenedalf](https://github.com/lczech/grenedalf)
- SAMtools, BCFtools, MultiQC 

