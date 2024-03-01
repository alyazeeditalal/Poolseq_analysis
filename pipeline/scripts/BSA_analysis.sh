#!/bin/bash

snps_Data='../5_freebayes_variants/pooled.continuous.freebayes.parallel.snpeff.snps.table.tsv'


#SNPs
Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_60 --highBulk_size=96 --lowBulk=DE17F2_15 --lowBulk_size=64 --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_60 --highBulk_size=96 --lowBulk=DE17F2_20 --lowBulk_size=77  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_60 --highBulk_size=96 --lowBulk=DE5F1_3 --lowBulk_size=64  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_60 --highBulk_size=96 --lowBulk=DE5F1_8 --lowBulk_size=96  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/



Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_75 --highBulk_size=40 --lowBulk=DE17F2_15 --lowBulk_size=64 --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_75 --highBulk_size=40 --lowBulk=DE17F2_20 --lowBulk_size=77  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_75 --highBulk_size=40 --lowBulk=DE5F1_3 --lowBulk_size=64  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL17F2_75 --highBulk_size=40 --lowBulk=DE5F1_8 --lowBulk_size=96  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/



Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_25 --highBulk_size=94 --lowBulk=DE17F2_15 --lowBulk_size=64 --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_25 --highBulk_size=94 --lowBulk=DE17F2_20 --lowBulk_size=77  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_25 --highBulk_size=94 --lowBulk=DE5F1_3 --lowBulk_size=64  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_25 --highBulk_size=94 --lowBulk=DE5F1_8 --lowBulk_size=96  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/


Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_30 --highBulk_size=40 --lowBulk=DE17F2_15 --lowBulk_size=64 --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_30 --highBulk_size=40 --lowBulk=DE17F2_20 --lowBulk_size=77  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_30 --highBulk_size=40 --lowBulk=DE5F1_3 --lowBulk_size=64  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/

Rscript scripts/BSA.R --rawData=$snps_Data --highBulk=AL5F1_30 --highBulk_size=40 --lowBulk=DE5F1_8 --lowBulk_size=96  --windowSiz=1e6  --output=../8_BSA_analysis/results/snps/ --plot_path=../8_BSA_analysis/plots/snps/