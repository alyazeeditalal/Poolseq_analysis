#!/bin/bash

#----------------------------------------------------------------
# This script will be used to convert the vcf file to a sync file 
#----------------------------------------------------------------
#location of config file
source config.sh

#----------------------------------------------------------------
# windowed pop_gen
#----------------------------------------------------------------

#Calculating Fst, allele frequency difference and diversity between all cohorts at different window sizes 

#Kofler fst 
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type sliding \
--method kofler \
--window-sliding-width $window_sliding_width \
--window-sliding-stride $window_sliding_stride \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indel_windowed_kofler \
--log-file $log_grenedalf/indel_windowed_kofler.log 

#unbiased-nei fst
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type sliding \
--method unbiased-nei \
--window-sliding-width $window_sliding_width \
--window-sliding-stride $window_sliding_stride \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indels_windowed_unbiased-nei \
--log-file $log_grenedalf/indels_windowed_unbiased-nei.log 


#Diversity
grenedalf diversity \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type sliding \
--window-sliding-width $window_sliding_width \
--window-sliding-stride $window_sliding_stride \
--filter-sample-min-coverage $min_cov \
--filter-sample-min-count $min_count \
--pool-sizes $pool_sizes \
--file-prefix indels_windowed_ \
--log-file $log_grenedalf/indels_windowed_diversity.log 

#----------------------------------------------------------------
# single pop_gen
#----------------------------------------------------------------
#frequency
grenedalf frequency \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--file-prefix indels_ \
--write-sample-counts \
--write-sample-coverage \
--write-sample-alt-freq \
--write-total-counts \
--write-total-coverage \
--write-total-frequency \
--log-file $log_grenedalf/indels_freq.log 

#Kofler fst 
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type single \
--method kofler \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indels_windowed_kofler \
--log-file $log_grenedalf/indels_single_kofler.log 

#unbiased-nei fst
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type single \
--method unbiased-nei \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indels_windowed_unbiased-nei \
--log-file $log_grenedalf/indels_single_unbiased-nei.log 

#----------------------------------------------------------------
# Gene_wise pop_gen
#----------------------------------------------------------------
#Kofler fst 
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type regions \
--method kofler \
--window-region-gff $genes_only \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indels_genes_kofler \
--log-file $log_grenedalf/indels_genes_kofler.log 

#unbiased-nei fst
grenedalf fst \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type regions \
--method unbiased-nei \
--window-region-gff $genes_only \
--filter-sample-min-coverage $min_cov \
--filter-total-only-biallelic-snps \
--pool-sizes $pool_sizes \
--file-prefix indels_genes_unbiased-nei \
--log-file $log_grenedalf/indels_genes_unbiased-nei.log 

#Diversity
grenedalf diversity \
--threads $t \
--vcf-path $indel_input \
--out-dir $out_dir \
--window-type regions \
--window-region-gff $genes_only \
--filter-sample-min-coverage $min_cov \
--filter-sample-min-count $min_count \
--pool-sizes $pool_sizes \
--file-prefix indels_genes_ \
--log-file $log_grenedalf/indels_genes_diversity.log 