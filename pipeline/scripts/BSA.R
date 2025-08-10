#!/usr/local/bin/Rscript

#install libraries
suppressMessages(library(tidyverse))
suppressMessages(library(docopt))
suppressMessages(library(QTLseqr))

# Define the usage message
usage <- "
Usage:
  script.R (--rawData <rawData>) (--highBulk <highBulk>) (--lowBulk <lowBulk>) (--lowBulk_size=<lbsize>) (--highBulk_size=<hbsize>) (--windowSize=<windowSize>) [--output=<output_path>] [--plot_path=<path_to_plots>]

Options:
  --rawData=<rawData>              Path to SNPs table 
  --highBulk=<highBulk>            High bulk sample name
  --lowBulk=<lowBulk>              Low bulk sample name
  --lowBulk_size=<lowbulksize>     Low bulk size (numeric)
  --highBulk_size=<highbulksize>   High bulk size (numeric)
  --windowSize=<windowSize>        Analysis window size (numeric)
  --output=<outputfile>            path to the results output directory (path)
  --plot_path=<path_to_plots>      path to the plots output directory (path)
"

#Parse command-line arguments using docopt
args <- docopt(usage)

# Convert to numeric - FIXED: correct variable assignment
hbsize <- as.numeric(args$highBulk_size)  # HIGH bulk size
lbsize <- as.numeric(args$lowBulk_size)   # LOW bulk size

#load the data 
data <- importFromGATK(args$rawData, 
                      highBulk = args$highBulk, 
                       lowBulk = args$lowBulk, 
                       chromList = NULL)

#print head of the data 
print(head(data))

#Filtering SNPs using default parameters
df_filt <-
  filterSNPs(
    SNPset = data,
    refAlleleFreq = 0.20,
    minTotalDepth = 100,
    maxTotalDepth = 400,
    minSampleDepth = 40,
    verbose=TRUE)

#Print head of the filtered data  
print(head(df_filt))

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis (
  SNPset = df_filt,
  windowSize = as.numeric(args$windowSize),
  popStruc = "F2",
  bulkSize = c(hbsize, lbsize),
  replications = 10000,
  intervals = c(95, 99)
)

#Run G prime stats
df_filt <- runGprimeAnalysis(SNPset = df_filt,
                              windowSize = as.numeric(args$windowSize),
                              outlierFilter = "deltaSNP", 
                              filterThreshold = 0.1)
print(head(df_filt))

#Get results 
results <- getQTLTable(SNPset = df_filt, method = "Gprime", alpha = 0.1, export = FALSE)

#print
print(head(results))

#Write table of the results 

#if output directory is not assigned use current working directory
if (is.null(args$output)) {
  output_dir <- getwd()
} else {
  output_dir <- args$output
}

#Assign an out file
raw_file <- file.path(output_dir, paste0(args$highBulk, "vs", args$lowBulk, "_", args$windowSize, "_WS_raw.csv"))

outfile <- file.path(output_dir, paste0(args$highBulk, "vs", args$lowBulk, "_alpha0.01.csv"))

#write table  
write.table(results, file = outfile, row.names = F, sep = ",") 

#write table  
write.table(df_filt, file = raw_file, row.names = F, sep = ",") 

##Producing and saving the plots 

#if plot directory is not assigned use current working directory - FIXED
if (is.null(args$plot_path)) {
  plot_dir <- getwd()
} else {
  plot_dir <- args$plot_path  # FIXED: was args$output
}

#Assign a SNP distribution file
SNP_dist <- file.path(plot_dir, paste0(args$highBulk, "vs", args$lowBulk, "_SNPs_dist.pdf"))

#Plotting the SNP/window distribution 
pdf(SNP_dist, width = 11, height = 8.5)  # Width > Height (landscape)
p1 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")
p1
dev.off()

#Assign a delta SNP distribution file
deltaSNP_dist <- file.path(plot_dir, paste0(args$highBulk, "vs", args$lowBulk, "_deltaSNP_dist.pdf"))

pdf(deltaSNP_dist, width = 11, height = 8.5)  # Width > Height (landscape)
p2 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)
p2
dev.off()

#Assign a Gprime distribution file
Gprime_dist <- file.path(plot_dir, paste0(args$highBulk, "vs", args$lowBulk, "_Gprime_dist.pdf"))

pdf(Gprime_dist, width = 11, height = 8.5)  # Width > Height (landscape)
p3 <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
p3
dev.off()

cat("BSA analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")
cat("Plots saved to:", plot_dir, "\n")