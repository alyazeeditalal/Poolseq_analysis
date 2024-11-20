#!/usr/local/bin/Rscript

#install libraries
suppressMessages(library(tidyverse))
suppressMessages(library(docopt))
suppressMessages(library(glue))
suppressMessages(library(QTLseqr))
suppressMessages(library(ggpubr))
# Define the usage message

usage <- "
Usage:
  script.R (--rawData <rawData>) (--highBulk <highBulk>) (--lowBulk <lowBulk>) (--lowBulk_size=<lbsize>) (--highBulk_size=<hbsize>) (--windowSize=<windowSize>) [--output=<output_path>] [--plot_path=<path_to_plots>]

Options:
  --rawData=<rawData>              Path to SNPs table 
  --highBulk=<highBulk>            High bulk sample name
  --lowBulk=<lowBulk>              Low bulk sample name
  --lowBulk_size=<lowbulksize>     Low bulk size (numeric)
  --highBulk_size=<highbulksize>   Hight bulk size (numeric)
  --windowSize=<windowSize>        Analaysis window size (numeric)
  --output=<outputfile>            path to the results output directory (path)
  --plot_path=<path_to_plots>      path to the plots potput directory (path)
"

#Parse command-line arguments using docopt
args <- docopt(usage)

# Convert hbsize and lbsize to numeric
hbsize <- as.numeric(args$lowBulk_size)
lbsize <- as.numeric(args$highBulk_size)

#load the data 
data <- importFromGATK(args$rawData, 
                      highBulk = args$highBulk, 
                       lowBulk = args$lowBulk, 
                       chromList = NULL)

#print head of the data 
print(head(data))

#save(data, file = "~/Documents/LSTM_Postdoc/FANGvsFUMOZ_pool_seq/All.data.RData")

#Filtering SNPs using default parameters as in the 
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
                              outlierFilter = "deltaSNP")
print(head(df_filt))

#Get results 
results <- getQTLTable(SNPset = df_filt, method = "Gprime", alpha = 0.01, export = FALSE)

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

#if output directory is not assigned use current working directory
if (is.null(args$plot_path)) {
  plot_dir <- getwd()
} else {
  plot_dir <- args$plot_path
}

#Assign a file name for the plots file
SNP_plot <- file.path(plot_dir, paste0(args$highBulk, "vs", args$lowBulk, "_SNPs_plot.pdf"))

#Plotting the SNP/window distribution 
p1 <- plotQTLStats(SNPset = df_filt, var = "nSNPs")

p2 <- plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

p3 <- plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)

p4 <- plotQTLStats(SNPset = df_filt, var = "negLog10Pval", plotThreshold = TRUE, q = 0.01)

#combining the plots into one window 
p_all <- ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4, labels = c("A", "B", "C", "D"), common.legend	= T)

#saving the final plot result
pdf(SNP_plot, width = 11, height = 8.5)  # Width > Height (landscape)
p_all
dev.off()

