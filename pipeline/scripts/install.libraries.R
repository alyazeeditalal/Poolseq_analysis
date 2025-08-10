#!/usr/local/bin/Rscript

# Simplified R library installer - assumes most packages installed via conda
cat("Checking R packages for BSA analysis...\n")

# Function to check if a library is installed
is_library_installed <- function(lib) {
  return(requireNamespace(lib, quietly = TRUE))
}

# Check conda-installed packages
conda_packages <- c("tidyverse", "ggplot2", "devtools", "docopt")

cat("Checking conda-installed packages:\n")
for (lib in conda_packages) {
  if (is_library_installed(lib)) {
    cat(paste("✓", lib, "available\n"))
  } else {
    cat(paste("✗", lib, "missing - install with: mamba install -c conda-forge r-", lib, "\n", sep=""))
  }
}

# Only install QTLseqr from GitHub (not available in conda)
cat("\nChecking QTLseqr (GitHub package):\n")
if (!is_library_installed("QTLseqr")) {
  cat("Installing QTLseqr from GitHub...\n")
  tryCatch({
    devtools::install_github("bmansfeld/QTLseqr")
    cat("✓ QTLseqr installed successfully\n")
  }, error = function(e) {
    cat("✗ QTLseqr installation failed:", e$message, "\n")
  })
} else {
  cat("✓ QTLseqr already available\n")
}

cat("\nR package check completed!\n")