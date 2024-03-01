#!/usr/local/bin/Rscript

# Libraries to install
libraries_to_install <- c("devtools", "tidyverse", "ggplot2",
                          "docopt")

# Function to check if a library is installed
is_library_installed <- function(lib) {
  return(requireNamespace(lib, quietly = TRUE))
}

# Install libraries if they don't exist
for (lib in libraries_to_install) {
  if (!is_library_installed(lib)) {
    install.packages(lib, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

# Install QTLseqr from GitHub
if (!is_library_installed("QTLseqr")) {
 
 
}

# Print installation status
for (lib in libraries_to_install) {
  if (is_library_installed(lib)) {
    cat(paste("Library", lib, "is installed.\n"))
  } else {
    cat(paste("Failed to install", lib, ".\n"))
  }
}
