# figure 2g: Primary vs ATDC5

# ==============================================================================
# Script Name: 01_fetch_and_prep_external_data.R
# Author: Wooseok Lee
# Date: 2026-03-03
# Description: Automatically downloads the GSE269585 comparative dataset from 
#              GEO, extracts the X37C target samples, maps gene symbols, and 
#              initializes the baseline Seurat object for downstream analysis.
# Dependencies: GEOquery, Seurat
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------
# Uncomment the following lines to install required dependencies if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
# if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
})

# Define global variables
GSE_ACCESSION <- "GSE269585"
TARGET_SAMPLES <- c("X37C_1", "X37C_2", "X37C_3", "X37C_4")
PROJECT_NAME <- "X37C_Comparison"

# Use relative paths for GitHub reproducibility (avoid absolute local paths like D:/...)
DATA_DIR <- "data/raw/geo_downloads"
OUTPUT_DIR <- "data/processed"

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Data Acquisition
# ------------------------------------------------------------------------------
message(sprintf("[%s] Fetching series-level data for %s...", Sys.time(), GSE_ACCESSION))

# Download supplementary files
file_info <- getGEOSuppFiles(GSE_ACCESSION, baseDir = DATA_DIR, makeDirectory = TRUE)
file_path <- rownames(file_info)[1]

if (is.null(file_path) || !file.exists(file_path)) {
  stop("Error: File download failed or file path is incorrect.")
}

# ------------------------------------------------------------------------------
# 3. Data Processing & Formatting
# ------------------------------------------------------------------------------
message(sprintf("[%s] Reading the master count matrix...", Sys.time()))
counts <- read.delim(gzfile(file_path), row.names = 1, header = TRUE)

# Verify target samples exist in the downloaded matrix
missing_samples <- setdiff(TARGET_SAMPLES, colnames(counts))
if (length(missing_samples) > 0) {
  stop(paste("Error: The following target samples were not found in the matrix:", 
             paste(missing_samples, collapse = ", ")))
}

message(sprintf("[%s] Subsetting target samples and formatting rownames...", Sys.time()))
# Extract only the 4 target count columns
clean_target_counts <- counts[, TARGET_SAMPLES]

# Map provided gene names to rownames, ensuring uniqueness
if ("gene_name" %in% colnames(counts)) {
  rownames(clean_target_counts) <- make.unique(as.character(counts$gene_name))
} else {
  warning("Column 'gene_name' not found. Retaining original Ensembl IDs as rownames.")
}

# ------------------------------------------------------------------------------
# 4. Seurat Object Initialization & Export
# ------------------------------------------------------------------------------
message(sprintf("[%s] Initializing Seurat object...", Sys.time()))
final_seurat <- CreateSeuratObject(counts = clean_target_counts, 
                                   project = PROJECT_NAME)

print(final_seurat)

# Save the object for the next script in the analytical pipeline
output_file <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_raw_seurat.rds"))
saveRDS(final_seurat, file = output_file)
message(sprintf("[%s] Success! Seurat object saved to: %s", Sys.time(), output_file))

# Optional: Print session info for reproducibility logging
# sessionInfo()