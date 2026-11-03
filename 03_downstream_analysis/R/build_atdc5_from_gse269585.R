# =============================================================================
# build_atdc5_from_gse269585.R -- Build atdc5.pub Seurat object from
#                                  GSE269585 (ATDC5 chondrogenic cell line)
#
# Wires up the data source for Figure 2f (plot_gene_divergence(), "Primary
# MAC vs. ATDC5" lollipop chart) -- previously defined in figure2.R but
# never called, since no ATDC5 dataset was available in this repo.
# archive/script.r documents the intended source (GSE269585) and originally
# fetched it live from GEO; a local raw count matrix
# (Rdata/RNAseq/GSE269585_Raw_Count.txt) is now available instead, so this
# script reads it directly rather than re-fetching from GEO.
#
# Only the 4 "37C_*" columns are used: per user confirmation, 37C is the
# only condition in GSE269585 (a temperature-response study, unrelated to
# this project's IL-1b axis) that is a valid baseline counterpart to the
# primary-MAC RNA-seq data (rna.pub) -- i.e. untreated/physiological-
# temperature ATDC5, comparable to rna.pub's control condition.
#
# Pipeline:
#   1. Read GSE269585_Raw_Count.txt; keep gene_name + the 4 37C_* columns.
#   2. Collapse gene_name duplicates (>1 Ensembl gene_id per symbol, 102
#      cases) by summing counts, same convention as
#      build_rna_pub_from_featurecounts.R.
#   3. Build the Seurat object (raw counts + LogNormalize "data" layer).
#   4. Save Rdata/pub/seurat_atdc5.pub.Rdata.
#
# Prerequisites: config.R
# Input:  Rdata/RNAseq/GSE269585_Raw_Count.txt
# Output: Rdata/pub/seurat_atdc5.pub.Rdata
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
})

GSE269585_RAW <- here::here("Rdata", "RNAseq", "GSE269585_Raw_Count.txt")
ATDC5_SAMPLES <- c("37C_1", "37C_2", "37C_3", "37C_4")

# =============================================================================
# STEP 1 -- Read + subset to the 37C baseline samples
# =============================================================================

raw <- read.delim(GSE269585_RAW, check.names = FALSE)
stopifnot("Expected columns missing from GSE269585_Raw_Count.txt" =
            all(c("gene_name", ATDC5_SAMPLES) %in% colnames(raw)))

message(sprintf("GSE269585_Raw_Count.txt: %d gene_id rows, %d distinct gene_name.",
                nrow(raw), length(unique(raw$gene_name))))

count_mat <- as.matrix(raw[, ATDC5_SAMPLES])
rownames(count_mat) <- raw$gene_name
storage.mode(count_mat) <- "integer"

message(sprintf("37C baseline assigned-read totals: %s",
                paste(sprintf("%s=%d", ATDC5_SAMPLES, colSums(count_mat)), collapse = ", ")))

# =============================================================================
# STEP 2 -- Collapse gene_name duplicates (>1 Ensembl gene_id per symbol)
# =============================================================================

dup_symbols <- unique(rownames(count_mat)[duplicated(rownames(count_mat))])
if (length(dup_symbols) > 0) {
  message(sprintf(
    "Collapsing %d gene symbols backed by >1 Ensembl gene_id (summing counts).",
    length(dup_symbols)
  ))
}
count_final <- rowsum(count_mat, group = rownames(count_mat))
storage.mode(count_final) <- "integer"
count_final <- count_final[, ATDC5_SAMPLES, drop = FALSE]  # rowsum reorders columns; restore

message(sprintf("Final gene-symbol x sample count matrix: %d genes x %d samples.",
                nrow(count_final), ncol(count_final)))

panel_present <- intersect(GENES_ALL, rownames(count_final))
message(sprintf("GEHTS 30-gene panel present in ATDC5 matrix: %d / %d (missing: %s)",
                length(panel_present), length(GENES_ALL),
                paste(setdiff(GENES_ALL, rownames(count_final)), collapse = ", ")))

# =============================================================================
# STEP 3 -- Build Seurat object
# =============================================================================

atdc5.pub <- CreateSeuratObject(counts = count_final, project = "atdc5.pub")
atdc5.pub$cell_line   <- "ATDC5"
atdc5.pub$temperature <- "37C"
atdc5.pub <- NormalizeData(atdc5.pub, verbose = FALSE)

print(atdc5.pub)

# =============================================================================
# STEP 4 -- Save
# =============================================================================

atdc5_pub_path <- here::here("Rdata", "pub", "seurat_atdc5.pub.Rdata")
save(atdc5.pub, file = atdc5_pub_path)
message("Saved: ", atdc5_pub_path)

save_session_info("build_atdc5_from_gse269585")
