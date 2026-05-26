# =============================================================================
# normalization.R — Seurat / SCTransform normalization for GEHTS data
#
# Inputs  (h5ad files written by  python main.py <dir> --save):
#   data_sorted/control&inflam/merged.h5ad
#   data_sorted/single_only/merged.h5ad
#   data_sorted/combi_only/merged.h5ad
#
# Outputs (written to ../Rdata/pub/ for use by 03_downstream_analysis):
#   mac.sct.RData, sin.sct.RData, cmb.sct.RData
#
# Note: sceasy uses reticulate to read h5ad files. If reticulate does not
# find the correct Python automatically, set RETICULATE_PYTHON before
# sourcing this script:
#   Sys.setenv(RETICULATE_PYTHON = "/path/to/your/conda/envs/gehts/bin/python")
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(sceasy)
  library(dplyr)
  library(harmony)
})

set.seed(123)

# =============================================================================
# Helper functions
# =============================================================================

extract_filename <- function(file_name) substr(file_name, 1, 6)

load_data <- function(file_path) {
  seurat_obj <- sceasy::convertFormat(file_path, from = "anndata", to = "seurat")
  Idents(seurat_obj) <- "id"
  seurat_obj$nCount_RNA   <- Matrix::colSums(seurat_obj@assays$RNA@counts)
  seurat_obj$nFeature_RNA <- Matrix::colSums(seurat_obj@assays$RNA@counts > 0)
  seurat_obj
}

preprocess_data <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 0 & nCount_RNA > 0)
  seurat_obj@meta.data$file <- sapply(seurat_obj@meta.data$file, extract_filename)
  seurat_obj
}

# =============================================================================
# Step 1. Control & inflammation reference (mac)
# =============================================================================

mac <- load_data(here("data_sorted", "control&inflam", "merged.h5ad"))
mac <- preprocess_data(mac)
mac$drug_name <- mac$drug_condition
mac$source    <- "mac"

# =============================================================================
# Step 2. Single-drug screen (sin)
# =============================================================================

drg <- load_data(here("data_sorted", "single_only", "merged.h5ad"))
drg <- preprocess_data(drg)
drg <- subset(drg, subset = drug_condition != 'inflammatory')

split_drg     <- strsplit(drg$drug_condition, "_", fixed = TRUE)
drg$drug_name <- sapply(split_drg, `[[`, 1)
drg$dose      <- sapply(split_drg, `[[`, 2)
drg$source    <- "drg"

# =============================================================================
# Step 3. Combination-drug screen (cmb)
# =============================================================================

drg2        <- load_data(here("data_sorted", "combi_only", "merged.h5ad"))
drg2$file   <- extract_filename(drg2$file)
drg2$source <- "cmb"

`%notin%` <- Negate(`%in%`)

# Separate true combination wells (both slots filled) from single-drug wells
keep_combi                 <- !grepl("^0&|&999$", drg2@meta.data$drug_condition)
drg2@meta.data$keep_combi <- keep_combi
cmb <- subset(drg2, subset = keep_combi)
cmb$drug_condition <- as.character(droplevels(as.factor(cmb$drug_condition)))

split_cmb  <- strsplit(cmb$drug_condition, "&", fixed = TRUE)
split_cmb1 <- lapply(split_cmb, function(x) strsplit(x[1], "_", fixed = TRUE)[[1]])
split_cmb2 <- lapply(split_cmb, function(x) strsplit(x[2], "_", fixed = TRUE)[[1]])
cmb$drug_name1 <- sapply(split_cmb1, `[[`, 1)
cmb$dose1      <- sapply(split_cmb1, `[[`, 2)
cmb$drug_name2 <- sapply(split_cmb2, `[[`, 1)
cmb$dose2      <- sapply(split_cmb2, `[[`, 2)

# Retrieve single-drug wells embedded in the combination batch
keep_single                  <- grepl("^0&", drg2@meta.data$drug_condition)
drg2@meta.data$keep_single  <- keep_single
single_cmb <- subset(drg2, subset = keep_single)

unwanted   <- c('0&0', '0&10', '0&27', '0&999')
single_cmb <- subset(single_cmb, subset = drug_condition %notin% unwanted)
single_cmb$drug_condition <- as.character(
  droplevels(as.factor(single_cmb$drug_condition))
)
split_sc              <- strsplit(single_cmb$drug_condition, "&", fixed = TRUE)
single_cmb$drug_condition <- sapply(split_sc, `[[`, 2)
split_sc2             <- lapply(single_cmb$drug_condition,
                                function(x) strsplit(x, "_", fixed = TRUE)[[1]])
single_cmb$drug_name  <- sapply(split_sc2, `[[`, 1)
single_cmb$dose       <- sapply(split_sc2, `[[`, 2)
single_cmb$source     <- "drg"

# =============================================================================
# Step 4. Merge all batches and apply SCTransform
# =============================================================================

mac.pub <- mac
sin.pub <- merge(drg, single_cmb)
cmb.pub <- cmb

iss_all <- Reduce(merge, list(mac.pub, cmb.pub, sin.pub))
iss_all$file <- extract_filename(iss_all$file)
iss_all <- subset(iss_all, subset = nFeature_RNA > 0 & nCount_RNA > 0)

iss_SCT <- SCTransform(
  iss_all,
  vst.flavor            = "v2",
  verbose               = FALSE,
  return.only.var.genes = FALSE,
  min_cells             = 3
)

# =============================================================================
# Step 5. Separate by source and save
# =============================================================================

mac.sct <- subset(iss_SCT, subset = source == "mac")
sin.sct <- subset(iss_SCT, subset = source == "drg")
cmb.sct <- subset(iss_SCT, subset = source == "cmb")

out_dir <- here("..", "Rdata", "pub")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save(mac.sct, file = file.path(out_dir, "mac.sct.RData"))
save(sin.sct, file = file.path(out_dir, "sin.sct.RData"))
save(cmb.sct, file = file.path(out_dir, "cmb.sct.RData"))

message("Saved mac.sct, sin.sct, cmb.sct to ", out_dir)
