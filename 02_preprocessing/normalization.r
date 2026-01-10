# GE-HTS project preprocessing pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee, SNU

library(Seurat)
library(SeuratDisk)
path_to_python <- "C:/Users/name/AppData/Local/yourconda/envs/envname"
Sys.setenv(RETICULATE_PYTHON=path_to_python)
library(sceasy)
library(dplyr)
library(harmony)
library(gridExtra)
library(ggplot2) # the latest version 3.5.0 released at 24.02.23 is desirable.

# For reproducibility
set.seed(123)


###################################################################################
### Define functions ###
extract_filename <- function(file_name) {
  substr(file_name, 1, 6)
}

# Function to load data
load_data <- function(file_path) {
  seurat_obj <- sceasy::convertFormat(file_path, from = "anndata", to = "seurat")
  Idents(seurat_obj) <-"id"
  # RNA counts per well
  total_counts <- Matrix::colSums(seurat_obj@assays$RNA@counts)
  seurat_obj$nCount_RNA <- total_counts
  # RNA types more than 0 per well
  nfeatures <- Matrix::colSums(seurat_obj@assays$RNA@counts>0)
  seurat_obj$nFeature_RNA <- nfeatures
  return(seurat_obj)
}

# Function to preprocess data
preprocess_data <- function(seurat_obj) {
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 0 & nCount_RNA > 0)
  seurat_obj@meta.data$file <- sapply(seurat_obj@meta.data$file, extract_filename)
  return(seurat_obj)
}
##################################################################

# Load raw gene count file into R variables
# set input file directory
f_dir1 <- "./data_sorted/control&inflam/merged.h5ad" 
f_dir2 <- "./data_sorted/single_only/merged.h5ad"
f_dir3 <- "./data_sorted/combi_only/merged.h5ad"

### Step 1. Control & inflammation only
mac <- load_data(f_dir1)
mac <- preprocess_data(mac)
mac$drug_name <- mac$drug_condition

### Step 2. Single drug file 
drg <- load_data(f_dir2)
drg <- preprocess_data(drg)
drg <- subset(drg, subset = drug_condition != 'inflammatory') # remove empty wells

# split drug_name and dose from drug_condition
split_conditions <- strsplit(drg$drug_condition, "_", fixed = TRUE)
drg$drug_name <- sapply(split_conditions, function(x) {return(x[1])})
drg$dose <- sapply(split_conditions, function(x) {return(x[2])})

### Step 3. combi drg
drg2 <- load_data(f_dir3)
drg2$file <- extract_filename(drg2$file)

# Processing meta data
# split with '&' first, then split with '_' second for two drugs
split_conditions <- strsplit(drg2$drug_condition, "&", fixed = TRUE)
split_conditions1 <- sapply(split_conditions, function(x) {return( strsplit(x[1], '_', fixed = TRUE))})
split_conditions2 <- sapply(split_conditions, function(x) {return( strsplit(x[2], '_', fixed = TRUE))})
drg2$drug_name1 <- sapply(split_conditions1, function(x) {return(x[1])})
drg2$dose1 <- sapply(split_conditions1, function(x) {return(x[2])})
drg2$drug_name2 <- sapply(split_conditions2, function(x) {return(x[1])})
drg2$dose2 <- sapply(split_conditions2, function(x) {return(x[2])})

# single dose drugs sepration from cmb and merging with single batch data
# Prepare single drugs per dose = 10, 0.1
keep_combi <- !grepl("^0&|&999$", drg2@meta.data$drug_condition)
# Update the metadata
drg2@meta.data$keep_combi <- keep_combi
# Subset the Seurat object to exclude these cells
cmb <- subset(drg2, subset = keep_combi)
cmb$drug_condition <- as.factor(cmb$drug_condition)
cmb@meta.data$drug_condition <- droplevels(cmb@meta.data$drug_condition)
cmb$drug_condition <- as.character(cmb$drug_condition)

# Retrieve single drugs from cmbsub
keep_single <- grepl("^0&", drg2@meta.data$drug_condition) # for retrieving single drug treated cells from combination batches
drg2@meta.data$keep_single <- keep_single 
single_cmb <- subset(drg2, subset = keep_single) # for retrieving single drug treated cells from combination batches
unwanted_drugs <- c('0&0', '0&10', '0&27', '0&999')
'%notin%' <- Negate('%in%') # define negate of %in%
single_cmb <- subset(single_cmb, subset = drug_condition %notin% unwanted_drugs) # subset some unwanted drug conditions
single_cmb@meta.data$drug_condition <- as.factor(single_cmb@meta.data$drug_condition)
single_cmb@meta.data$drug_condition <- droplevels(single_cmb@meta.data$drug_condition) # remove levels with 0 value
single_cmb@meta.data$drug_condition <- as.character(single_cmb@meta.data$drug_condition) # make sure the levels are characters
split_conditions <- strsplit(single_cmb$drug_condition, "&", fixed = TRUE) # split 0& from the single drug
single_cmb$drug_condition <- sapply(split_conditions, function(x) {return(x[2])}) # return only the single drug name
split_conditions <- sapply(single_cmb$drug_condition, function(x) {return( strsplit(x[1], '_', fixed = TRUE))}) # split drug name and drug dose
single_cmb$drug_name <- sapply(split_conditions, function(x) {return(x[1])}) # assign drug name
single_cmb$dose <- sapply(split_conditions, function(x) {return(x[2])}) # assign drug dose

# split with '&' first, then split with '_' second for two drugs
split_conditions <- strsplit(cmb$drug_condition, "&", fixed = TRUE)
split_conditions1 <- sapply(split_conditions, function(x) {return( strsplit(x[1], '_', fixed = TRUE))})
split_conditions2 <- sapply(split_conditions, function(x) {return( strsplit(x[2], '_', fixed = TRUE))})
cmb$drug_name1 <- sapply(split_conditions1, function(x) {return(x[1])})
cmb$dose1 <- sapply(split_conditions1, function(x) {return(x[2])})
cmb$drug_name2 <- sapply(split_conditions2, function(x) {return(x[1])})
cmb$dose2 <- sapply(split_conditions2, function(x) {return(x[2])})

# datasets before normalization
mac.pub <- mac
sin.pub <- merge(drg, single_cmb)
cmb.pub <- cmb

### Step 4. Normalization
# merge all datasets
mac.pub$source <- 'mac'
cmb.pub$source <- 'cmb'
sin.pub$source <- 'drg'

mac.pub$source <- 'ref'
cmb.pub$source <- 'combi'
sin.pub$source <- 'single'

iss_all <- merge(mac.pub,cmb.pub)
iss_all <- merge(iss_all,sin.pub)
iss_all$file <- extract_filename(iss_all$file)
iss_all <- subset(iss_all, subset = nFeature_RNA > 0 & nCount_RNA > 0)

iss_SCT <- SCTransform(iss_all, vst.flavor='v2',verbose=FALSE, return.only.var.genes = FALSE, min_cells = 3) # This data is the evenly normalizaed across all data we need in major analysis.

# Now seprate data by their treatment
mac.sct <- subset(iss_SCT, subset= source == 'mac')
sin.sct <- subset(iss_SCT, subset= source == 'drg')
cmb.sct <- subset(iss_SCT, subset= source == 'cmb')
