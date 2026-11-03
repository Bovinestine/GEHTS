# =============================================================================
# build_rna_pub_from_featurecounts.R -- Rebuild rna.pub from raw featureCounts
#
# This script rebuilds rna.pub directly from 6 raw featureCounts outputs added to
# Rdata/RNAseq/ (ctrl[23]-control, il1b[23]-il1b -- n=3 control + n=3 il1b).
# These are the same 6 samples underlying the already-checked-in
# "control vs il1b_GeneDifferential_vst_count.csv" (identical column-name
# convention: "<sample>-control" / "<sample>-il1b"), but here we start from
# raw counts instead of an externally VST-transformed matrix, so
# normalization is fully reproducible in this repo.
#
# Pipeline:
#   1. Read + merge the 6 raw featureCounts outputs (Ensembl gene IDs,
#      GRCm38.75 annotation).
#   2. Map Ensembl gene ID -> mouse gene symbol via msigdbr's mouse-native
#      gene tables (offline; no org.Mm.eg.db/biomaRt/internet access in this
#      environment). Many-to-one Ensembl-ID collisions on the same symbol
#      are collapsed by summing counts.
#   3. Build the Seurat object (raw counts + LogNormalize "data" layer),
#      cells named control/control.1/control.2/inflammatory/inflammatory.1/
#      inflammatory.2 -- matching the cell-naming convention already used by
#      rna.ref (seurat_rna.ref.Rdata) elsewhere in this pipeline.
#   4. QC only: DESeq2 vst() + sample PCA (does not feed back into rna.pub).
#   5. Save Rdata/pub/seurat_rna.pub.Rdata (previous file backed up once to
#      seurat_rna.pub.legacy_n2.Rdata alongside it).
#
# Prerequisites: config.R
# Inputs: Rdata/RNAseq/{ctrl,ctrl2,ctrl3,il1b,il1b2,il1b3}-{control,il1b}
# Outputs: Rdata/pub/seurat_rna.pub.Rdata
#          03_downstream_analysis/output/rnaseq_build/rna_pub_featurecounts_QC_pca.pdf
#          03_downstream_analysis/output/rnaseq_build/rna_pub_featurecounts_methods_note.md
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
})

RNASEQ_RAW_DIR <- here::here("Rdata", "RNAseq")
OUTPUT_BUILD   <- here::here("03_downstream_analysis", "output", "rnaseq_build")
dir.create(OUTPUT_BUILD, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# STEP 1 -- Read + merge raw featureCounts outputs
# =============================================================================

sample_files <- c(
  control        = "ctrl-control",
  control.1      = "ctrl2-control",
  control.2      = "ctrl3-control",
  inflammatory   = "il1b-il1b",
  inflammatory.1 = "il1b2-il1b",
  inflammatory.2 = "il1b3-il1b"
)

read_featurecounts <- function(path) {
  df <- read.delim(path, skip = 1, header = TRUE, check.names = FALSE)
  counts <- df[[ncol(df)]]
  names(counts) <- df$Geneid
  counts
}

count_list <- lapply(file.path(RNASEQ_RAW_DIR, sample_files), read_featurecounts)
names(count_list) <- names(sample_files)

gene_ids <- Reduce(intersect, lapply(count_list, names))
stopifnot("All 6 featureCounts files must share the same gene universe" =
            length(gene_ids) == length(count_list[[1]]))

count_mat <- do.call(cbind, lapply(count_list, function(x) x[gene_ids]))
colnames(count_mat) <- names(sample_files)
rownames(count_mat) <- gene_ids
storage.mode(count_mat) <- "integer"

message(sprintf("Merged featureCounts matrix: %d Ensembl genes x %d samples.",
                nrow(count_mat), ncol(count_mat)))
message(sprintf("Assigned-read totals: %s",
                paste(sprintf("%s=%d", colnames(count_mat), colSums(count_mat)), collapse = ", ")))

# =============================================================================
# STEP 2 -- Ensembl gene ID -> mouse gene symbol (offline, via msigdbr)
# =============================================================================
# msigdbr's mouse-native ("MM") gene tables, pooled across every MSigDB
# collection, are used as an offline Ensembl-ID -> gene-symbol map since
# org.Mm.eg.db/biomaRt are not available and this environment has no
# internet access. Coverage is therefore MSigDB-curated-gene-limited, not
# whole-genome; genes without a curated symbol are dropped (see methods
# note for exact counts).

msig_all <- msigdbr(species = "Mus musculus", db_species = "MM")
id_map <- unique(as.data.frame(msig_all[, c("ensembl_gene", "gene_symbol")]))
id_map <- id_map[!is.na(id_map$ensembl_gene) & id_map$ensembl_gene != "", ]

mapped <- id_map[id_map$ensembl_gene %in% rownames(count_mat), ]
n_genes_in    <- nrow(count_mat)
n_genes_mapped <- length(unique(mapped$ensembl_gene))
n_genes_dropped <- n_genes_in - n_genes_mapped

message(sprintf(
  "Ensembl ID -> gene symbol (msigdbr, mouse-native, all collections): %d / %d genes mapped; %d dropped (no curated MSigDB symbol).",
  n_genes_mapped, n_genes_in, n_genes_dropped
))

count_mapped <- count_mat[mapped$ensembl_gene, , drop = FALSE]
rownames(count_mapped) <- mapped$gene_symbol

dup_symbols <- unique(rownames(count_mapped)[duplicated(rownames(count_mapped))])
if (length(dup_symbols) > 0) {
  message(sprintf(
    "Collapsing %d gene symbols backed by >1 Ensembl ID (summing counts across IDs).",
    length(dup_symbols)
  ))
}

count_final <- rowsum(count_mapped, group = rownames(count_mapped))
storage.mode(count_final) <- "integer"
count_final <- count_final[, names(sample_files), drop = FALSE]  # rowsum reorders columns; restore

# msigdbr reports current MGI-approved symbols, which have drifted from a
# few of the GEHTS chip panel's legacy probe-design symbols (GENES_ALL,
# config.R). Alias known renames back to the panel symbol so the rebuilt
# matrix stays directly comparable to mac.sct/macSebastian without editing
# the panel definition itself (which is fixed to the physical chip design).
LEGACY_SYMBOL_ALIASES <- c(Hprt1 = "Hprt")
alias_hits <- intersect(rownames(count_final), names(LEGACY_SYMBOL_ALIASES))
if (length(alias_hits) > 0) {
  new_names <- LEGACY_SYMBOL_ALIASES[alias_hits]
  stopifnot("Alias target symbol must not already exist in the matrix" =
              !any(new_names %in% setdiff(rownames(count_final), alias_hits)))
  rownames(count_final)[match(alias_hits, rownames(count_final))] <- new_names
  message(sprintf("Aliased legacy panel symbol(s): %s",
                  paste(sprintf("%s -> %s", alias_hits, new_names), collapse = ", ")))
}

message(sprintf("Final gene-symbol x sample count matrix: %d genes x %d samples.",
                nrow(count_final), ncol(count_final)))

# Sanity check: 30-gene GEHTS panel should be (almost) fully present
panel_present <- intersect(GENES_ALL, rownames(count_final))
message(sprintf("GEHTS 30-gene panel present in rebuilt matrix: %d / %d (missing: %s)",
                length(panel_present), length(GENES_ALL),
                paste(setdiff(GENES_ALL, rownames(count_final)), collapse = ", ")))

# =============================================================================
# STEP 3 -- Build Seurat object
# =============================================================================

condition <- ifelse(grepl("^control", colnames(count_final)), COND_CONTROL, COND_DISEASE)

rna.pub <- CreateSeuratObject(counts = count_final, project = "rna.pub")
rna.pub$drug_condition <- condition
rna.pub <- NormalizeData(rna.pub, verbose = FALSE)

print(rna.pub)

# =============================================================================
# STEP 4 -- QC only: DESeq2 vst() + sample PCA (does not feed back into rna.pub)
# =============================================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_final,
  colData   = data.frame(condition = condition, row.names = colnames(count_final)),
  design    = ~condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message(sprintf("DESeq2 QC filter: %d / %d genes retained (count >= 10 in >= 3/6 samples).",
                sum(keep), length(keep)))

vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)), scale. = FALSE)
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2])
pca_df <- data.frame(pca$x[, 1:2], sample = colnames(vsd), condition = condition)

qc_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(show.legend = FALSE, size = 3.5, seed = SEED) +
  scale_color_manual(values = setNames(c(COL_CONTROL, COL_DISEASE),
                                       c(COND_CONTROL, COND_DISEASE))) +
  labs(x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       color = "Condition",
       title = "rna.pub rebuild -- featureCounts raw-count PCA (QC only)",
       subtitle = sprintf("%d genes (count >= 10 in >= 3/6 samples), VST", sum(keep))) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey30"))

ggsave(file.path(OUTPUT_BUILD, "rna_pub_featurecounts_QC_pca.pdf"), plot = qc_pca,
       width = 6, height = 5, dpi = 300)

# =============================================================================
# STEP 5 -- Save (back up the previous n=2 object once)
# =============================================================================

rna_pub_path    <- here::here("Rdata", "pub", "seurat_rna.pub.Rdata")
rna_pub_backup  <- here::here("Rdata", "pub", "seurat_rna.pub.legacy_n2.Rdata")

if (file.exists(rna_pub_path) && !file.exists(rna_pub_backup)) {
  file.copy(rna_pub_path, rna_pub_backup)
  message("Backed up previous rna.pub to: ", rna_pub_backup)
}

save(rna.pub, file = rna_pub_path)
message("Saved: ", rna_pub_path)

# =============================================================================
# Methods note
# =============================================================================

deseq2_version_str   <- as.character(utils::packageVersion("DESeq2"))
msigdbr_version_str  <- as.character(utils::packageVersion("msigdbr"))

methods_note <- sprintf(
"## Methods note -- rna.pub rebuild from raw featureCounts

rna.pub (bulk RNA-seq, primary murine chondrocytes, control vs. IL-1beta,
n = 3 independent samples per condition) was rebuilt from 6 raw
featureCounts (v2.0.0; GRCm38.75 annotation, paired-end) outputs -- the same
6 samples underlying the already-checked-in external VST matrix
(\"control vs il1b_GeneDifferential_vst_count.csv\"), quantified here from
raw alignments instead of an externally VST-transformed matrix. Counts were
merged on %d shared Ensembl gene IDs across all 6 samples, then mapped to
mouse gene symbols using msigdbr's mouse-native (db_species = \"MM\") gene
tables pooled across every MSigDB collection (msigdbr v%s) -- used offline
in place of org.Mm.eg.db/biomaRt, which were unavailable in this
environment. %d / %d genes had a curated symbol (%d dropped, no match);
%d gene symbols backed by more than one Ensembl ID had their counts summed.
%d legacy GEHTS-panel symbol(s) with a newer MGI-approved name in msigdbr
(Hprt1 -> Hprt) were aliased back to the panel's original symbol so the
matrix stays directly comparable to the chip panel (GENES_ALL, config.R)
without editing the panel definition itself.
The resulting %d-gene x 6-sample count matrix is the rna.pub Seurat object's
RNA assay (raw counts + Seurat::NormalizeData() log-normalized data layer);
cells are named control/control.1/control.2/inflammatory/inflammatory.1/
inflammatory.2, matching the cell-naming convention already used by
rna.ref. As a QC check only (not used downstream), DESeq2 (v%s) vst() was
applied to genes with count >= 10 in >= 3/6 samples (%d / %d genes) and
visualized by PCA.
",
  nrow(count_mat), msigdbr_version_str,
  n_genes_mapped, n_genes_in, n_genes_dropped,
  length(dup_symbols),
  length(alias_hits),
  nrow(count_final),
  deseq2_version_str, sum(keep), length(keep)
)
writeLines(methods_note, file.path(OUTPUT_BUILD, "rna_pub_featurecounts_methods_note.md"))
message("Methods note written to: ", file.path(OUTPUT_BUILD, "rna_pub_featurecounts_methods_note.md"))

save_session_info("build_rna_pub_from_featurecounts")
