# =============================================================================
# gsea_panel_modules_external_deg_vst.R — fgsea + panel modules from an
#                                          external, near-genome-wide VST
#                                          matrix (FINAL / primary GO:BP
#                                          module-derivation pipeline)
#
# EXTERNAL_DEG_VST_CSV ("control vs il1b_GeneDifferential_vst_count.csv") is
# the DEG data source used for the GO_MODULES ultimately checked into
# config.R -- see README.md "Module definitions". It holds ~14.5k genes
# (near-genome-wide, NOT pre-filtered to a "significant" subset), so the
# GSEA-needs-a-full-background caveat that applies to a pre-filtered subset
# does not apply here.
#
# Pipeline:
#   1. Load the CSV; identify 3 control- and 3 il1b-condition columns; QC PCA.
#   2. Panel gene table: direct VST fold changes for the 30-gene panel.
#   3. Rank all ~14.5k genes by Welch's t-statistic (il1b vs. control).
#   4. fgsea multilevel on mouse GO:BP (msigdbr) -- identical settings
#      (GO_MINSIZE, GO_MAXSIZE, PADJ_CUTOFF, SEED) across every GSEA variant
#      in this repo, for comparability.
#   5. Panel modules: leading edge intersected with GENES_ALL.
#   6. Methods note.
#   7. Diagnostic (previously check_cartilage_suppression_vst.R): is the
#      anabolic/cartilage program actually suppressed in this reference, or
#      is the single significant anabolic term (chondrocyte development) an
#      FDR-threshold artifact? Reuses the Step 3/4 ranking and fgsea result
#      directly -- no re-derivation, so there is nothing to reconcile against
#      a separately-run CSV.
#
# Prerequisites: config.R, utils.R
# Input: EXTERNAL_DEG_VST_CSV (Rdata/pub/"control vs il1b_GeneDifferential_
#        vst_count.csv")
# Outputs (03_downstream_analysis/output/gsea/):
#   external_deg_vst_QC_pca.pdf
#   external_deg_vst_panel_gene_foldchanges.csv
#   external_deg_vst_gsea_panel_modules.csv
#   GO_MODULES_FROM_EXTERNAL_DEG_VST.R
#   external_deg_vst_methods_note.md
#   external_deg_vst_cartilage_keyword_terms.csv
#   external_deg_vst_top30_negative_terms.csv
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
  library(ggrepel)
})

options(width = 200)

# =============================================================================
# STEP 1 — Load + QC
# =============================================================================

deg_raw <- read.csv(EXTERNAL_DEG_VST_CSV, row.names = 1, check.names = FALSE)

ctrl_cols <- grep("-control$", colnames(deg_raw), value = TRUE)
il1b_cols <- grep("-il1b$",    colnames(deg_raw), value = TRUE)
stopifnot(
  "Expected 3 '-control'-suffixed columns" = length(ctrl_cols) == 3,
  "Expected 3 '-il1b'-suffixed columns"    = length(il1b_cols) == 3
)
message(sprintf(
  "External VST CSV: %d genes (near-genome-wide). control = %s | il1b = %s",
  nrow(deg_raw), paste(ctrl_cols, collapse = ", "), paste(il1b_cols, collapse = ", ")
))

vst_mat <- as.matrix(deg_raw[, c(ctrl_cols, il1b_cols)])

pca <- prcomp(t(vst_mat), scale. = FALSE)
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2])
pca_df <- data.frame(pca$x[, 1:2], sample = rownames(pca$x),
                     condition = ifelse(rownames(pca$x) %in% ctrl_cols,
                                        COND_CONTROL, COND_DISEASE))

qc_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(show.legend = FALSE, size = 3.5, seed = SEED) +
  scale_color_manual(values = setNames(c(COL_CONTROL, COL_DISEASE),
                                       c(COND_CONTROL, COND_DISEASE))) +
  labs(x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       color = "Condition",
       title = "External VST CSV — PCA (QC only)",
       subtitle = sprintf("%d genes, near-genome-wide", nrow(vst_mat))) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey30"))

ggsave(file.path(OUTPUT_GSEA, "external_deg_vst_QC_pca.pdf"), plot = qc_pca,
       width = 6, height = 5, dpi = 300)

# =============================================================================
# STEP 2 — Panel gene fold changes (direct, on VST scale) — 30-gene panel check
# =============================================================================

panel_present <- intersect(GENES_ALL, rownames(vst_mat))
panel_missing <- setdiff(GENES_ALL, rownames(vst_mat))
if (length(panel_missing) > 0) {
  message("Panel genes not in this VST file: ", paste(panel_missing, collapse = ", "))
}

category <- rep(NA_character_, length(panel_present))
category[panel_present %in% GENES_ANABOLIC]     <- "Anabolic"
category[panel_present %in% GENES_CATABOLIC]    <- "Catabolic"
category[panel_present %in% GENES_HOUSEKEEPING] <- "Housekeeping"

panel_vst <- vst_mat[panel_present, , drop = FALSE]
panel_fc <- data.frame(
  gene          = panel_present,
  category      = category,
  ctrl_mean_vst = rowMeans(panel_vst[, ctrl_cols, drop = FALSE]),
  il1b_mean_vst = rowMeans(panel_vst[, il1b_cols, drop = FALSE]),
  stringsAsFactors = FALSE
)
panel_fc$log2fc_vst <- panel_fc$il1b_mean_vst - panel_fc$ctrl_mean_vst
panel_fc <- cbind(panel_fc, as.data.frame(panel_vst))
panel_fc <- panel_fc[order(panel_fc$category, -panel_fc$log2fc_vst), ]

write.csv(panel_fc, file.path(OUTPUT_GSEA, "external_deg_vst_panel_gene_foldchanges.csv"),
         row.names = FALSE)
message(sprintf(
  "Panel genes present in the VST file: %d / %d -- external_deg_vst_panel_gene_foldchanges.csv",
  length(panel_present), length(GENES_ALL)
))

# =============================================================================
# STEP 3 — Rank — Welch's two-sample t-statistic (il1b vs. control), n = 3/group
# =============================================================================
# n = 3 per group supports a genuine per-gene variance estimate, unlike
# rna.ref's n = 2 (which is why DESeq2 relies on empirical-Bayes dispersion
# shrinkage there). Computed manually (not via t.test() per gene) so genes
# with zero within-group variance produce Inf/NaN rather than an error --
# those are dropped explicitly below instead of crashing the loop.

mean_ctrl <- rowMeans(vst_mat[, ctrl_cols, drop = FALSE])
mean_il1b <- rowMeans(vst_mat[, il1b_cols, drop = FALSE])
var_ctrl  <- apply(vst_mat[, ctrl_cols, drop = FALSE], 1, var)
var_il1b  <- apply(vst_mat[, il1b_cols, drop = FALSE], 1, var)
n_ctrl <- length(ctrl_cols)
n_il1b <- length(il1b_cols)

se     <- sqrt(var_il1b / n_il1b + var_ctrl / n_ctrl)
t_stat <- (mean_il1b - mean_ctrl) / se
names(t_stat) <- rownames(vst_mat)

n_nonfinite <- sum(!is.finite(t_stat))
if (n_nonfinite > 0) {
  message(sprintf(
    "Dropped %d / %d genes with a non-finite Welch t-statistic (zero replicate variance in one or both groups).",
    n_nonfinite, length(t_stat)
  ))
}
ranks <- sort(t_stat[is.finite(t_stat)], decreasing = TRUE)

message(sprintf("Ranked vector: %d genes (near-genome-wide background).", length(ranks)))
if (length(ranks) < 1000) {
  warning("Ranked vector has fewer than 1000 genes -- unexpectedly small for this input.")
}

# =============================================================================
# STEP 4 — fgsea multilevel on mouse GO:BP + panel modules from leading edge
# =============================================================================

gsea_out <- run_gsea_panel_analysis(ranks)

go_bp_list        <- gsea_out$go_bp_list
fgsea_all         <- gsea_out$fgsea_res    # UNFILTERED -- every tested GO:BP term (Step 7)
n_pathways_tested <- gsea_out$n_pathways_tested
n_significant     <- gsea_out$n_significant
n_panel_overlap   <- gsea_out$n_panel_overlap
panel_res         <- gsea_out$panel_res
GO_MODULES_FROM_EXTERNAL_DEG_VST <- gsea_out$panel_modules

message(sprintf("fgsea: %d / %d GO:BP terms reached minSize=%d; %d significant (BH padj < %.2f).",
                n_pathways_tested, length(go_bp_list), GO_MINSIZE, n_significant, PADJ_CUTOFF))
message(sprintf("Panel-overlap filter: %d / %d significant terms retained.",
                n_panel_overlap, n_significant))

csv_out <- data.frame(
  pathway                   = panel_res$pathway,
  NES                       = panel_res$NES,
  padj                      = panel_res$padj,
  size                      = panel_res$size,
  panel_overlap_full        = panel_res$panel_overlap_full,
  panel_overlap_leadingedge = panel_res$panel_overlap_leadingedge,
  stringsAsFactors = FALSE
)
csv_out <- csv_out[order(csv_out$padj), ]
write.csv(csv_out, file.path(OUTPUT_GSEA, "external_deg_vst_gsea_panel_modules.csv"),
         row.names = FALSE)

write_go_modules_file(
  GO_MODULES_FROM_EXTERNAL_DEG_VST,
  file.path(OUTPUT_GSEA, "GO_MODULES_FROM_EXTERNAL_DEG_VST.R"),
  list_name = "GO_MODULES_FROM_EXTERNAL_DEG_VST",
  header_lines = c(
    "# Source: external, near-genome-wide VST matrix (NOT rna.ref, NOT DESeq2 --",
    "# DESeq2 cannot be re-run on already-VST-transformed values). Ranked by a",
    "# per-gene Welch's t-statistic (il1b vs. control, n=3/group) -> fgsea",
    "# multilevel on mouse GO:BP (msigdbr). Each entry is one significant term",
    paste0("# (BH padj < ", PADJ_CUTOFF, ") whose leading-edge genes intersect the"),
    "# 30-gene panel (GENES_ALL, config.R), formatted like GO_MODULES in config.R.",
    "# No biological-relevance filter was applied. This is the primary,",
    "# fully-used GSEA input (see config.R's GO_MODULES / README.md).",
    paste0("# Significant terms (padj < ", PADJ_CUTOFF, "): ", n_significant,
          " | Retained after panel-overlap filter: ", n_panel_overlap)
  )
)

# =============================================================================
# STEP 5 — Methods note
# =============================================================================

fgsea_version_str   <- as.character(utils::packageVersion("fgsea"))
msigdbr_version_str <- as.character(utils::packageVersion("msigdbr"))

methods_note <- sprintf(
"## Methods note — GSEA on an external, near-genome-wide VST matrix

An external collaborator's VST-transformed count matrix (%s; %d genes, near-genome-wide -- not pre-filtered to a significant subset) was used to rank genes for GSEA. Because this file holds already variance-stabilized values rather than raw counts, DESeq2 could not be re-fit (DESeqDataSetFromMatrix()/DESeq() require raw integer counts to estimate dispersions); a per-gene Welch's two-sample t-statistic (il1b vs. control, n = 3 independent samples per condition) was computed directly on the VST values instead and used to rank %d genes with a finite statistic (%d genes with zero within-group variance in one or both conditions were excluded). fgsea::fgseaMultilevel() (fgsea v%s) was run against mouse Gene Ontology Biological Process gene sets (msigdbr v%s, collection C5, subcollection GO:BP) with minSize = %d, maxSize = %d, eps = 0, and a fixed seed (%d). %d / %d GO:BP terms reached the minSize overlap threshold; %d were significant (BH padj < %.2f), of which %d had leading-edge overlap with the 30-gene GEHTS panel (GO_MODULES_FROM_EXTERNAL_DEG_VST). This ranking covers a near-genome-wide background, so the standard GSEA full-background assumption is reasonably satisfied. This is the sole DEG data source used for the GO_MODULES ultimately checked into config.R.
",
  basename(EXTERNAL_DEG_VST_CSV), nrow(deg_raw), length(ranks), n_nonfinite,
  fgsea_version_str, msigdbr_version_str,
  GO_MINSIZE, GO_MAXSIZE, SEED,
  n_pathways_tested, length(go_bp_list), n_significant, PADJ_CUTOFF, n_panel_overlap
)
writeLines(methods_note, file.path(OUTPUT_GSEA, "external_deg_vst_methods_note.md"))

message("\nSteps 1-5 done. See ", OUTPUT_GSEA, " for external_deg_vst_* outputs.")

# =============================================================================
# STEP 6 — Diagnostic: is the anabolic/cartilage program actually suppressed
#          in this reference, or is the significant table's single anabolic
#          term (chondrocyte development) an FDR-threshold artifact hiding a
#          real-but-weak cartilage-suppression signal (or the absence of one)?
#
# DIAGNOSTIC ONLY. Does not emit a module definition file. Reuses `ranks`,
# `vst_mat`, `mean_ctrl`/`mean_il1b`/`var_ctrl`/`var_il1b`, `fgsea_all`, and
# `go_bp_list` from Steps 3-4 directly -- nothing here is re-derived, so there
# is no separate run to reconcile against.
# =============================================================================

# -----------------------------------------------------------------------------
# 6a — Gene level first: where do the anabolic panel genes sit in the ranked
#      list? (load-bearing check)
# -----------------------------------------------------------------------------
# Rank 1 = most UP in inflammatory (highest t-statistic); rank = total_ranked
# = most DOWN in inflammatory (most suppressed). A gene suppressed by
# inflammation should sit near total_ranked, not near the middle.

anabolic_check_genes <- c("Acan", "Col2a1", "Sox9", "Matn1", "Matn3", "Ucma",
                          "Col27a1", "Pth1r", "Ccnd3", "Gadd45g", "Gm26633")
stopifnot("anabolic_check_genes must match GENES_ANABOLIC (config.R) as a set" =
            setequal(anabolic_check_genes, GENES_ANABOLIC))

# Per-gene log2FC (VST-scale mean difference) and a Welch p-value/padj added
# ON TOP OF the existing ranking t-statistic for reporting only -- this does
# NOT alter the ranking metric, gene sets, seed, or fgsea settings used above.
log2fc_all <- mean_il1b - mean_ctrl
df_welch <- (var_il1b / n_il1b + var_ctrl / n_ctrl)^2 /
  ((var_il1b / n_il1b)^2 / (n_il1b - 1) + (var_ctrl / n_ctrl)^2 / (n_ctrl - 1))
pval_all <- 2 * pt(-abs(t_stat), df = df_welch)

de_table <- data.frame(
  gene       = names(t_stat),
  t_stat     = as.numeric(t_stat),
  log2FC_vst = as.numeric(log2fc_all[names(t_stat)]),
  pvalue     = as.numeric(pval_all),
  stringsAsFactors = FALSE
)
de_table <- de_table[is.finite(de_table$t_stat), ]     # matches `ranks` exactly
de_table$padj <- p.adjust(de_table$pvalue, method = "BH")

rank_position <- setNames(seq_along(ranks), names(ranks))
de_table$rank_position <- rank_position[de_table$gene]
de_table$total_ranked  <- length(ranks)

present_genes <- intersect(anabolic_check_genes, de_table$gene)
missing_genes <- setdiff(anabolic_check_genes, de_table$gene)
if (length(missing_genes) > 0) {
  absent_from_file <- !(missing_genes %in% rownames(vst_mat))
  message("Anabolic panel genes NOT in the ranked list: ",
          paste(sprintf("%s (%s)", missing_genes,
                        ifelse(absent_from_file, "absent from VST file", "excluded, non-finite t-stat")),
                collapse = ", "))
}

gene_table <- de_table[match(present_genes, de_table$gene),
                       c("gene", "t_stat", "log2FC_vst", "pvalue", "padj",
                         "rank_position", "total_ranked")]
gene_table[c("t_stat", "log2FC_vst")] <- lapply(gene_table[c("t_stat", "log2FC_vst")], round, 4)
gene_table[c("pvalue", "padj")]       <- lapply(gene_table[c("pvalue", "padj")], signif, 4)

cat("\n=== STEP 6a: Anabolic panel gene positions in the ranked list ===\n")
cat("(rank_position 1 = most UP in inflammatory; rank_position = total_ranked = most DOWN)\n")
print(gene_table, row.names = FALSE)

# -----------------------------------------------------------------------------
# 6b — Term level, unfiltered and unbiased by direction: every GO:BP term
#      matching cartilage/matrix/skeletal keywords, no padj cutoff
# -----------------------------------------------------------------------------

keywords <- c("CARTILAGE", "CHONDRO", "COLLAGEN", "EXTRACELLULAR_MATRIX",
             "EXTERNAL_ENCAPSULATING", "PROTEOGLYCAN", "GLYCOSAMINOGLYCAN",
             "SKELET", "OSSIFICATION", "CONNECTIVE_TISSUE", "GROWTH_PLATE",
             "TISSUE_HOMEOSTASIS")
keyword_pattern <- paste(keywords, collapse = "|")

keyword_hits <- fgsea_all[grepl(keyword_pattern, fgsea_all$pathway, ignore.case = TRUE), ]

# fgsea's own leadingEdge already reflects genes from the BOTTOM of the ranked
# list for negative-NES terms (and the top for positive-NES terms) -- used
# directly here, never re-derived from `ranks`, so direction is never mislabeled.
kw_overlap_full <- lapply(keyword_hits$pathway, function(term) sort(intersect(go_bp_list[[term]], GENES_ALL)))
kw_overlap_le   <- lapply(seq_len(nrow(keyword_hits)), function(i) sort(intersect(keyword_hits$leadingEdge[[i]], GENES_ALL)))

keyword_hits$panel_overlap_full        <- vapply(kw_overlap_full, paste, character(1), collapse = "; ")
keyword_hits$panel_overlap_leadingedge <- vapply(kw_overlap_le,   paste, character(1), collapse = "; ")

keyword_csv <- keyword_hits[, c("pathway", "NES", "pval", "padj", "size",
                                "panel_overlap_full", "panel_overlap_leadingedge")]
keyword_csv <- keyword_csv[order(keyword_csv$NES), ]   # ascending, positive-NES terms kept

write.csv(keyword_csv, file.path(OUTPUT_GSEA, "external_deg_vst_cartilage_keyword_terms.csv"),
         row.names = FALSE)
message(sprintf("STEP 6b: %d cartilage/matrix/skeletal-keyword terms found (%d negative NES, %d positive/zero NES) -- external_deg_vst_cartilage_keyword_terms.csv",
                nrow(keyword_csv), sum(keyword_csv$NES < 0), sum(keyword_csv$NES >= 0)))

# -----------------------------------------------------------------------------
# 6c — Unbiased view: the 30 most negatively enriched terms overall, no
#      keyword filter, no significance filter
# -----------------------------------------------------------------------------

neg30 <- fgsea_all[order(fgsea_all$NES), ]
neg30 <- head(neg30, 30)

neg30_full <- lapply(neg30$pathway, function(term) sort(intersect(go_bp_list[[term]], GENES_ALL)))
neg30_le   <- lapply(seq_len(nrow(neg30)), function(i) sort(intersect(neg30$leadingEdge[[i]], GENES_ALL)))

neg30$panel_overlap_full        <- vapply(neg30_full, paste, character(1), collapse = "; ")
neg30$panel_overlap_leadingedge <- vapply(neg30_le,   paste, character(1), collapse = "; ")

neg30_csv <- neg30[, c("pathway", "NES", "pval", "padj", "size",
                       "panel_overlap_full", "panel_overlap_leadingedge")]
write.csv(neg30_csv, file.path(OUTPUT_GSEA, "external_deg_vst_top30_negative_terms.csv"),
         row.names = FALSE)
message("STEP 6c: top-30 most negative terms written -- external_deg_vst_top30_negative_terms.csv")

# 6d — Leading-edge direction sanity note: confirmed by construction in 6b-6c
# above. keyword_hits$leadingEdge and neg30$leadingEdge are taken directly
# from fgsea::fgseaMultilevel()'s own output, which computes the leading edge
# relative to each term's own enrichment direction (bottom of the ranked list
# for negative NES, top for positive NES). No leading-edge genes were
# re-derived from `ranks` in this script, so there is no risk of mislabeling
# direction.

# -----------------------------------------------------------------------------
# 6e — Decision rule and verdict (does not choose modules)
# -----------------------------------------------------------------------------

kw_neg <- keyword_hits[keyword_hits$NES < 0, ]
kw_neg$panel_overlap_full_n <- ifelse(kw_neg$panel_overlap_full == "", 0,
                                      lengths(strsplit(kw_neg$panel_overlap_full, "; ")))

# Terms already padj < 0.05 (e.g. the known "chondrocyte development" hit)
# are already in external_deg_vst_gsea_panel_modules.csv -- surfaced
# separately here so they are never mislabeled "not defensible" by Branch B.
already_sig   <- kw_neg[!is.na(kw_neg$padj) & kw_neg$padj < 0.05, ]
kw_neg_notsig <- kw_neg[is.na(kw_neg$padj) | kw_neg$padj >= 0.05, ]

branch_a <- kw_neg_notsig[!is.na(kw_neg_notsig$padj) & kw_neg_notsig$padj <= 0.25 &
                          kw_neg_notsig$panel_overlap_full_n >= 2, ]

cat("\n=== STEP 6e: Verdict ===\n")
cat(sprintf("Cartilage/matrix/skeletal keyword terms tested: %d (%d negative NES, %d positive/zero NES).\n",
           nrow(keyword_hits), nrow(kw_neg), nrow(keyword_hits) - nrow(kw_neg)))
if (nrow(already_sig) > 0) {
  cat(sprintf("Already FDR-significant (padj < 0.05, negative NES) -- already in external_deg_vst_gsea_panel_modules.csv, not a new finding here (%d term[s]):\n",
             nrow(already_sig)))
  print(already_sig[order(already_sig$padj), c("pathway", "NES", "padj", "panel_overlap_full_n")], row.names = FALSE)
}

if (nrow(keyword_hits) == 0) {
  cat("VERDICT: no cartilage/matrix/skeletal keyword terms were tested at all against this ranked list.\n",
      "Cannot assess suppression from GO:BP terms here -- rely on the Step 6a gene-level table.\n", sep = "")
} else if (nrow(branch_a) > 0) {
  cat("VERDICT: BRANCH A -- candidate cartilage/matrix term(s) beyond the already-significant one(s)\n",
      "above trended down but did not clear FDR. Reportable as exploratory candidates with padj\n",
      "stated openly:\n", sep = "")
  print(branch_a[order(branch_a$padj), c("pathway", "NES", "padj", "panel_overlap_full_n", "panel_overlap_full")],
       row.names = FALSE)
} else if (nrow(kw_neg_notsig) > 0) {
  cat("VERDICT: BRANCH B -- beyond the already-significant term(s) above, remaining cartilage/matrix\n",
      "terms trend negative but none clear the candidate bar (padj in [0.05, 0.25] AND >=2 panel\n",
      "genes in full overlap). Not defensible as an additional data-derived module. Recommend\n",
      "keeping curated anabolic modules (GO_MODULES in config.R) for anything beyond what is\n",
      "already significant.\n", sep = "")
  cat("Negative-NES, non-significant keyword terms found (none qualify for Branch A):\n")
  print(kw_neg_notsig[order(kw_neg_notsig$padj), c("pathway", "NES", "padj", "panel_overlap_full_n")],
       row.names = FALSE)
} else if (nrow(already_sig) == 0) {
  cat("VERDICT: BRANCH C -- every cartilage/matrix/skeletal keyword term tested has positive or\n",
      "zero NES. This reference does NOT suppress the cartilage/anabolic program at the GO:BP\n",
      "term level. Curated anabolic modules are the only honest option for this reference, and\n",
      "this asymmetry (catabolic/inflammatory terms suppress cleanly, anabolic/cartilage terms\n",
      "do not) must be stated explicitly in the methods.\n", sep = "")
  cat("Positive/zero-NES keyword terms found:\n")
  print(keyword_hits[order(-keyword_hits$NES), c("pathway", "NES", "padj")], row.names = FALSE)
} else {
  cat("VERDICT: every negative-NES cartilage/matrix keyword term is already FDR-significant (see\n",
      "above) -- no additional borderline (padj 0.05-0.25) candidates exist beyond what is already\n",
      "captured in external_deg_vst_gsea_panel_modules.csv.\n", sep = "")
}

message("\nStep 6 diagnostic complete -- no module definition file written by this step.")

save_session_info("gsea_panel_modules_external_deg_vst")
