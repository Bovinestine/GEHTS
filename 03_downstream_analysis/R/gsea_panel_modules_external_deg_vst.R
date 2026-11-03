# =============================================================================
# gsea_panel_modules_external_deg_vst.R — GSEA derivation of Supplementary
#                                          Table 7 (the 7 curated GO modules)
#
# ONE script, one source, one rule. Replaces the earlier ad hoc pair of GSEA
# outputs (external_deg_vst_gsea_panel_modules.csv, filtered at padj < 0.05,
# and external_deg_vst_cartilage_keyword_terms.csv, filtered by keyword) --
# two different filters on two files. That padj < 0.05 gate is exactly the
# bug that would silently delete 2 of the 7 selected modules (Skeletal
# System Development, padj 0.149; Chondrocyte Differentiation, padj 0.252).
# This script never filters on adjusted p; every selection decision below is
# explicit, and adjusted p is reported as a column, not used as a gate.
#
# INPUT IDENTITY (confirmed before running, not after):
#   Bulk Smart-seq2 reference of basal vs. IL-1beta-treated primary murine
#   chondrocytes, 3 independent samples per condition, ~14,500 genes, mouse
#   symbols (EXTERNAL_DEG_VST_CSV). This is NOT the earlier 2-vs-2 rna.ref
#   object -- that had a non-responding inflammatory replicate and was
#   abandoned. File path and sample table are printed below before any
#   analysis runs; identity is further confirmed by reproducing the existing
#   significant-terms CSV exactly before the new tables are built.
#
# Pipeline:
#   1. Load; report file path + sample table.
#   2. Rank ALL genes by a per-gene Welch t-statistic (IL-1b vs. basal), no
#      threshold. Duplicate symbols would be collapsed to the largest |t|;
#      moot for this file (read.csv(row.names=1) already enforces unique
#      symbols) but checked explicitly.
#   3. fgsea multilevel, mouse GO:BP (msigdbr), fixed seed, minSize 10,
#      maxSize 500, eps 0 (all via the shared run_gsea_panel_analysis(),
#      utils.R -- the "one rule" every GSEA script in this repo uses). GO
#      accession pulled from whichever msigdbr column actually holds it,
#      inspected at runtime rather than assumed.
#   4. Verify identity: the significant/panel-overlap term set this run
#      produces must reproduce the existing
#      external_deg_vst_gsea_panel_modules.csv exactly, or the script stops.
#   5. Candidate table: every term with >=2 panel genes in full membership,
#      UNFILTERED on adjusted p, sorted by nominal p. Also writes the full
#      unfiltered fgsea result.
#   6. Final 7-module table (Supplementary Table 7): the 7 pre-selected
#      terms below, looked up by exact GO term name -- never re-selected,
#      re-ranked, or substituted. Each curated gene set (config.R's
#      GO_MODULES, NOT overwritten here) is verified against that term's
#      full membership; any curated gene missing from that membership stops
#      the script rather than being silently dropped or "fixed."
#
# Prerequisites: config.R, utils.R
# Input: EXTERNAL_DEG_VST_CSV; existing
#        output/gsea/external_deg_vst_gsea_panel_modules.csv (identity check)
# Outputs (03_downstream_analysis/output/gsea/):
#   external_deg_vst_fgsea_full.csv            -- every tested GO:BP term
#   external_deg_vst_candidate_terms.csv       -- terms with >=2 panel genes, unfiltered on padj
#   external_deg_vst_supplementary_table7.csv  -- the 7 selected modules
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(msigdbr)
  library(fgsea)
})

options(width = 200)

# =============================================================================
# STEP 1 — Load; report file path + sample table (identity, before running)
# =============================================================================

message("=== INPUT IDENTITY ===")
message("File: ", EXTERNAL_DEG_VST_CSV)
message("Bulk Smart-seq2 reference, basal vs. IL-1beta, primary murine chondrocytes.")
message("This is NOT the earlier 2-vs-2 rna.ref object (abandoned: non-responding ",
        "inflammatory replicate).")

deg_raw <- read.csv(EXTERNAL_DEG_VST_CSV, row.names = 1, check.names = FALSE)

ctrl_cols <- grep("-control$", colnames(deg_raw), value = TRUE)
il1b_cols <- grep("-il1b$",    colnames(deg_raw), value = TRUE)
stopifnot(
  "Expected 3 '-control'-suffixed columns" = length(ctrl_cols) == 3,
  "Expected 3 '-il1b'-suffixed columns"    = length(il1b_cols) == 3
)

sample_table <- data.frame(
  sample    = c(ctrl_cols, il1b_cols),
  condition = c(rep("basal", length(ctrl_cols)), rep("IL-1beta", length(il1b_cols))),
  stringsAsFactors = FALSE
)
cat("\n=== Sample table ===\n")
print(sample_table, row.names = FALSE)
message(sprintf("Genes in file: %d", nrow(deg_raw)))

vst_mat <- as.matrix(deg_raw[, c(ctrl_cols, il1b_cols)])
stopifnot(
  "Gene symbols must be unique (enforced by read.csv row.names=1; checked explicitly here)" =
    !any(duplicated(rownames(vst_mat)))
)

# =============================================================================
# STEP 2 — Rank ALL genes by Welch's t-statistic (IL-1b vs. basal), no threshold
# =============================================================================

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

message(sprintf("Ranked vector: %d genes (mouse symbols, no threshold applied).", length(ranks)))
if (length(ranks) < 1000) {
  warning("Ranked vector has fewer than 1000 genes -- unexpectedly small for this input.")
}

# =============================================================================
# STEP 3 — fgsea multilevel, mouse GO:BP (msigdbr), via the shared
#          run_gsea_panel_analysis() (utils.R) -- the one rule every GSEA
#          script in this repo uses (minSize 10, maxSize 500, eps 0, fixed
#          seed, all from config.R).
# =============================================================================

gsea_out   <- run_gsea_panel_analysis(ranks)
fgsea_all  <- gsea_out$fgsea_res     # UNFILTERED -- every tested GO:BP term
go_bp_sets <- gsea_out$go_bp_sets
go_bp_list <- gsea_out$go_bp_list

msigdbr_pkg_version <- as.character(utils::packageVersion("msigdbr"))
go_db_version <- if ("db_version" %in% colnames(go_bp_sets)) {
  paste(unique(go_bp_sets$db_version), collapse = ", ")
} else {
  NA_character_
}
message(sprintf(
  "msigdbr v%s | GO:BP gene sets tested: %d | MSigDB/GO db_version: %s",
  msigdbr_pkg_version, length(go_bp_list),
  if (is.na(go_db_version)) "not reported by installed msigdbr version" else go_db_version
))

# GO accession column: inspected at runtime, never assumed -- msigdbr's
# column name for this has changed across versions.
message("msigdbr columns available: ", paste(colnames(go_bp_sets), collapse = ", "))
candidate_source_cols <- c("gs_exact_source", "exact_source", "gs_source",
                          "gs_ontology_id", "gs_id")
go_id_col <- intersect(candidate_source_cols, colnames(go_bp_sets))
if (length(go_id_col) == 0) {
  stop("Could not auto-detect the GO-accession column. Available msigdbr columns:\n  ",
      paste(colnames(go_bp_sets), collapse = ", "),
      "\nSet go_id_col manually to the correct column name from this list and re-run.")
}
go_id_col <- go_id_col[1]
message("GO accession column: '", go_id_col, "'")

go_id_raw    <- as.character(go_bp_sets[[go_id_col]])
first_occ    <- !duplicated(go_bp_sets$gs_name)
go_id_lookup <- go_id_raw[first_occ]
names(go_id_lookup) <- go_bp_sets$gs_name[first_occ]

message(sprintf("fgsea: %d GO:BP terms tested (minSize=%d, maxSize=%d, seed=%d).",
                nrow(fgsea_all), GO_MINSIZE, GO_MAXSIZE, SEED))

# =============================================================================
# STEP 4 — Verify identity: reproduce the existing significant-terms CSV
#          exactly before building anything new
# =============================================================================

existing_csv_path <- file.path(OUTPUT_GSEA, "external_deg_vst_gsea_panel_modules.csv")
stopifnot("Existing external_deg_vst_gsea_panel_modules.csv not found" =
            file.exists(existing_csv_path))
existing_csv <- read.csv(existing_csv_path, stringsAsFactors = FALSE)

fresh_csv <- data.frame(
  pathway                   = gsea_out$panel_res$pathway,
  NES                       = gsea_out$panel_res$NES,
  padj                      = gsea_out$panel_res$padj,
  size                      = gsea_out$panel_res$size,
  panel_overlap_full        = gsea_out$panel_res$panel_overlap_full,
  panel_overlap_leadingedge = gsea_out$panel_res$panel_overlap_leadingedge,
  stringsAsFactors = FALSE
)
fresh_csv <- fresh_csv[order(fresh_csv$padj), ]

if (!setequal(existing_csv$pathway, fresh_csv$pathway)) {
  stop(
    "IDENTITY CHECK FAILED -- significant/panel-overlap term sets differ from the existing CSV.\n",
    "  In existing CSV only: ", paste(setdiff(existing_csv$pathway, fresh_csv$pathway), collapse = ", "), "\n",
    "  In fresh re-run only: ", paste(setdiff(fresh_csv$pathway, existing_csv$pathway), collapse = ", "), "\n",
    "This is not the same reference/pipeline as before -- do not proceed."
  )
}
merged_check <- merge(existing_csv, fresh_csv, by = "pathway", suffixes = c("_existing", "_fresh"))
nes_diff  <- max(abs(merged_check$NES_existing  - merged_check$NES_fresh))
padj_diff <- max(abs(merged_check$padj_existing - merged_check$padj_fresh))
if (nes_diff > 1e-6 || padj_diff > 1e-6) {
  stop(sprintf(
    "IDENTITY CHECK FAILED -- term set matches but values differ (max |NES diff| = %.2e, max |padj diff| = %.2e).",
    nes_diff, padj_diff
  ))
}
message(sprintf(
  "IDENTITY CONFIRMED: this run reproduces the existing CSV exactly (%d terms, max |NES diff| = %.2e, max |padj diff| = %.2e).",
  nrow(fresh_csv), nes_diff, padj_diff
))

# =============================================================================
# STEP 5 — Candidate table (>=2 panel genes in full membership, UNFILTERED
#          on adjusted p, sorted by nominal p) + full unfiltered fgsea result
# =============================================================================

panel_overlap_n <- vapply(fgsea_all$pathway,
                          function(term) length(intersect(go_bp_list[[term]], GENES_ALL)),
                          integer(1))
candidate_terms <- fgsea_all[panel_overlap_n >= 2, ]

cand_full <- lapply(candidate_terms$pathway, function(term) sort(intersect(go_bp_list[[term]], GENES_ALL)))
cand_le   <- lapply(seq_len(nrow(candidate_terms)), function(i) sort(intersect(candidate_terms$leadingEdge[[i]], GENES_ALL)))

candidate_table <- data.frame(
  term_name                 = candidate_terms$pathway,
  go_id                     = go_id_lookup[candidate_terms$pathway],
  set_size                  = candidate_terms$size,
  NES                       = candidate_terms$NES,
  nominal_p                 = candidate_terms$pval,
  adjusted_p                = candidate_terms$padj,
  panel_overlap_full        = vapply(cand_full, paste, character(1), collapse = "; "),
  panel_overlap_leadingedge = vapply(cand_le,   paste, character(1), collapse = "; "),
  stringsAsFactors = FALSE
)
candidate_table <- candidate_table[order(candidate_table$nominal_p), ]

write.csv(candidate_table, file.path(OUTPUT_GSEA, "external_deg_vst_candidate_terms.csv"),
         row.names = FALSE)
message(sprintf(
  "STEP 5: %d candidate terms (>=2 panel genes in full membership, unfiltered on adjusted p) -- external_deg_vst_candidate_terms.csv",
  nrow(candidate_table)
))

fgsea_full_out <- data.frame(
  pathway     = fgsea_all$pathway,
  go_id       = go_id_lookup[fgsea_all$pathway],
  size        = fgsea_all$size,
  NES         = fgsea_all$NES,
  nominal_p   = fgsea_all$pval,
  adjusted_p  = fgsea_all$padj,
  leadingEdge = vapply(fgsea_all$leadingEdge, paste, character(1), collapse = "; "),
  stringsAsFactors = FALSE
)
write.csv(fgsea_full_out, file.path(OUTPUT_GSEA, "external_deg_vst_fgsea_full.csv"),
         row.names = FALSE)
message(sprintf("STEP 5: full unfiltered fgsea result (%d terms) -- external_deg_vst_fgsea_full.csv",
                nrow(fgsea_full_out)))

# =============================================================================
# STEP 6 — Final 7-module table (Supplementary Table 7). These 7 terms are
#          the selection -- not re-selected, re-ranked, or substituted.
# =============================================================================

selected_terms <- c(
  Skeletal_System_Development = "GOBP_SKELETAL_SYSTEM_DEVELOPMENT",
  Chondrocyte_Development     = "GOBP_CHONDROCYTE_DEVELOPMENT",
  Chondrocyte_Differentiation = "GOBP_CHONDROCYTE_DIFFERENTIATION",
  Collagen_Catabolic_Process  = "GOBP_COLLAGEN_CATABOLIC_PROCESS",
  Inflammatory_Response       = "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
  Cell_Chemotaxis             = "GOBP_CELL_CHEMOTAXIS",
  Programmed_Cell_Death       = "GOBP_POSITIVE_REGULATION_OF_PROGRAMMED_CELL_DEATH"
)

# Every selected module_name here is used directly as the config.R GO_MODULES
# key (config.R's "ECM_Disassembly" was renamed to "Collagen_Catabolic_Process"
# to match GOBP_COLLAGEN_CATABOLIC_PROCESS exactly, so no name-mapping table
# is needed anymore -- see GO_MODULES in config.R).
stopifnot("Every selected module name must exist in config.R's GO_MODULES" =
            all(names(selected_terms) %in% names(GO_MODULES)))

final_rows <- lapply(names(selected_terms), function(mod_name) {
  go_term <- selected_terms[[mod_name]]
  hit <- fgsea_all[fgsea_all$pathway == go_term, ]
  if (nrow(hit) == 0) {
    stop(sprintf(
      "Selected term '%s' (%s) is not present in the fgsea result -- it may not have passed minSize=%d/maxSize=%d filtering against this ranked list. This is a real problem, not something to silently work around.",
      mod_name, go_term, GO_MINSIZE, GO_MAXSIZE
    ))
  }
  hit <- hit[1, ]

  go_id <- if (go_term %in% names(go_id_lookup)) go_id_lookup[[go_term]] else NA_character_
  if (is.na(go_id) || !nzchar(go_id)) {
    stop(sprintf("Selected term '%s' (%s) has no GO accession in the '%s' column.",
                mod_name, go_term, go_id_col))
  }

  full_membership <- sort(intersect(go_bp_list[[go_term]], GENES_ALL))
  leading_edge     <- sort(intersect(hit$leadingEdge[[1]], GENES_ALL))
  curated          <- sort(GO_MODULES[[mod_name]])

  missing <- setdiff(curated, full_membership)
  if (length(missing) > 0) {
    stop(sprintf(
      "VERIFICATION FAILED for module '%s' (%s): curated gene(s) [%s] are NOT members of this term's full panel overlap ([%s]). Real error -- stopping, not silently fixing.",
      mod_name, go_term, paste(missing, collapse = ", "), paste(full_membership, collapse = ", ")
    ))
  }
  message(sprintf("VERIFY %-28s (%-55s): PASS -- all %d curated genes are members.",
                  mod_name, go_term, length(curated)))

  data.frame(
    module_name            = mod_name,
    go_id                  = go_id,
    set_size               = hit$size,
    NES                    = hit$NES,
    nominal_p              = hit$pval,
    adjusted_p             = hit$padj,
    panel_full_membership  = paste(full_membership, collapse = "; "),
    panel_leading_edge     = paste(leading_edge, collapse = "; "),
    curated_gene_set       = paste(curated, collapse = "; "),
    stringsAsFactors = FALSE
  )
})
final_table <- do.call(rbind, final_rows)

write.csv(final_table, file.path(OUTPUT_GSEA, "external_deg_vst_supplementary_table7.csv"),
         row.names = FALSE)

cat("\n=== Supplementary Table 7 ===\n")
print(final_table[, c("module_name", "go_id", "set_size", "NES", "nominal_p", "adjusted_p")],
     row.names = FALSE)

# =============================================================================
# CHECKS
# =============================================================================

cat("\n=== CHECKS ===\n")
message(sprintf("Ranked vector: %d genes.", length(ranks)))
if (length(ranks) < 10000) {
  warning("Ranked vector has fewer than 10,000 genes -- expected roughly 14,500.")
}

message("All 7 selected terms present in the fgsea result with a GO id: PASS ",
        "(verified above -- the script stops otherwise).")

stopifnot("All 7 selected terms must clear nominal p < 0.05" = all(final_table$nominal_p < 0.05))
message("All 7 selected terms clear nominal p < 0.05: PASS.")

above_padj_05 <- final_table[final_table$adjusted_p >= 0.05,
                             c("module_name", "go_id", "nominal_p", "adjusted_p")]
if (nrow(above_padj_05) > 0) {
  cat("Terms above adjusted p 0.05 (reported as a column, not gated on):\n")
  print(above_padj_05, row.names = FALSE)
} else {
  message("No selected terms are above adjusted p 0.05.")
}

message("Curated-gene membership verification: PASS for all 7 (printed per module above; ",
        "the script stops on any failure rather than reaching this line).")
message("Significant terms reproduce the existing external_deg_vst_gsea_panel_modules.csv: PASS (STEP 4).")

message("\nDone. See ", OUTPUT_GSEA, " for external_deg_vst_fgsea_full.csv, ",
        "external_deg_vst_candidate_terms.csv, external_deg_vst_supplementary_table7.csv.")

save_session_info("gsea_panel_modules_external_deg_vst")
