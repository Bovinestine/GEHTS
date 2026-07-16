# =============================================================================
# export_combo_test_data.R — Export SCT residuals + drug labels to CSV for
# figure5g_test.ipynb (Python combination-scoring test on the 7-module set)
#
# figure5g.ipynb's production inputs (mac_cmb0.1.sct.h5ad, mac_sin0.1.sct.h5ad)
# are not present anywhere in this repo -- built by an external, undocumented
# step (per README.md, the Python notebooks "must be run separately"). This
# script is a self-contained substitute: it exports the SAME cell populations
# (mac basal/IL-1b + cmb combos at 0.1 uM; mac basal/IL-1b + sin singles at
# 0.1 uM) directly from the .Rdata objects the R test harness already uses,
# as CSV, so the Python test notebook has real data to run against.
#
# Prerequisites: config.R
# Outputs (03_downstream_analysis/output/test/):
#   mac_cmb0.1_scores.csv -- mac (basal/IL-1b) + cmb (0.1 uM) cells. Combo
#     drug_name is the SORTED "min&max" form (pmin/pmax of drug_name1/2) --
#     matching test_new_modules.R's std_combi convention exactly, NOT the
#     well's raw drug1/drug2 order (which is fixed per combo but arbitrary
#     across combos, and does not match the sorted names in top5_combos.json).
#   mac_sin0.1_scores.csv -- mac (basal/IL-1b) + sin (0.1 uM) cells
#   top5_combos.json -- the same Bliss-synergy top-5 combo ranking
#     test_new_modules.R computes for the Fig 5f dot plot, so the Python
#     notebook can select one of the SAME combos for an apples-to-apples
#     comparison instead of picking an arbitrary, unrelated one.
#   Columns: cell, drug_name, <30 GENES_ALL panel gene residuals>
# =============================================================================

source(here::here("03_downstream_analysis", "R", "config.R"))
source(here::here("03_downstream_analysis", "R", "utils.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

.load_rdata <- function(path, description) {
  tryCatch(load(path, envir = .GlobalEnv),
           error = function(e) stop("Failed to load ", description, "\n  Path: ", path,
                                     "\n  Error: ", conditionMessage(e)))
  message("Loaded: ", description)
}
.load_rdata(file.path(DATA_DIR, "seurat_mac.sct.Rdata"), "mac.sct - primary MAC (SCT)")
.load_rdata(file.path(DATA_DIR, "seurat_sin.sct.Rdata"), "sin.sct - single-drug screen (SCT)")
.load_rdata(file.path(DATA_DIR, "seurat_cmb.sct.Rdata"), "cmb.sct - combo-drug screen (SCT)")
.load_rdata(file.path(PRED_DIR, "prob_cmb0.1_rf.Rdata"), "pred_combo_0.1 - RF combo predictions, 0.1 uM")
.load_rdata(file.path(PRED_DIR, "prob_sng0.1_rf.Rdata"), "pred_single_0.1 - RF single predictions, 0.1 uM")
prob_cmb0.1_rf <- pred_combo_0.1
prob_sng0.1_rf <- pred_single_0.1


mac.sct$drug_name <- mac.sct$drug_condition

export_scores <- function(seurat_obj, drug_name_vec, out_path) {
  mat <- GetAssayData(seurat_obj, assay = "SCT", layer = "scale.data")
  df  <- as.data.frame(t(as.matrix(mat[GENES_ALL, , drop = FALSE])))
  df$cell      <- rownames(df)
  df$drug_name <- drug_name_vec
  df <- df[, c("cell", "drug_name", GENES_ALL)]
  write.csv(df, out_path, row.names = FALSE)
  message("Wrote: ", out_path, " (", nrow(df), " cells)")
}

# --- mac (basal/IL-1b) + cmb (0.1 uM combos) ---
cmb0.1.sct <- subset(cmb.sct, subset = dose1 == 0.1)
# Sorted "min&max" form -- matches test_new_modules.R's std_combi exactly,
# so drug_name here is directly comparable/joinable to top5_combos.json.
cmb0.1.sct$drug_name <- paste0(
  pmin(as.character(cmb0.1.sct$drug_name1), as.character(cmb0.1.sct$drug_name2)), "&",
  pmax(as.character(cmb0.1.sct$drug_name1), as.character(cmb0.1.sct$drug_name2))
)
mac_cmb0.1.sct <- merge(mac.sct, cmb0.1.sct, add.cell.ids = c("mac", "cmb"))
export_scores(mac_cmb0.1.sct, unname(mac_cmb0.1.sct$drug_name),
             file.path(OUTPUT_TEST, "mac_cmb0.1_scores.csv"))

# --- mac (basal/IL-1b) + sin (0.1 uM singles) ---
sin0.1.sct <- subset(sin.sct, subset = dose == 0.1)
mac_sin0.1.sct <- merge(mac.sct, sin0.1.sct, add.cell.ids = c("mac", "sin"))
export_scores(mac_sin0.1.sct, unname(as.character(mac_sin0.1.sct$drug_name)),
             file.path(OUTPUT_TEST, "mac_sin0.1_scores.csv"))

# --- Top-5 Bliss synergy combos -- identical ranking to test_new_modules.R,
#     so the notebook can pick from the SAME combos Fig 5f actually plots ---
eff_obs_mean <- vapply(prob_cmb0.1_rf, mean, numeric(1))
splits <- strsplit(names(eff_obs_mean), "&")
drug1  <- trimws(vapply(splits, `[`, character(1), 1))
drug2  <- trimws(vapply(splits, `[`, character(1), 2))
sng_mean <- vapply(prob_sng0.1_rf, mean, numeric(1))
bliss_df <- data.frame(
  std_combi = paste0(pmin(drug1, drug2), "&", pmax(drug1, drug2)),
  drug1 = drug1, drug2 = drug2,
  E_obs = as.numeric(eff_obs_mean),
  Synergy = as.numeric(eff_obs_mean) - bliss_predict(sng_mean[drug1], sng_mean[drug2])
) %>%
  filter(!is.na(Synergy), drug1 != drug2) %>%
  group_by(std_combi) %>%
  summarise(E_obs = mean(E_obs), Synergy = mean(Synergy), .groups = "drop") %>%
  mutate(Rank_Sum = rank(-E_obs, ties.method = "min") + rank(-Synergy, ties.method = "min")) %>%
  arrange(Rank_Sum, desc(E_obs))
top5_combos <- head(bliss_df$std_combi, 5)
write_json(top5_combos, file.path(OUTPUT_TEST, "top5_combos.json"))
message("Wrote: ", file.path(OUTPUT_TEST, "top5_combos.json"), " -- ",
        paste(top5_combos, collapse = "; "))
