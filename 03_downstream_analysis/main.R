# =============================================================================
# main.R — Master entry point for all GEHTS downstream analyses
#
# Usage:
#   From the project root:
#     Rscript 03_downstream_analysis/main.R
#   Or interactively:
#     source(here::here("03_downstream_analysis", "main.R"))
#
# Order of execution:
#   1. config.R            — global constants, gene panels, paths, palette
#   2. utils.R             — shared plotting helpers
#   3. data_loading.R      — load all Seurat / RData objects
#   4. scoring_methods.R   — train ML models, save CV metrics to PRED_DIR/
#   5. prediction_score.R  — generate efficacy predictions
#   6. figure1.R           — Figure 1 (platform comparison)
#   7. figure2.R           — Figure 2 (heatmap, UMAP, volcano, divergence)
#   8. figure2_supplement.R — Ext. Data Fig. 5a, 5b, 6, 7
#   9. figure3.R           — Figure 3b–3c (heatmap, UMAP)
#  10. figure3_supplement.R — Ext. Data Fig. 3 (boxplots), Ext. Data Fig. 8 (UMAPs)
#  11. figure4.R           — Figure 4b–4f
#  12. figure4_supplement.R — Ext. Data (CV comparison, GO violins)
#  13. figure5.R           — Figure 5c–5f (Bliss, synergy scatter, dot plot)
#
# Python notebooks (run separately — see README.md):
#   figures/figure5_supplement.ipynb  — Extended Figure 5 (class-sum barplot)
#   figures/figure5g.ipynb            — Figure 5g (per-module synergy)
#
# Outputs:
#   03_downstream_analysis/output/main/          ← main figure PDFs
#   03_downstream_analysis/output/supplementary/ ← supplementary PDFs
#   03_downstream_analysis/session_info/          ← sessionInfo() per script
#
# Expected runtime: ~30–60 min on a machine with ≥16 GB RAM
# Tested on R 4.3.x; see session_info/ for exact package versions
# =============================================================================

suppressPackageStartupMessages(library(here))

message("=== GEHTS downstream analysis — starting ===")
message("Working directory: ", getwd())
message("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# ---------------------------------------------------------------------------
# Core infrastructure
# ---------------------------------------------------------------------------
source(here("03_downstream_analysis", "R", "config.R"))
source(here("03_downstream_analysis", "R", "utils.R"))
source(here("03_downstream_analysis", "R", "data_loading.R"))

# ---------------------------------------------------------------------------
# ML model training and prediction
# ---------------------------------------------------------------------------
source(here("03_downstream_analysis", "R", "scoring_methods.R"))
source(here("03_downstream_analysis", "R", "prediction_score.R"))

# ---------------------------------------------------------------------------
# Figure scripts (run in figure order)
# ---------------------------------------------------------------------------
figure_scripts <- c(
  "figure1.R",
  "figure2.R",
  "figure2_supplement.R",
  "figure3.R",
  "figure3_supplement.R",
  "figure4.R",
  "figure4_supplement.R",
  "figure5.R"
)

for (script in figure_scripts) {
  path <- here("03_downstream_analysis", "figures", script)
  message("\n--- Sourcing: ", script, " ---")
  source(path)
}

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
message("\n=== All R figures generated ===")
message("Main figures:          ", OUTPUT_MAIN)
message("Supplementary figures: ", OUTPUT_SUP)
message("Session info logs:     ", SESSION_DIR)
message("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("\nNote: Python notebooks (figure5_supplement.ipynb, figure5g.ipynb)")
message("      must be run separately — see README.md for instructions.")
