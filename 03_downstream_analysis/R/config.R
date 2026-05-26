# =============================================================================
# config.R — Central parameter store for GEHTS downstream analysis
#
# Source this file at the top of every script:
#   source(here::here("03_downstream_analysis", "R", "config.R"))
#
# All gene panels, GO module definitions, analysis constants, color palettes,
# and file paths are defined here once to ensure consistency across figures.
# =============================================================================

library(here)
library(ggsci)

# -----------------------------------------------------------------------------
# 1. Gene panels (GEHTS 30-gene chip)
# -----------------------------------------------------------------------------
GENES_ANABOLIC <- c(
  "Acan", "Sox9", "Col2a1", "Matn1", "Matn3", "Ucma",
  "Ccnd3", "Gadd45g", "Pth1r", "Gm26633", "Col27a1"
)

GENES_CATABOLIC <- c(
  "Mmp3", "Mmp13", "Il6", "Il17b", "Adamts5", "Igfbp3",
  "Ccl2", "Cxcl5", "Cxcl1", "Fosl2", "Tlr2", "Tnfrsf1b"
)

GENES_HOUSEKEEPING <- c(
  "Hprt", "Actb", "Gapdh", "B2m", "Ubc", "Ppia", "Rpl23"
)

GENES_ALL <- c(GENES_ANABOLIC, GENES_CATABOLIC, GENES_HOUSEKEEPING)

# -----------------------------------------------------------------------------
# 2. GO module definitions (6 biological pathway gene sets)
#    Used by: radar charts, trajectory plots, bubble plots, violin plots
# -----------------------------------------------------------------------------
GO_MODULES <- list(
  Chondrocyte_Development  = c("Matn1", "Acan", "Sox9", "Col27a1"),
  Cartilage_Condensation   = c("Col2a1", "Acan", "Sox9"),
  NO_Synthase_Upregulation = c("Tlr2", "Ccl2"),
  ECM_Disassembly          = c("Adamts5", "Mmp13"),
  Inflammatory_response    = c("Tlr2", "Ccl2", "Il17b", "Tnfrsf1b",
                               "Cxcl1", "Il6", "Cxcl5", "Fosl2"),
  Tissue_homeostasis       = c("Sox9", "Col2a1", "Pth1r", "Fosl2")
)

# Anabolic and catabolic axis membership for 2D GO coordinate plots
GO_ANABOLIC_TERMS  <- c("Chondrocyte_Development",
                         "Cartilage_Condensation",
                         "Tissue_homeostasis")
GO_CATABOLIC_TERMS <- c("NO_Synthase_Upregulation",
                         "ECM_Disassembly",
                         "Inflammatory_response")

# Cellular state modules used in Figure 5f dot plot
STATE_MODULES <- list(
  State_Viability = c("Hprt", "Actb", "Gapdh", "B2m", "Ubc", "Ppia", "Rpl23"),
  State_Stress    = c("Gadd45g", "Igfbp3"),
  State_Prolif    = c("Ccnd3")
)

# -----------------------------------------------------------------------------
# 3. Analysis parameters
# -----------------------------------------------------------------------------
SEED        <- 123     # global RNG seed for all ML models and dimension reductions
LOG2FC_CAP  <- 3.0    # symmetric cap applied to Log2FC values in volcano/dot plots
PSEUDOCOUNT <- 0.01   # added before log2 to avoid log(0)
N_FOLDS_CV  <- 5      # number of cross-validation folds
DOSE_LEVELS <- c(0.1, 10)  # μM doses screened

# Baseline condition labels used consistently across all functions
COND_CONTROL     <- "control"
COND_DISEASE     <- "inflammatory"

# -----------------------------------------------------------------------------
# 4. File paths (all relative via here::here() — no absolute paths)
# -----------------------------------------------------------------------------
DATA_DIR    <- here::here("Rdata", "pub")
PRED_DIR    <- here::here("Rdata", "pub", "prediction")
OUTPUT_MAIN <- here::here("03_downstream_analysis", "output", "main")
OUTPUT_SUP  <- here::here("03_downstream_analysis", "output", "supplementary")
SESSION_DIR <- here::here("03_downstream_analysis", "session_info")

# CSV data file for Figure 1
PLATFORM_CSV <- here::here("03_downstream_analysis",
                            "Transcriptomic_Screening_Platform_Comparison.csv")

# External reference RNA-seq
SEBASTIAN_RDS <- here::here("RNA-seq", "macCtrl_sebastian.rds")

# -----------------------------------------------------------------------------
# 5. Color palettes
# -----------------------------------------------------------------------------
PALETTE_NPG <- pal_npg("nrc")(10)

# Locked condition colors (used consistently across all figures)
COL_CONTROL     <- PALETTE_NPG[2]   # NPG Blue
COL_DISEASE     <- PALETTE_NPG[1]   # NPG Red
COL_DOSE_LOW    <- "#A6CEE3"        # Sky blue  (0.1 μM)
COL_DOSE_HIGH   <- "#1F78B4"        # Steel blue (10 μM)

# Gene-type annotation colors (heatmap row annotations)
COL_ANABOLIC    <- "forestgreen"
COL_CATABOLIC   <- "orange"
COL_HOUSEKEEP   <- "grey60"

# Chip region colors (edge effect figures)
COL_CENTER      <- "#66c2a5"
COL_EDGE        <- "#fc8d62"
COL_CORNER      <- "#8da0cb"

# -----------------------------------------------------------------------------
# 6. Helper: ensure output directories exist
# -----------------------------------------------------------------------------
.ensure_dirs <- function() {
  dirs <- c(OUTPUT_MAIN, OUTPUT_SUP, SESSION_DIR)
  invisible(lapply(dirs, function(d) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }))
}
.ensure_dirs()

# -----------------------------------------------------------------------------
# 7. Helper: save sessionInfo() to session_info/
# -----------------------------------------------------------------------------
save_session_info <- function(script_name) {
  out_file <- file.path(SESSION_DIR,
                        paste0(gsub("\\.R$", "", script_name), "_session.txt"))
  writeLines(capture.output(sessionInfo()), out_file)
  message("Session info saved to: ", out_file)
}
