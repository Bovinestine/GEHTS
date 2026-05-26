# =============================================================================
# data_loading.R — Load all pre-computed Seurat objects and ML model outputs
#
# Prerequisites: config.R must be sourced first (defines DATA_DIR, PRED_DIR).
# All paths are relative to the project root via here::here().
# =============================================================================

library(SeuratDisk)

.load_rdata <- function(path, description) {
  tryCatch(
    load(path, envir = .GlobalEnv),
    error = function(e) stop("Failed to load ", description, "\n  Path: ", path,
                             "\n  Error: ", conditionMessage(e))
  )
  message("Loaded: ", description)
}

.load_rds <- function(path, description) {
  tryCatch(
    readRDS(path),
    error = function(e) stop("Failed to load ", description, "\n  Path: ", path,
                             "\n  Error: ", conditionMessage(e))
  )
}

# -----------------------------------------------------------------------------
# 1. Raw (un-normalized) Seurat objects
# -----------------------------------------------------------------------------
.load_rdata(file.path(DATA_DIR, "seurat_mac.pub.Rdata"), "mac.pub — primary MAC (raw)")
.load_rdata(file.path(DATA_DIR, "seurat_sin.pub.Rdata"), "sin.pub — single-drug screen (raw)")
.load_rdata(file.path(DATA_DIR, "seurat_cmb.pub.Rdata"), "cmb.pub — combo-drug screen (raw)")
.load_rdata(file.path(DATA_DIR, "seurat_rna.pub.Rdata"), "rna.pub — bulk RNA-seq (raw)")

# -----------------------------------------------------------------------------
# 2. SCTransform-normalized Seurat objects
# -----------------------------------------------------------------------------
.load_rdata(file.path(DATA_DIR, "seurat_mac.sct.Rdata"), "mac.sct — primary MAC (SCT)")
.load_rdata(file.path(DATA_DIR, "seurat_sin.sct.Rdata"), "sin.sct — single-drug screen (SCT)")
.load_rdata(file.path(DATA_DIR, "seurat_cmb.sct.Rdata"), "cmb.sct — combo-drug screen (SCT)")
.load_rdata(file.path(DATA_DIR, "seurat_rna.sct.Rdata"), "rna.sct — bulk RNA-seq (SCT)")

# -----------------------------------------------------------------------------
# 3. External RNA-seq reference (Sebastian 2021, healthy MAC)
# -----------------------------------------------------------------------------
macSebastian <- tryCatch(
  readRDS(SEBASTIAN_RDS),
  error = function(e) {
    warning("macSebastian not loaded — file not found: ", SEBASTIAN_RDS)
    NULL
  }
)
if (!is.null(macSebastian)) {
  macSebastian$drug_condition <- COND_CONTROL
  message("Loaded: macSebastian — Sebastian 2021 healthy MAC RNA-seq")
}

# -----------------------------------------------------------------------------
# 4. ML efficacy predictions (Random Forest, dose 0.1 and 10 μM)
# -----------------------------------------------------------------------------
.load_rdata(file.path(PRED_DIR, "prob_sng10_rf.Rdata"),  "pred_single_10   — RF predictions, single drug, 10 μM")
.load_rdata(file.path(PRED_DIR, "prob_sng0.1_rf.Rdata"), "pred_single_0.1  — RF predictions, single drug, 0.1 μM")
.load_rdata(file.path(PRED_DIR, "prob_cmb0.1_rf.Rdata"), "pred_combo_0.1   — RF predictions, combo, 0.1 μM")

# -----------------------------------------------------------------------------
# 5. Trained prediction models
# -----------------------------------------------------------------------------
en_model  <- .load_rds(file.path(PRED_DIR, "Elastic Net_model.rds"),  "en_model  — Elastic Net")
rf_model  <- .load_rds(file.path(PRED_DIR, "Random Forest_model.rds"), "rf_model  — Random Forest")
xgb_model <- .load_rds(file.path(PRED_DIR, "XGBoost_model.rds"),      "xgb_model — XGBoost")
svr_model <- .load_rds(file.path(PRED_DIR, "SVR_model.rds"),          "svr_model — SVR")
gpr_model <- .load_rds(file.path(PRED_DIR, "GPR_model.rds"),          "gpr_model — GPR")

# -----------------------------------------------------------------------------
# 6. Universal radar chart scale (pre-computed for Figure 4e)
# -----------------------------------------------------------------------------
.load_rdata(file.path(DATA_DIR, "uni_scale.Rdata"), "uni_scale — universal radar chart scale")

message("\nData loading complete.")
