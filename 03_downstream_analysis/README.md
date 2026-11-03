# GEHTS Downstream Analysis

Statistical modeling, machine-learning efficacy prediction, synergy analysis, and figure generation for the GEHTS (Gene Expression High-Throughput Screening) study.

---

## System requirements

| Item | Requirement |
|------|-------------|
| Operating system | Windows 10/11, macOS 12+, or Ubuntu 20.04+ |
| R version | ≥ 4.3.0 (tested on 4.4.3) |
| RAM | ≥ 16 GB (32 GB recommended for full dataset) |
| Python | ≥ 3.9 (Fig. 5g, Extended Fig. 5, and the CSV/test tooling) |

R packages used (see `session_info/` for exact versions after running):
`Seurat`, `ggplot2`, `dplyr`, `glmnet`, `randomForest`, `xgboost`, `e1071`, `kernlab`,
`caret`, `Metrics`, `ggrepel`, `ggsci`, `cowplot`, `gridExtra`, `viridis`, `reshape2`,
`fmsb`, `svglite`, `scales`, `tidyr`, `pheatmap`, `ggsignif`, `jsonlite`, `here`

GSEA module-derivation scripts (`R/gsea_panel_modules*.R`, see [Module definitions](#module-definitions)) additionally need `msigdbr`, `fgsea`, `BiocParallel` — not required for the main figure pipeline.

`R/build_rna_pub_from_featurecounts.R` (see [rna.pub / atdc5.pub provenance](#rnapub--atdc5pub-provenance-fig-2c-2f)) additionally needs `DESeq2` and `msigdbr` — provenance only, not required for the main figure pipeline.

Python packages: `numpy`, `matplotlib`, `scipy`, `pandas`, `jupyter`. `scanpy` is needed only by `figure5_supplement.ipynb` (Extended Fig. 5, still h5ad-based); `figure5g.ipynb` and `figure5g_test.ipynb` read CSV and do not need it.

---

## Installation

### Recommended: conda environment

A verified environment spec is at the project root:

```bash
conda env create -f environment.yml
conda activate gehts
```

### Manual R installation

```r
install.packages(c(
  "here", "dplyr", "ggplot2", "ggrepel", "ggsci", "cowplot", "gridExtra",
  "viridis", "reshape2", "fmsb", "svglite", "scales", "tidyr", "pheatmap",
  "ggsignif", "jsonlite", "glmnet", "randomForest", "xgboost", "e1071",
  "kernlab", "caret", "Metrics"
))
# Seurat (follow https://satijalab.org/seurat/articles/install_additional)
install.packages("Seurat")
```

### Manual Python installation

```bash
pip install numpy matplotlib scipy pandas jupyter
pip install scanpy   # only needed for figure5_supplement.ipynb
```

---

## Data availability

Input data files (`.RData`, `.rds`, `.h5ad`) are deposited at:
- **GEO accession**: [to be added upon acceptance]
- **Zenodo DOI**: [to be added upon acceptance]

Download and place the files under `../Rdata/pub/` relative to this directory:

```
GEHTS/
├── Rdata/
│   ├── pub/
│   │   ├── seurat_mac.sct.Rdata / seurat_sin.sct.Rdata / seurat_cmb.sct.Rdata
│   │   ├── seurat_mac.pub.Rdata / seurat_sin.pub.Rdata / seurat_cmb.pub.Rdata
│   │   ├── seurat_rna.pub.Rdata
│   │   ├── seurat_atdc5.pub.Rdata ← ATDC5 cell line (GSE269585 37C baseline), Fig. 2f
│   │   ├── mac_sin0.1_scores.csv / mac_cmb0.1_scores.csv 
│   │   └── prediction/           ← trained models + cached predictions
│   │       ├── cv_results_per_fold.csv
│   │       ├── Eleastic_net_models.rds / GPR_model.rds / Random_forest_model.rds / SVR_model.rds / XGBoost_model.rds
│   │       └── pred_single_10_rf.Rdata / pred_single_0.1_rf.Rdata / pred_combo_0.1_rf.Rdata
│   └── RNAseq/
│       ├── GSE269585_Raw_count
│       ├── ctrl-control / ctrl2-condtrol / ctrl3-control
│       └── il1b-il1b / il1b2-il1b / il1b3-il1b
└── 03_downstream_analysis/   ← you are here
```

---

## Module definitions

`R/config.R`'s `GO_MODULES` (7 gene modules) is the single source of truth for every module-score figure (Fig. 4e, Fig. 5f, Fig. 5g). It was derived by running `fgsea` against an external basal-vs-IL-1β reference — see `R/gsea_panel_modules*.R` and their outputs in `output/gsea/` for the full derivation and methods notes.

- `GO_ANABOLIC_TERMS` / `GO_CATABOLIC_TERMS`: 3 modules each, assigned to the anabolic/catabolic axes used in 2D coordinate plots and the anabolic-minus-catabolic index.
- The 7th module (`Programmed_Cell_Death`) is a state/triage module: displayed in every module figure but not assigned to either axis.
- Every module figure scores from **SCT `scale.data` (Pearson residuals)**, not the log-normalized `data` slot — residual scores are centered near zero and go negative, so figures built from them use a plain baseline *difference* on a diverging palette, not log2FC.
- `config.R` exports `GO_MODULES`/`GO_ANABOLIC_TERMS`/`GO_CATABOLIC_TERMS` to `R/go_modules.json` on every source (via `jsonlite`). `figure5g.ipynb` and `figure5g_test.ipynb` read that JSON instead of hand-copying gene lists, so the R and Python pipelines cannot drift apart.

Fig. 5e (`plot_biological_alignment()`) is the one figure in this set that does **not** use `GO_MODULES` — its Anabolic/Catabolic score is the mean expression (SCT `data` slot) of the raw `GENES_ANABOLIC`/`GENES_CATABOLIC` gene panels, the same panels Fig. 4c uses.

---

## rna.pub / atdc5.pub provenance (Fig. 2c, 2f, Fig. 2f's supplementary companion)

`rna.pub` (bulk RNA-seq, control vs. IL-1β, n=3/condition) and `atdc5.pub` (ATDC5 cell line, GSE269585 37C baseline, n=4) are built from raw counts by two one-time, not-part-of-the-main-run scripts (same category as `R/gsea_panel_modules*.R` — provenance, re-run only if the source data changes):

- **`R/build_rna_pub_from_featurecounts.R`** — merges 6 raw featureCounts outputs (`Rdata/RNAseq/{ctrl,ctrl2,ctrl3,il1b,il1b2,il1b3}-{control,il1b}`; GRCm38.75, Ensembl gene IDs), maps Ensembl ID → mouse gene symbol via `msigdbr`'s mouse-native gene tables (offline stand-in for org.Mm.eg.db/biomaRt, unavailable in this environment; ~82% of genes get a curated symbol), and writes `Rdata/pub/seurat_rna.pub.Rdata`. Cells are named `control`/`control.1`/`control.2`/`inflammatory`/`inflammatory.1`/`inflammatory.2`. QC (DESeq2 vst + PCA) and a methods note are written to `output/rnaseq_build/`. The previous 2-replicate `rna.pub` is backed up once to `Rdata/pub/seurat_rna.pub.legacy_n2.Rdata`.
- **`R/build_atdc5_from_gse269585.R`** — reads `Rdata/RNAseq/GSE269585_Raw_Count.txt`, keeps only the 4 `37C_*` columns (the only condition in that dataset that is a valid baseline counterpart to `rna.pub`'s control condition — 40C/44C are an unrelated temperature-response axis), and writes `Rdata/pub/seurat_atdc5.pub.Rdata`.

Figure 2c's `SCTransform(rna.pub, ...)` call uses `min_cells = 1` (not Seurat's default of 5): `rna.pub`'s "cells" are bulk RNA-seq samples, not single cells, so a gene detected in only 1-4 of 6 samples (e.g. `Il6`, near-absent in control) is real biology, not single-cell dropout, and the default filter would otherwise silently drop panel genes.

Figure 2f (`plot_gene_divergence()`, utils.R) uses **mac.sct** (GE-HTS chip, control wells) as the primary side vs. `atdc5.pub`, both CLR-normalized from raw RNA-assay counts. `figure2_supplement.R` carries an unnumbered companion supplementary figure using **rna.pub** (bulk RNA-seq, control replicates) as the primary side instead, so the comparison is reported on both platforms — same shared function, different primary data source, output to `output/supplementary/Sfig_rna_vs_atdc5_lollipop.pdf`.

---

## Running the analysis

### R figures (all main and supplementary)

Run from the **project root** (`GEHTS/`):

```r
Rscript 03_downstream_analysis/main.R
```

Or interactively in R:

```r
source(here::here("03_downstream_analysis", "main.R"))
```

Expected runtime: 30–60 minutes on a machine with 16 GB RAM.

### Python notebooks

**Fig. 5g** — per-module synergy barplot (run from `03_downstream_analysis/figures/`). Needs the CSV exports from `R/export_combo_test_data.R` first:

```bash
Rscript 03_downstream_analysis/R/export_combo_test_data.R
cd 03_downstream_analysis/figures
jupyter nbconvert --to notebook --execute --inplace figure5g.ipynb
```

**Extended Figure 5** — class-summed gene expression barplot (run from `03_downstream_analysis/`, needs the h5ad files):

```bash
cd 03_downstream_analysis
jupyter nbconvert --to notebook --execute figures/figure5_supplement.ipynb
```

---

## Output

| Directory | Contents |
|-----------|----------|
| `output/main/` | Main figure PDFs (Figures 1–5) |
| `output/supplementary/` | Extended Data / supplementary PDFs |
| `output/gsea/` | GSEA module-derivation outputs (QC plots, panel fold-changes, methods notes) — provenance for `GO_MODULES`, not part of the main figure run |
| `output/test/` | Test-harness outputs (`R/test_new_modules.R`, `R/export_combo_test_data.R`, `figures/figure5g_test.ipynb`) — see [Testing & validation](#testing--validation); never overwrites `output/main`/`output/supplementary` |
| `session_info/` | `sessionInfo()` snapshots per R script; `pip freeze` per Python notebook |

---

## Figure-to-script mapping

| Figure | Script | Key function(s) |
|--------|--------|-----------------|
| Fig. 1 | `figures/figure1.R` | `plot_platform_comparison()` |
| Fig. 2b | `figures/figure2.R` | `create_heatmap_by_expression_ps()` |
| Fig. 2c | `figures/figure2.R` | `create_correlation_plots()` |
| Fig. 2d | `figures/figure2.R` | `plot_dim_reduction()` |
| Fig. 2e | `figures/figure2.R` | `generate_volcano_plot()` |
| Fig. 2f | `figures/figure2.R` | `plot_gene_divergence()` (defined in `R/utils.R`) |
| Ext. Data Fig. 5a | `figures/figure2_supplement.R` | `draw_boxplots_for_gene_sums_with_common_legend()` — class-summed anabolic/catabolic gene counts |
| Ext. Data Fig. 5b | `figures/figure2_supplement.R` | `get_spatial_coords()`, edge-effect boxplot |
| Ext. Data Fig. 5c | `figures/figure2_supplement.R` | spatial uniformity heatmap |
| Ext. Data Fig. 6 | `figures/figure2_supplement.R` | `draw_boxplots_for_genes_with_common_legend()` |
| Ext. Data Fig. 7 | `figures/figure2_supplement.R` | `plot_clr_heatmap()` |
| Fig. 3b | `figures/figure3.R` | `create_heatmap_by_expression_ps()` |
| Fig. 3c | `figures/figure3.R` | `plot_main_umap()` |
| Ext. Data Fig. 8 | `figures/figure3_supplement.R` | `plot_highlighted_umaps()` |
| Ext. Data Fig. 9 | `figures/figure3_supplement.R` | `draw_boxplots_for_genes_with_common_legend()` |
| Fig. 4b | `figures/figure4.R` | `plot_efficacy_barplot()` |
| Fig. 4c | `figures/figure4.R` | `plot_phenotypic_screening()` |
| Fig. 4d | `figures/figure4.R` | `plot_automated_coexpression()` |
| Fig. 4e | `figures/figure4.R` | `plot_dotplot_diff()` (utils.R) — 7 GO modules, Δ residual vs. IL-1β, dot size = Wilcoxon/BH significance |
| Ext. Data (CV) | `figures/figure4_supplement.R` | 5-fold CV comparison boxplot |
| Ext. Data (violins) | `figures/figure4_supplement.R` | `plot_supplementary_violins()` |
| Fig. 5c | `figures/figure5.R` | `generate_bliss_pdf()` |
| Fig. 5d | `figures/figure5.R` | `plot_bliss_synergy()` |
| Fig. 5e | `figures/figure5.R` | `plot_biological_alignment()` — GENES_ANABOLIC/GENES_CATABOLIC gene-panel means, not GO_MODULES |
| Fig. 5f | `figures/figure5.R` | `plot_dotplot_diff()` (utils.R) — same function as Fig. 4e |
| Fig. 5g | `figures/figure5g.ipynb` | Python: per-module synergy barplot, reads CSV via `export_combo_test_data.R` |
| Ext. Data Fig. 10 | `figures/figure5_supplement.ipynb` | Python: `sem_of_class_sum()` |

---

## Code structure

```
03_downstream_analysis/
├── main.R                          ← master entry point (R only)
├── R/
│   ├── config.R                    ← gene panels, GO_MODULES, palette, paths, seeds;
│   │                                  exports R/go_modules.json on every source
│   ├── utils.R                     ← shared plotting helpers, incl. residual-based
│   │                                  module scoring (compute_module_scores_long(),
│   │                                  plot_dotplot_diff())
│   ├── data_loading.R              ← load all data objects
│   ├── build_rna_pub_from_featurecounts.R ← rebuilds rna.pub from raw
│   │                                  featureCounts (provenance, not part of
│   │                                  the main figure run — see "rna.pub /
│   │                                  atdc5.pub provenance" above)
│   ├── build_atdc5_from_gse269585.R ← builds atdc5.pub from GSE269585 (same,
│   │                                  provenance only)
│   ├── scoring_methods.R           ← ML model training (5-fold CV)
│   ├── prediction_score.R          ← efficacy prediction
│   ├── go_modules.json             ← GO_MODULES, exported by config.R (R + Python
│   │                                  single source of truth)
│   ├── gsea_panel_modules*.R       ← fgsea derivation of GO_MODULES (provenance,
│   │                                  not part of the main figure run);
│   │                                  gsea_panel_modules_external_deg_vst.R is
│   │                                  the final pipeline (only DEG source used
│   │                                  for the checked-in GO_MODULES) and
│   │                                  includes its own cartilage-suppression
│   │                                  diagnostic as its last step
│   └── export_combo_test_data.R    ← exports CSV (SCT residuals + drug_name) for
│                                      figure5g.ipynb / figure5g_test.ipynb, since the
│                                      h5ad inputs aren't in this repo
├── figures/
│   ├── figure1.R
│   ├── figure2.R
│   ├── figure2_supplement.R        ← Ext. Data Fig. 5a, 5b, 5c, 6, 7
│   ├── figure3.R
│   ├── figure3_supplement.R        ← Ext. Data Fig. 3, 8
│   ├── figure4.R
│   ├── figure4_supplement.R        ← Ext. Data CV comparison, GO violins
│   ├── figure5.R
│   ├── figure5g.ipynb              ← Python: Fig. 5g per-module synergy (CSV-based)
│   └── figure5_supplement.ipynb    ← Python: Extended Fig. 5 class-sum barplot (h5ad-based)
├── output/
│   ├── main/                       ← generated PDFs (main figures)
│   ├── supplementary/              ← generated PDFs (supplementary)
│   ├── gsea/                       ← GO_MODULES derivation outputs (provenance)
│   └── test/                       ← test-harness outputs (never touches main/supplementary)
└── session_info/                   ← R sessionInfo() and pip freeze logs

```

---

## Testing & validation

Two independent pieces of tooling validate the module set and keep R/Python in sync, without ever touching `output/main` or `output/supplementary`:

- **`R/test_new_modules.R`** — sources `config.R`/`utils.R`, then runs, in order: (1) definition checks on `GO_MODULES` (every gene in the 30-gene panel, no anabolic/catabolic cross-contamination, no housekeeping gene in a module, every gene present in `scale.data` of `mac.sct`/`sin.sct`/`cmb.sct`); (2) a direction gate — scores all modules on the real chip (`mac.sct`, basal vs. IL-1β) and asserts anabolic modules move down / catabolic and state modules move up, since `GO_MODULES` was derived from an *external* dataset, so this is a genuine out-of-sample check; (3) regenerates Fig. 4e/5f/Supp. Fig. 2 with the current module set into `output/test/`. Fails loudly (prints a pass/fail table, stops) on any check failure — this is what caught the original `Response_To_IL1` mis-composition (see `output/gsea/` methods notes and `config.R`'s inline comment on that module).
- **`R/export_combo_test_data.R`** + **`figures/figure5g_test.ipynb`** — exports per-cell SCT residuals to CSV (since the h5ad inputs `figure5g.ipynb` was originally written for aren't available in this repo) and re-derives Fig. 5g's per-module synergy statistic from them, cross-checked numerically against the R Fig. 5f dot plot's values for the same combination.

---

## Reproducibility notes

- All R random seeds are fixed via `SEED <- 123` in `R/config.R`.
- ML cross-validation uses `set.seed(SEED + fold_index)` per fold so fold assignment is deterministic.
- Running `main.R` twice from a clean R session should produce byte-identical output PDFs.
- `figure5g.ipynb`'s bootstrap synergy test uses `np.random.default_rng(42)` (not the legacy `np.random.seed`/`np.random.choice` API).
- The `session_info/` directory captures the exact R and Python package versions used for each script.
