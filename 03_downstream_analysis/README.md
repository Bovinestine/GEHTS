# GEHTS Downstream Analysis

Statistical modeling, machine-learning efficacy prediction, synergy analysis, and figure generation for the GEHTS (Gene Expression High-Throughput Screening) study.

---

## System requirements

| Item | Requirement |
|------|-------------|
| Operating system | Windows 10/11, macOS 12+, or Ubuntu 20.04+ |
| R version | ≥ 4.3.0 |
| RAM | ≥ 16 GB (32 GB recommended for full dataset) |
| Python | ≥ 3.9 (Figures 5g and Extended Figure 5) |

R packages used (see `session_info/` for exact versions after running):
`Seurat`, `ggplot2`, `dplyr`, `glmnet`, `randomForest`, `xgboost`, `e1071`, `kernlab`,
`caret`, `Metrics`, `ggrepel`, `ggsci`, `cowplot`, `gridExtra`, `viridis`, `reshape2`,
`fmsb`, `svglite`, `scales`, `tidyr`, `pheatmap`, `ggsignif`, `here`

Python packages (Figures 5g and Extended Figure 5): `scanpy`, `numpy`, `matplotlib`, `scipy`, `pandas`

---

## Installation

Install R package dependencies (estimated time: 10–20 min on a typical machine):

```r
install.packages(c(
  "here", "dplyr", "ggplot2", "ggrepel", "ggsci", "cowplot", "gridExtra",
  "viridis", "reshape2", "fmsb", "svglite", "scales", "tidyr", "pheatmap",
  "ggsignif", "glmnet", "randomForest", "xgboost", "e1071", "kernlab",
  "caret", "Metrics"
))
# Seurat (follow https://satijalab.org/seurat/articles/install_additional)
install.packages("Seurat")
```

For Python notebooks (Figures 5g and Extended Figure 5):

```bash
pip install scanpy numpy matplotlib scipy pandas
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
│   └── pub/
│       ├── mac.sct.RData
│       ├── sin.sct.RData
│       ├── cmb.sct.RData
│       ├── primary.sct.RData
│       ├── atdc5.sct.RData
│       ├── mac_cmb0.1.sct.h5ad
│       └── mac_sin0.1.sct.h5ad
└── 03_downstream_analysis/   ← you are here
```

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

### Python notebooks (Figures 5g and Extended Figure 5)

**Extended Figure 5** — class-summed gene expression barplot (run from `03_downstream_analysis/`):

```bash
cd 03_downstream_analysis
jupyter nbconvert --to notebook --execute figures/figure5_supplement.ipynb
```

**Figure 5g** — per-module synergy barplot (run from `03_downstream_analysis/figures/`):

```bash
cd 03_downstream_analysis/figures
jupyter nbconvert --to notebook --execute figure5g.ipynb
```

---

## Output

| Directory | Contents |
|-----------|----------|
| `output/main/` | Main figure PDFs (Figures 1–5) |
| `output/supplementary/` | Extended Data / supplementary PDFs |
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
| Fig. 2f | `figures/figure2.R` | `plot_gene_divergence()` |
| Ext. Data Fig. 5a | `figures/figure2_supplement.R` | `get_spatial_coords()`, edge-effect boxplot |
| Ext. Data Fig. 5b | `figures/figure2_supplement.R` | spatial uniformity heatmap |
| Ext. Data Fig. 6 | `figures/figure2_supplement.R` | `draw_boxplots_for_genes_with_common_legend()` |
| Ext. Data Fig. 7 | `figures/figure2_supplement.R` | `plot_clr_heatmap()` |
| Fig. 3b | `figures/figure3.R` | `create_heatmap_by_expression_ps()` |
| Fig. 3c | `figures/figure3.R` | `plot_main_umap()` |
| Ext. Data Fig. 8 | `figures/figure3_supplement.R` | `plot_highlighted_umaps()` |
| Ext. Data Fig. 9 | `figures/figure3_supplement.R` | `draw_boxplots_for_genes_with_common_legend()` |
| Fig. 4b | `figures/figure4.R` | `plot_efficacy_barplot()` |
| Fig. 4c | `figures/figure4.R` | `plot_phenotypic_screening()` |
| Fig. 4d | `figures/figure4.R` | `plot_highlighted_umaps()` |
| Fig. 4e | `figures/figure4.R` | `plot_automated_coexpression()` |
| Fig. 4f | `figures/figure4.R` | `plot_radar_fast()` |
| Ext. Data (CV) | `figures/figure4_supplement.R` | 5-fold CV comparison boxplot |
| Ext. Data (violins) | `figures/figure4_supplement.R` | `plot_supplementary_violins()` |
| Fig. 5c | `figures/figure5.R` | `generate_bliss_pdf()` |
| Fig. 5d | `figures/figure5.R` | `plot_bliss_synergy()` |
| Fig. 5e | `figures/figure5.R` | `plot_biological_alignment()` |
| Fig. 5f | `figures/figure5.R` | `plot_dotplot_log2fc()` |
| Fig. 5g | `figures/figure5g.ipynb` | Python: per-module synergy barplot |
| Ext. Data Fig. 10 | `figures/figure5_supplement.ipynb` | Python: `sem_of_class_sum()` |

---

## Code structure

```
03_downstream_analysis/
├── main.R                          ← master entry point (R only)
├── R/
│   ├── config.R                    ← gene panels, palette, paths, seeds
│   ├── utils.R                     ← shared plotting helpers
│   ├── data_loading.R              ← load all data objects
│   ├── scoring_methods.R           ← ML model training (5-fold CV)
│   └── prediction_score.R         ← efficacy prediction
├── figures/
│   ├── figure1.R
│   ├── figure2.R
│   ├── figure2_supplement.R        ← Ext. Data Fig. 5a, 5b, 6, 7
│   ├── figure3.R
│   ├── figure3_supplement.R        ← Ext. Data Fig. 3, 8
│   ├── figure4.R
│   ├── figure4_supplement.R        ← Ext. Data CV comparison, GO violins
│   ├── figure5.R
│   ├── figure5g.ipynb              ← Python: Fig. 5g per-module synergy
│   └── figure5_supplement.ipynb   ← Python: Extended Fig. 5 class-sum barplot
├── output/
│   ├── main/                       ← generated PDFs (main figures)
│   └── supplementary/              ← generated PDFs (supplementary)
├── session_info/                   ← R sessionInfo() and pip freeze logs
└── archive/                        ← superseded scripts (not part of pipeline)
```

---

## Reproducibility notes

- All R random seeds are fixed via `SEED <- 123` in `R/config.R`.
- ML cross-validation uses `set.seed(SEED + fold_index)` per fold so fold assignment is deterministic.
- Running `main.R` twice from a clean R session should produce byte-identical output PDFs.
- The `session_info/` directory captures the exact R and Python package versions used for each script.
