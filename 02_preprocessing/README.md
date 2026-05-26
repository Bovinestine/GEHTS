# 02_preprocessing

Python preprocessing pipeline and R normalization for the GEHTS single-cell dataset.

## Prerequisites

**Python** (≥ 3.9) with the following packages:

```
anndata
numpy
pandas
scipy
scanpy
scikit-learn
matplotlib
seaborn
statannotations
```

**R** with the following packages:

```r
install.packages(c("Seurat", "harmony", "dplyr", "here"))
# sceasy (reads h5ad into Seurat):
remotes::install_github("cellgeni/sceasy")
```

---

## Step 1 — Python preprocessing

Converts raw sorted CSV gene-count files into an AnnData object (`.h5ad`).

Run once for each of the three input directories:

```bash
python main.py data_sorted/control&inflam  --save
python main.py data_sorted/single_only     --save
python main.py data_sorted/combi_only      --save
```

Each call writes `merged.h5ad` into the specified directory.

`--save` is optional; omit it to inspect the AnnData in memory without writing to disk.

---

## Step 2 — R normalization (SCTransform)

Reads the three `merged.h5ad` files, merges them, applies SCTransform, and writes
three Seurat objects used by `03_downstream_analysis/`.

```r
# From within 02_preprocessing/:
source("normalization.R")

# Or from the repo root:
Rscript 02_preprocessing/normalization.R
```

**Output** (written to `../Rdata/pub/`):

| File | Contents |
|------|----------|
| `mac.sct.RData` | Control & inflammation reference (`source == "mac"`) |
| `sin.sct.RData` | Single-drug screen (`source == "drg"`) |
| `cmb.sct.RData` | Combination-drug screen (`source == "cmb"`) |

> **reticulate / Python path:** `sceasy` uses reticulate to read `.h5ad` files.
> If R does not find the correct Python automatically, set `RETICULATE_PYTHON`
> before sourcing the script:
> ```r
> Sys.setenv(RETICULATE_PYTHON = "/path/to/your/conda/envs/gehts/bin/python")
> source("normalization.R")
> ```

---

## Customisation

All project-level constants live in `config.py`.

| Constant | Description |
|----------|-------------|
| `GENES_ANABOLIC` | 11-gene anabolic panel |
| `GENES_CATABOLIC` | 12-gene catabolic panel |
| `GENES_HOUSEKEEPING` | 7-gene housekeeping panel |
| `DRUG_NAMES` | Ordered list of 14 drug names |
| `LOOKUP_DOSE_1` / `_01` / `_10` | DLP barcode → drug name mapping per dose tier |

Gene panels in `config.py` must stay in sync with
`03_downstream_analysis/R/config.R`.

To use a different barcode set, edit `DRUG_NAMES` and the three `LOOKUP_DOSE_*`
lists. If the filename convention for dose detection also changes, update
`_get_lookup_table()` in `processing.py`.

---

## File overview

```
02_preprocessing/
├── main.py           — CLI entry point
├── config.py         — gene panels, drug names, barcode lookup tables
├── processing.py     — process_data(), normalize_by_ubc()
├── visualization.py  — heatmap, PCA, t-SNE, UMAP, boxplot, correlation plots
└── normalization.R   — load h5ad → SCTransform → save Seurat objects
```
