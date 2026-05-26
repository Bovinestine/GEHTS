# =============================================================================
# processing.py — Data loading and normalization for the GEHTS pipeline
#
# Main entry points:
#   process_data(input_dir, save=False) -> AnnData
#   normalize_by_ubc(adata)            -> AnnData
# =============================================================================

from pathlib import Path
import os
from glob import glob

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from config import (
    DRUG_NAMES,
    GENES_ANABOLIC, GENES_CATABOLIC, GENES_HOUSEKEEPING,
    LOOKUP_DOSE_1, LOOKUP_DOSE_01, LOOKUP_DOSE_10,
)

# Build lookup dictionaries once at import time
_LOOKUP_DICT_1  = {num: DRUG_NAMES[i] + '_1'   for i, num in enumerate(LOOKUP_DOSE_1)}
_LOOKUP_DICT_01 = {num: DRUG_NAMES[i] + '_0.1' for i, num in enumerate(LOOKUP_DOSE_01)}
_LOOKUP_DICT_10 = {num: DRUG_NAMES[i] + '_10'  for i, num in enumerate(LOOKUP_DOSE_10)}


def _get_lookup_table(filename: str) -> dict | None:
    category = filename[7:9]
    if category == '01':
        return _LOOKUP_DICT_1
    elif category == '10':
        return _LOOKUP_DICT_10
    elif category == '-1':
        return _LOOKUP_DICT_01
    return None


def _replace_drug_number(row: pd.Series, col: str) -> object:
    lookup = _get_lookup_table(row['file'])
    if lookup:
        return lookup.get(row[col], row[col])
    return row[col]


def process_data(input_dir, save: bool = False) -> AnnData:
    """Load sorted CSV files from input_dir, annotate barcodes, and return AnnData.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing sorted raw gene-count CSV files.
    save : bool
        If True, write the merged AnnData to <input_dir>/merged.h5ad.
    """
    input_dir = Path(input_dir)
    files = glob(str(input_dir / '*.csv'))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {input_dir}")

    df_list = []
    for file in files:
        df = pd.read_csv(file)
        df['file'] = os.path.basename(file)
        df_list.append(df)
    df = pd.concat(df_list, ignore_index=True)

    if '*drug no' not in df.columns:
        df['*drug no'] = np.nan

    metadata_cols = ['file', 'GeneName', '*drug no', '*drug1 no', '*drug2 no']
    existing_meta = [c for c in metadata_cols if c in df.columns]

    adata = AnnData(
        X   = df.loc[:, ~df.columns.isin(existing_meta)],
        obs = df[existing_meta],
    )

    # Replace barcode numbers with drug name strings
    for col in ['*drug no', '*drug1 no', '*drug2 no']:
        if col in adata.obs.columns:
            adata.obs[col] = adata.obs.apply(
                lambda row, c=col: _replace_drug_number(row, c), axis=1
            )

    # Combine drug1 + drug2 into a single '*drug no' for combination wells
    if '*drug1 no' in adata.obs.columns and '*drug2 no' in adata.obs.columns:
        combined = (
            adata.obs['*drug1 no'].astype(str)
            + '&'
            + adata.obs['*drug2 no'].astype(str)
        )
        adata.obs.loc[adata.obs['*drug no'].isnull(), '*drug no'] = combined
        adata.obs['*drug1 no'] = adata.obs['*drug1 no'].astype(str)
        adata.obs['*drug2 no'] = adata.obs['*drug2 no'].astype(str)

    adata.obs['*drug no'] = adata.obs['*drug no'].astype(str).astype('category')

    # Distinguish control (drug_no == 0, filename contains 'control' or 'c.csv')
    # from inflammatory baseline (drug_no == 0, all other files)
    is_control = (adata.obs['*drug no'] == '0') & (
        adata.obs['file'].str.lower().str.contains(r'c\.csv$|control\.csv$', regex=True)
    )
    adata.obs['drug_condition'] = np.where(
        is_control,
        'control',
        np.where(adata.obs['*drug no'] == '0', 'inflammatory',
                 adata.obs['*drug no'].astype(str))
    )

    adata.obs_names = adata.obs['GeneName']
    adata.obs.rename(columns={'GeneName': 'well_no'}, inplace=True)

    # Remove uninformative probe entries
    mask = ~adata.var.index.isin(['NNNN', 'Homomer'])
    adata = adata[:, mask]

    # Annotate gene categories: 1=anabolic, 2=catabolic, 3=housekeeping, 0=other
    adata.var['gene_idx'] = np.where(
        adata.var.index.isin(GENES_ANABOLIC), 1,
        np.where(adata.var.index.isin(GENES_CATABOLIC), 2,
                 np.where(adata.var.index.isin(GENES_HOUSEKEEPING), 3, 0))
    )

    adata.obs['id'] = adata.obs['file'].astype(str) + adata.obs['well_no'].astype(str)

    if save:
        out_path = input_dir / 'merged.h5ad'
        adata.write_h5ad(str(out_path))
        print(f"Saved: {out_path}")

    return adata


def normalize_by_ubc(adata: AnnData) -> AnnData:
    """Normalize each cell's counts by its Ubc housekeeping count.

    Cells with zero Ubc expression are dropped.
    """
    ubc_index = adata.var.index.get_loc('Ubc')
    ubc_counts = adata.X[:, ubc_index]

    nonzero_mask = np.asarray(ubc_counts != 0).ravel()
    adata = adata[nonzero_mask].copy()

    ubc_index = adata.var.index.get_loc('Ubc')
    adata.X = adata.X.astype(np.float64)
    adata.X = adata.X / adata.X[:, ubc_index][:, None]

    return adata
