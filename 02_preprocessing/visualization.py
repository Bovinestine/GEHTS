# =============================================================================
# visualization.py — Plotting utilities for the GEHTS preprocessing pipeline
#
# All functions accept an AnnData object (or a DataFrame for legacy heatmaps).
# Functions are grouped by type:
#   1. Heatmaps
#   2. Dimensionality reduction plots
#   3. QC / correlation plots
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from scipy import sparse
from anndata import AnnData
import scanpy as sc
import itertools
from statannotations.Annotator import Annotator

from config import GENES_ANABOLIC, GENES_CATABOLIC, GENES_HOUSEKEEPING

# =============================================================================
# 1. Heatmaps
# =============================================================================

def plot_heatmap(adata: AnnData) -> None:
    """Log2-scale heatmap of all genes in AnnData (cells × genes)."""
    adata.var.sort_values(['gene_idx'], inplace=True)
    sorted_order = np.argsort(adata.var.index)
    X = adata.X[:, sorted_order]
    log_X = np.log2(X + 1)

    plt.figure(figsize=(10, 10))
    sns.heatmap(log_X, cmap='viridis',
                xticklabels=adata.obs.index,
                yticklabels=adata.var.index)
    plt.title('Heatmap of Gene Counts (log2)')
    plt.tight_layout()
    plt.show()


def draw_heatmap_separated_lin(
    data: pd.DataFrame,
    anabolic_genes=GENES_ANABOLIC,
    catabolic_genes=GENES_CATABOLIC,
    housekeeping_genes=GENES_HOUSEKEEPING,
) -> None:
    """Linear-scale heatmap with three gene-category panels, proportional row heights."""
    if data.index.name != 'Unnamed: 0':
        data = data.set_index('Unnamed: 0')

    ana  = data.loc[anabolic_genes]
    cata = data.loc[catabolic_genes]
    hk   = data.loc[housekeeping_genes]
    total_rows = len(ana) + len(cata) + len(hk)

    fig = plt.figure(figsize=(8, total_rows))
    gs  = fig.add_gridspec(3, 1, height_ratios=[len(ana), len(cata), len(hk)])

    for ax, panel_data, label in zip(
        [fig.add_subplot(gs[i]) for i in range(3)],
        [ana, cata, hk],
        ['Anabolic Genes', 'Catabolic Genes', 'Housekeeping Genes'],
    ):
        sns.heatmap(panel_data, cmap='viridis', ax=ax)
        ax.yaxis.set_label_position('right')
        ax.set_ylabel(label)
        ax.set_yticks(range(len(panel_data)))
        ax.set_yticklabels(panel_data.index, rotation=0)
        ax.set_xticks([])

    fig.tight_layout(pad=1.5)
    plt.subplots_adjust(hspace=1)
    plt.show()


def draw_heatmap_separated_wohk(
    data: pd.DataFrame,
    anabolic_genes=GENES_ANABOLIC,
    catabolic_genes=GENES_CATABOLIC,
) -> None:
    """Linear-scale heatmap with anabolic and catabolic panels only (no housekeeping)."""
    if data.index.name != 'Unnamed: 0':
        data = data.set_index('Unnamed: 0')

    ana  = data.loc[anabolic_genes]
    cata = data.loc[catabolic_genes]

    plt.figure(figsize=(9, 12))
    plt.rcParams.update({'font.size': 8})

    plt.subplot(2, 1, 1)
    sns.heatmap(ana, cmap='viridis')
    plt.title('Anabolic Genes Expression Heatmap')
    plt.ylabel('Anabolic Genes')
    plt.xticks([])

    plt.subplot(2, 1, 2)
    sns.heatmap(cata, cmap='viridis')
    plt.title('Catabolic Genes Expression Heatmap')
    plt.ylabel('Catabolic Genes')
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(hspace=1)
    plt.show()


def draw_heatmap_lin(data: pd.DataFrame) -> None:
    """Linear-scale heatmap for an arbitrary gene × sample DataFrame."""
    if data.index.name != 'Unnamed: 0':
        data = data.set_index('Unnamed: 0')

    plt.figure(figsize=(8, 9))
    plt.rcParams.update({'font.size': 8})
    sns.heatmap(data, cmap='viridis')
    plt.title('Gene Expression Heatmap (Linear Scale)')
    plt.tight_layout(pad=3.0)
    plt.show()


# =============================================================================
# 2. Dimensionality reduction plots
# =============================================================================

def plot_pca(adata: AnnData) -> AnnData:
    """Normalize, log-transform, scale, and plot PCA colored by drug_condition."""
    sc.pp.normalize_total(adata, target_sum=1e3)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    sc.pl.pca(adata, color='drug_condition')
    plt.tight_layout()
    return adata


def plot_tsne(adata: AnnData, random_state: int = 0) -> AnnData:
    """Run t-SNE on PCA embedding and plot colored by drug_condition."""
    tsne = TSNE(random_state=random_state)
    adata.obsm['X_tsne'] = tsne.fit_transform(adata.obsm['X_pca'])
    sc.pl.tsne(adata, color='drug_condition')
    plt.tight_layout()
    return adata


def plot_umap(adata: AnnData) -> AnnData:
    """Compute UMAP neighborhood graph and plot colored by drug_condition."""
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='drug_condition')
    plt.tight_layout()
    return adata


def draw_trajectory(input_data) -> None:
    """Draw a force-directed graph trajectory from an AnnData object or h5ad path."""
    import os
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80, frameon=False,
                                  figsize=(3, 3), facecolor='white')

    if isinstance(input_data, str):
        adata = sc.read_h5ad(input_data)
    elif isinstance(input_data, AnnData):
        adata = input_data
    else:
        raise ValueError("input_data must be a file path (str) or AnnData object.")

    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=10)
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata, color='drug_condition', legend_loc='on data')


# =============================================================================
# 3. QC / boxplots / correlation plots
# =============================================================================

def _get_high_variance_genes(adata: AnnData, top_n: int = 2):
    gene_vars = np.var(adata.X, axis=0)
    if sparse.issparse(adata.X):
        gene_vars = np.asarray(gene_vars).ravel()
    top_idx = np.argsort(gene_vars)[-top_n:]
    return adata.var.index[top_idx]


def plot_boxplots(adata: AnnData, top_n: int = 3) -> None:
    """Boxplots of normalized counts for the top-N highest-variance genes."""
    top_genes = _get_high_variance_genes(adata, top_n=top_n)
    for gene in top_genes:
        plot_df = pd.DataFrame({
            'file': adata.obs['file'],
            'drug_no': adata.obs['*drug no'],
            'normalized_count': adata[:, gene].X.flatten().astype(np.float64),
        })
        x, y = 'file', 'normalized_count'
        plt.figure(figsize=(12, 6))
        box = sns.boxplot(x=x, y=y, data=plot_df)
        plt.xticks(rotation=90)
        pairs = list(itertools.combinations(plot_df['file'].unique(), 2))
        ann = Annotator(box, pairs, data=plot_df, x=x, y=y)
        ann.configure(test='Mann-Whitney', text_format='star', loc='inside')
        ann.apply_and_annotate()
        plt.title(f'Normalized Counts — {gene}')
        plt.tight_layout()
        plt.show()


def draw_correlation_matrix(
    data: pd.DataFrame,
    title: str = 'Correlation Matrix',
    fontsize: int = 12,
    cmap: str = 'coolwarm',
) -> None:
    """Annotated heatmap of the pairwise Pearson correlation matrix."""
    corr = data.corr()
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap=cmap, cbar=True)
    plt.title(title, fontsize=fontsize)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()


def plot_correlation_differences(
    data1: pd.DataFrame,
    data2: pd.DataFrame,
    column1: str,
    column2: str,
) -> None:
    """Bar plot of (corr with column1) − (corr with column2) for each column in data1."""
    if data1.index.name != 'Unnamed: 0':
        data1 = data1.set_index('Unnamed: 0')
    if data2.index.name != 'Unnamed: 0':
        data2 = data2.set_index('Unnamed: 0')

    diff = (
        data1.apply(lambda x: x.corr(data2[column1]))
        - data1.apply(lambda x: x.corr(data2[column2]))
    ).dropna()

    plt.figure(figsize=(10, 6))
    diff.plot(kind='bar')
    plt.title(f'Δ Correlation ({column1} − {column2})')
    plt.ylabel('Difference in Pearson r')
    plt.xlabel('Columns in data1')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()


def plot_correlation_bar_plots(
    data1: pd.DataFrame,
    data2: pd.DataFrame,
    col1: str,
    col2: str,
) -> None:
    """Side-by-side bar plots of correlation coefficients against two reference columns."""
    r1 = data1.apply(lambda x: x.corr(data2[col1])).dropna()
    r2 = data1.apply(lambda x: x.corr(data2[col2])).dropna()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    r1.plot(kind='bar', ax=axes[0])
    axes[0].set_title(f"Correlations with '{col1}'")
    axes[0].set_ylabel('Pearson r')
    axes[0].set_xlabel('Columns in data1')
    axes[0].tick_params(axis='x', rotation=45)

    r2.plot(kind='bar', ax=axes[1], color='orange')
    axes[1].set_title(f"Correlations with '{col2}'")
    axes[1].set_ylabel('Pearson r')
    axes[1].set_xlabel('Columns in data1')
    axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.show()
