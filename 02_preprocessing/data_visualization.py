# data_visualization.py

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
import numpy as np
import pandas as pd
from anndata import AnnData

def pd_heatmap(df):
    # Heatmap
    plt.figure(figsize=(10,10))
    gene_names = df['gene'].unique()
    sample_names = df['file'].unique()
    heatmap_data = df.pivot_table(index='gene', columns='file', values='normalized_count')
    sns.heatmap(np.log10(heatmap_data+1), cmap='viridis', cbar_kws={'label': 'Log10 Normalized Count'})
    plt.title('Heatmap of Normalized Gene Counts')
    plt.show()

def adata_heatmap(adata):
    # Convert the normalized count matrix into a DataFrame
    df = pd.DataFrame(data=adata.X,
                      index=adata.obs.index,
                      columns=adata.var.index)
    
    # Convert to log scale
    df_log = np.log10(df+1)

    # Draw the heatmap
    plt.figure(figsize=(10,10))
    sns.heatmap(df_log, cmap='viridis')
    plt.title('Heatmap of Normalized Gene Counts')
    plt.show()

def plot_tsne(df):
    # t-SNE plot
    X = df['normalized_count'].values
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(X.reshape(-1, 1))

    plt.figure(figsize=(10,10))
    sns.scatterplot(x=X_tsne[:,0], y=X_tsne[:,1], hue=df['gene_group'])
    plt.title('t-SNE Plot')
    plt.show()
