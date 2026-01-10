# GE-HTS project pipeline for in situ sequencing data analysis
# author: Nathan Wooseok Lee
# Date of start: 230801
# Date of update: 231214
# env: conda activate gehts

from matplotlib import axes
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
import os
from glob import glob
from tkinter import filedialog
# from tkinter import Tk
from anndata import AnnData
from scipy import sparse
from scipy.stats import kruskal
from statannotations.Annotator import Annotator
import itertools
import re
import scanpy as sc

def process_data(save=None):

    anabolic = ['Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Il17b','Pth1r','Gm26633','Col27a1']
    inflammatory= ['Mmp3','Mmp13','Il6','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b']
    Housekeeping = ['Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23']



    lookuptable1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14] # dose 1 um, 14 is Xav
    lookuptable01 = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 12] # dose 0.1 um, 12 is Xav
    lookuptable10 = [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 3] # dose 10 um, 3 is Xav

    drug_name = ['rapamycin', 'SB431542', 'SANT-1', 'LY294002', 'SB525334', 'KU0063794', 'ALK5 inhibitor IV', 'MK-2206 dihydrochloride', 'SB203580', 'pamapimod', 'BMS-345541', 'JNK inhibitor V', 'CAPE', 'XAV']    
    
    # Create individual lookup dictionaries
    lookup_dict1 = {num: drug_name[i] + '_1' for i, num in enumerate(lookuptable1)}
    lookup_dict01 = {num: drug_name[i] + '_0.1' for i, num in enumerate(lookuptable01)}
    lookup_dict10 = {num: drug_name[i] + '_10' for i, num in enumerate(lookuptable10)}
    
    # Function to select the appropriate lookup table based on filename
    def __get_lookup_table__(filename):
        category = filename[7:9]
        if category == '01':
            return lookup_dict1
        elif category == '10':
            return lookup_dict10
        elif category == '-1':
            return lookup_dict01
        else:
            return None
        

    # Function to replace drug numbers with drug names based on the correct lookup table
    def __replace_drug_numbers__(row, col):
        lookup_table = __get_lookup_table__(row['file'])
        if lookup_table:
            return lookup_table.get(row[col], row[col])
        else:
            return row[col]


    # Create a tkinter root and hide it
    # root = Tk()
    # root.withdraw()

    # Open a file dialog and get the directory path
    dir_path = filedialog.askdirectory()

    # Gather all CSV files in the directory
    files = glob(os.path.join(dir_path, '*.csv'))

    df_list = []
    for file in files:
        df = pd.read_csv(file)
        df['file'] = os.path.basename(file)
        df_list.append(df)
        
    df = pd.concat(df_list, ignore_index=True)
 
    # Ensure '*drug no' column exists
    if '*drug no' not in df.columns:
        df['*drug no'] = np.nan

    # List of potential columns
    potential_columns = ['file', 'GeneName', '*drug no', '*drug1 no', '*drug2 no']

    # Determine which of these columns actually exist in the dataframe
    existing_columns = [col for col in potential_columns if col in df.columns]

    # Now create AnnData object using only the existing columns
    adata = AnnData(X=df.loc[:, ~df.columns.isin(existing_columns)], 
                    obs=df[existing_columns])


    # Apply the replacement for '*drug no', '*drug1 no', and '*drug2 no'
    drug_columns = ['*drug no', '*drug1 no', '*drug2 no']
    for col in drug_columns:
        if col in adata.obs.columns:
            # Apply the function row-wise
            adata.obs[col] = adata.obs.apply(lambda row: __replace_drug_numbers__(row, col), axis=1)


    # Check if '*drug1 no' and '*drug2 no' exist and create a combined column
    if '*drug1 no' in adata.obs.columns and '*drug2 no' in adata.obs.columns:
        combined_drug = adata.obs['*drug1 no'].astype(str) + '&' + adata.obs['*drug2 no'].astype(str)
        # Use the combined column to update '*drug no' where it's null
        adata.obs.loc[adata.obs['*drug no'].isnull(), '*drug no'] = combined_drug

        adata.obs['*drug1 no'] = adata.obs['*drug1 no'].astype(str)
        adata.obs['*drug2 no'] = adata.obs['*drug2 no'].astype(str)

    # Convert '*drug no' to a categorical column
    adata.obs['*drug no'] = adata.obs['*drug no'].astype(str) # before conversion to category, make uniform of mix of int and str.  
    adata.obs['*drug no'] = adata.obs['*drug no'].astype('category')

    # Creating the special sub-cases for '*drug no' == 0 for control and treatment control cases
    adata.obs['drug_no_subcase'] = np.where(
        (adata.obs['*drug no'] == '0') & adata.obs['file'].str.lower().str.contains('c\.csv$|control\.csv$', regex=True),
        'control',
        np.where(
            (adata.obs['*drug no'] == '0'),
            'inflammatory',
            ''  # We use empty string for all other cases
        )
    )
    # Creating a combined 'drug_condition' column
    adata.obs['drug_condition'] = np.where(
        adata.obs['drug_no_subcase'] == '',
        adata.obs['*drug no'].astype(str),  # Normal '*drug no' categories
        adata.obs['drug_no_subcase']  # Special sub-cases
    )    

    adata.obs_name = adata.obs['GeneName']
    adata.obs.rename(columns={"GeneName": "well_no"}, inplace=True)
    # adata.obs['drug_no'] = adata[:, '*drug no']
    # adata.obs['file'] = adata[:, 'file']
    # adata.obs_names = adata[:, "GeneName"]

    # remove NNNN and homomer from var
    mask = ~adata.var.index.isin(['NNNN', 'Homomer'])
    adata = adata[:, mask]

    # 1: anabolic, 2: inflammatory, 3: housekeeping, 0: others
    adata.var['gene_idx'] = np.where(adata.var.index.isin(anabolic), 1, 
                                        np.where(adata.var.index.isin(inflammatory), 2, 
                                                    np.where(adata.var.index.isin(Housekeeping), 3, 0)))
    
    adata.obs['id']=adata.obs['file'].astype(str)+adata.obs['well_no'].astype(str)

    if save == True:
        adata.write_h5ad(os.path.join(dir_path, 'merged.h5ad'))
    return adata


def normalize_by_ubc(adata):
    ubc_index = adata.var.index.get_loc('Ubc')
    ubc_counts = adata.X[:, ubc_index]

    # Identify cells with non-zero Ubc counts
    nonzero_ubc_cells = ubc_counts != 0

    if sparse.issparse(adata.X):
        nonzero_ubc_cells = np.asarray(nonzero_ubc_cells).ravel()

    # Subset the AnnData object to include only cells with non-zero Ubc counts
    adata = adata[nonzero_ubc_cells]

    # Normalize by Ubc counts
    adata.X = adata.X.astype(np.float64)
    adata.X = adata.X/adata.X[:, ubc_index][:, None]

    return adata

#def RECscore(adata):
    

def draw_heatmap(adata):
    # Sort var by 'gene_idx' and then by index (gene names)
    adata.var.sort_values(['gene_idx'], inplace=True)

    # Sort X by the new order of var
    sorted_order = np.argsort(adata.var.index)
    adata.X = adata.X[:, sorted_order]
    
    # Convert the count matrix to log scale
    log_counts = np.log2(adata.X + 1)
        
    # Draw the heatmap
    plt.figure(figsize=(10,10))
    sns.heatmap(log_counts, cmap='viridis', xticklabels=adata.obs.index, yticklabels=adata.var.index)
    plt.title('Heatmap of Normalized Gene Counts')
    plt.show()

def plot_boxplots(adata):
    # Find the genes with the highest variance
    top_genes = get_high_variance_genes(adata, top_n=3)
    
    for gene in top_genes:
        
        # Create a DataFrame for plotting
        plot_df = pd.DataFrame({
            'file': adata.obs['file'],
            'drug_no': adata.obs['*drug no'],
            'normalized_count': adata[:,gene].X.flatten().astype(np.float64)
        })
        # print(adata[:,gene].X)
        # print(type(adata[:,gene].X))
        x='file'
        y='normalized_count'
        print(plot_df['normalized_count'])
        # Create boxplot
        plt.figure(figsize=(12, 6))
        box = sns.boxplot(x=x, y=y, data=plot_df, hue=None)
        plt.xticks(rotation=90)

        # Calculate and annotate p-values
        # groups = [plot_df.loc[(plot_df['file'] == group), 'normalized_count'] for group in plot_df['file'].unique()]
        # _, p_val = kruskal(*groups)

        file_list = plot_df['file'].unique()
        pairs = list(itertools.combinations(file_list,2))
        annotator = Annotator(box, pairs, data=plot_df,x=x,y=y)
        annotator.configure(test='Mann-Whitney', text_format='star', loc='inside')
        annotator.apply_and_annotate()
        # plt.text(len(groups) / 2, max(group.max() for group in groups), f'p={p_val:.3f}', horizontalalignment='center')
        
        plt.title(f'Normalized Counts Boxplot for {gene}')
        plt.tight_layout()
        plt.show()

def get_high_variance_genes(adata, top_n=2):
    # Compute the variance for each gene
    gene_vars = np.var(adata.X, axis=0)
    if sparse.issparse(adata.X):
        gene_vars = np.asarray(gene_vars).ravel()

    # Get the indices of the top_n genes with the highest variance
    top_gene_indices = np.argsort(gene_vars)[-top_n:]

    return adata.var.index[top_gene_indices]


def plot_pca(adata):
    # Perform PCA on the data.
    # You might want to normalize and scale your data before performing PCA
    # First, normalize and log-transform the data
    sc.pp.normalize_total(adata, target_sum=1e3)
    sc.pp.log1p(adata)
    # Scale the datapoints for PCA, which is presumption of PCA method.
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)

    # Plot the results
    sc.pl.pca(adata, color='drug_condition')
    plt.tight_layout()

    return adata

def plot_tsne(adata, random_state=0):

    # Perform t-SNE on the PCA results
    tsne = TSNE(random_state=random_state)
    tsne_results = tsne.fit_transform(adata.obsm['X_pca'])

    # Add t-SNE results to the AnnData object
    adata.obsm['X_tsne'] = tsne_results

    # Plot the t-SNE results
    sc.pl.tsne(adata, color='drug_condition')
    plt.tight_layout()

    return adata


def plot_umap(adata):
    # Compute the UMAP
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)

    # Plot UMAP
    sc.pl.umap(adata, color='drug_condition')
    plt.tight_layout()

    return adata


def draw_trajectory(input_data):
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()

    if isinstance(input_data, str):
        adata = sc.read_h5ad(input_data)
        file_name = os.path.basename(input_data)
        # results_file = os.path.join('./trajectory', file_name)

    elif isinstance(input_data, anndata.AnnData):
        adata = input_data

    else:
        raise ValueError("Input should be either a string file address or an AnnData object.")

    sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures
    sc.pp.neighbors(adata, n_neighbors=4, n_pcs=10)
    sc.tl.draw_graph(adata)
    sc.pl.draw_graph(adata, color='paul15_clusters', legend_loc='on data')

anabolic_genes = ['Acan','Sox9','Col2a1','Matn1','Matn3','Ucma','Ccnd3','Gadd45g','Pth1r','Gm26633','Col27a1']
inflammatory_genes = ['Mmp3','Mmp13','Il6','Il17b','Adamts5','Igfbp3','Ccl2','Cxcl5','Cxcl1','Fosl2','Tlr2','Tnfrsf1b']
housekeeping_genes = ['Hprt','Actb','Gapdh','B2m','Ubc','Ppia','Rpl23']

def draw_heatmap_separated(data, anabolic_genes=anabolic_genes, inflammatory_genes=inflammatory_genes, housekeeping_genes=housekeeping_genes):
    # Set index to gene names if not already set
    if data.index.name != 'Unnamed: 0':
        data.set_index('Unnamed: 0', inplace=True)
    # Filter data for each gene category
    anabolic_data = data.loc[anabolic_genes]
    inflammatory_data = data.loc[inflammatory_genes]
    housekeeping_data = data.loc[housekeeping_genes]
    # Log transform the data
    log_anabolic_data = np.log(anabolic_data + 1)
    log_inflammatory_data = np.log(inflammatory_data + 1)
    log_housekeeping_data = np.log(housekeeping_data + 1)
    # Plotting the heatmaps
    plt.figure(figsize=(8, 18))
    plt.rcParams.update({'font.size': 8})
    # Anabolic heatmap (Linear)
    plt.subplot(6, 1, 1)
    sns.heatmap(anabolic_data, cmap="viridis")
    plt.title("Anabolic Genes Expression Heatmap (Linear Scale)")
    plt.ylabel("Anabolic Genes")
    # Anabolic heatmap (Log)
    plt.subplot(6, 1, 2)
    sns.heatmap(log_anabolic_data, cmap="viridis")
    plt.title("Anabolic Genes Expression Heatmap (Log Scale)")
    plt.ylabel("Anabolic Genes")
    # Inflammatory heatmap (Linear)
    plt.subplot(6, 1, 3)
    sns.heatmap(inflammatory_data, cmap="viridis")
    plt.title("Inflammatory Genes Expression Heatmap (Linear Scale)")
    plt.ylabel("Inflammatory Genes")
    # Inflammatory heatmap (Log)
    plt.subplot(6, 1, 4)
    sns.heatmap(log_inflammatory_data, cmap="viridis")
    plt.title("Inflammatory Genes Expression Heatmap (Log Scale)")
    plt.ylabel("Inflammatory Genes")
    # Housekeeping heatmap (Linear)
    plt.subplot(6, 1, 5)
    sns.heatmap(housekeeping_data, cmap="viridis")
    plt.title("Housekeeping Genes Expression Heatmap (Linear Scale)")
    plt.ylabel("Housekeeping Genes")
    # Housekeeping heatmap (Log)
    plt.subplot(6, 1, 6)
    sns.heatmap(log_housekeeping_data, cmap="viridis")
    plt.title("Housekeeping Genes Expression Heatmap (Log Scale)")
    plt.ylabel("Housekeeping Genes")
    plt.tight_layout(pad=3.0)
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    """
    Draws heatmaps for gene expression data.

    Parameters:
    data (DataFrame): The gene expression dataset.
    anabolic_genes (list): List of anabolic genes.
    inflammatory_genes (list): List of inflammatory genes.
    housekeeping_genes (list): List of housekeeping genes.
    """

def draw_heatmap_separated_lin(data, anabolic_genes=anabolic_genes, inflammatory_genes=inflammatory_genes, housekeeping_genes=housekeeping_genes):
    # Set index to gene names if not already set
    if data.index.name != 'Unnamed: 0':
        data.set_index('Unnamed: 0', inplace=True)
    # Filter data for each gene category
    anabolic_data = data.loc[anabolic_genes]
    inflammatory_data = data.loc[inflammatory_genes]
    housekeeping_data = data.loc[housekeeping_genes]
    # Calculate the total height of the figure
    total_height = (len(anabolic_data) + len(inflammatory_data) + len(housekeeping_data))
    plt.figure(figsize=(8, total_height))
    # Create a figure with adjusted height
    grid_spec = plt.GridSpec(3, 1, height_ratios=[len(anabolic_data), len(inflammatory_data), len(housekeeping_data)])
    # Anabolic heatmap (Linear)
    ax1 = plt.subplot(grid_spec[0])
    sns.heatmap(anabolic_data, cmap="viridis")
    # plt.title("Anabolic Genes Expression Heatmap (Linear Scale)")
    ax1.yaxis.set_label_position("right")
    plt.ylabel("Anabolic Genes")
    plt.yticks(rotation=0)
    plt.xticks([])    
    # Inflammatory heatmap (Linear)
    ax2 = plt.subplot(grid_spec[1])
    sns.heatmap(inflammatory_data, cmap="viridis")
    # plt.title("Inflammatory Genes Expression Heatmap (Linear Scale)")
    ax2.yaxis.set_label_position("right")
    plt.ylabel("Catabolic Genes")
    plt.yticks(rotation=0)
    plt.xticks([])
    # Housekeeping heatmap (Linear)
    ax3 = plt.subplot(grid_spec[2])
    sns.heatmap(housekeeping_data, cmap="viridis")
    # plt.title("Housekeeping Genes Expression Heatmap (Linear Scale)")
    ax3.yaxis.set_label_position("right")
    plt.ylabel("Housekeeping Genes")
    plt.yticks(rotation=0)
    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(hspace=1)
    plt.show()


def draw_heatmap_separated_wohk(data, anabolic_genes=anabolic_genes, inflammatory_genes=inflammatory_genes):
    # Set index to gene names if not already set
    if data.index.name != 'Unnamed: 0':
        data.set_index('Unnamed: 0', inplace=True)
    # Filter data for each gene category
    anabolic_data = data.loc[anabolic_genes]
    inflammatory_data = data.loc[inflammatory_genes]
    # Plotting the heatmaps
    plt.figure(figsize=(9, 12))
    plt.rcParams.update({'font.size': 8})
    # Anabolic heatmap (Linear)
    plt.subplot(2, 1, 1)
    sns.heatmap(anabolic_data, cmap="viridis")
    plt.title("Anabolic Genes Expression Heatmap (Linear Scale)")
    plt.ylabel("Anabolic Genes")
    plt.xticks([])
    # Inflammatory heatmap (Linear)
    plt.subplot(2, 1, 2)
    sns.heatmap(inflammatory_data, cmap="viridis")
    plt.title("Inflammatory Genes Expression Heatmap (Linear Scale)")
    plt.ylabel("Inflammatory Genes")
    plt.xticks(rotation=45, ha='right') # change if you handle few sample points.
    # Housekeeping heatmap (Linear)
    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(hspace=1)
    plt.show()

def draw_heatmap_lin(data):
    # Set index to gene names if not already set
    if data.index.name != 'Unnamed: 0':
        data.set_index('Unnamed: 0', inplace=True)
    # Plotting the heatmaps
    plt.figure(figsize=(8, 9))
    plt.rcParams.update({'font.size': 8})
    # Anabolic heatmap (Linear)
    sns.heatmap(data, cmap="viridis")
    plt.title("Genes Expression Heatmap (Linear Scale)")
    # plt.xticks([])
    plt.tight_layout(pad=3.0)
    plt.subplots_adjust(hspace=0.5)
    plt.show()

def draw_correlation_matrix(data, title="Correlation Matrix", fontsize=12, cmap='coolwarm'):
    # Compute the correlation matrix
    correlation_matrix = data.corr()
    # Plotting the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap=cmap, cbar=True)
    plt.title(title, fontsize=fontsize)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    """
    Draws a heatmap for the correlation matrix of the given data.

    Parameters:
    data (DataFrame): The dataset for which the correlation matrix will be computed and plotted.
    title (str): Title for the heatmap.
    fontsize (int): Font size for the heatmap labels and title.
    cmap (str): Colormap for the heatmap.
    """


def plot_correlation_differences(data1, data2, column1, column2):
    if data1.index.name != 'Unnamed: 0':
        data1.set_index('Unnamed: 0', inplace=True)
        
    if data2.index.name != 'Unnamed: 0':
        data2.set_index('Unnamed: 0', inplace=True)
    
    correlations_column1 = data1.apply(lambda x: x.corr(data2[column1]))
    correlations_column2 = data1.apply(lambda x: x.corr(data2[column2]))
    # Calculate the differences in correlations
    correlation_differences = correlations_column1 - correlations_column2
    # Filter out NaN values
    correlation_differences_filtered = correlation_differences.dropna()
    # Plotting the difference bar plot
    plt.figure(figsize=(10, 6))
    correlation_differences_filtered.plot(kind='bar')
    plt.title(f"Difference in Correlation Coefficients ({column1} - {column2})")
    plt.ylabel('Difference in Correlation Coefficient')
    plt.xlabel('Columns in Data1')
    plt.xticks(rotation=45, ha='right')
    plt.show()
    """
    Plots the difference in correlation coefficients (column1 - column2) for each column in data1 with two columns in data2.

    Parameters:
    - data1: DataFrame containing multiple columns to correlate.
    - data2: DataFrame containing the columns to correlate against.
    - column1: The name of the first column in data2 to correlate against.
    - column2: The name of the second column in data2 to correlate against.
    """

def plot_correlation_bar_plots(data1, data2, col1, col2):
    # Extract the specified columns from data2
    column1 = data2[col1]
    column2 = data2[col2]
    # Calculate the correlations for each column in data1 with the specified columns in data2
    correlations_col1 = data1.apply(lambda x: x.corr(column1))
    correlations_col2 = data1.apply(lambda x: x.corr(column2))
    # Filter out NaN and no correlation (zero correlation) values
    correlations_col1_filtered = correlations_col1.dropna()
    correlations_col2_filtered = correlations_col2.dropna()
    # Plotting the bar plots
    plt.figure(figsize=(14, 6))
    # Bar plot for the first specified column
    plt.subplot(1, 2, 1)
    correlations_col1_filtered.plot(kind='bar')
    plt.title(f"Correlations with '{col1}'")
    plt.ylabel('Correlation Coefficient')
    plt.xlabel('Columns in Data1 Dataset')
    plt.xticks(rotation=45, ha='right')
    # Bar plot for the second specified column
    plt.subplot(1, 2, 2)
    correlations_col2_filtered.plot(kind='bar', color='orange')
    plt.title(f"Correlations with '{col2}'")
    plt.ylabel('Correlation Coefficient')
    plt.xlabel('Columns in Data1 Dataset')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    """
    Plots bar plots of correlation coefficients between columns of two datasets.

    Parameters:
    - data1: DataFrame of the first dataset (for correlation calculation).
    - data2: DataFrame of the second dataset (to correlate against).
    - col1: The column name in data2 for the first correlation comparison.
    - col2: The column name in data2 for the second correlation comparison.
    """