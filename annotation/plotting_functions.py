import numpy as np
import pandas as pd
import scanpy as sc
import decoupler as dc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from itertools import combinations
import anndata
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import random

def get_transcript_lengths(gtf_file_path):
    transcript_lengths = {}  # Dictionary to store transcript lengths
    with open(gtf_file_path) as gtf_file:
        for line in gtf_file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'transcript':
                attributes = parts[8].split('; ')
                transcript_id = attributes[1].split('"')[1]
                gene_id = attributes[0].split('"')[1]
                transcript_length = int(parts[4]) - int(parts[3]) + 1
                if transcript_id not in transcript_lengths:
                    transcript_lengths[transcript_id] = transcript_length
    return transcript_lengths

def plot_genes_by_counts(adata, category_column='plate', figsize=(10, 7)):
    adata_obs_df = adata.obs
    adata_obs_df.reset_index(drop=True, inplace=True)

    # Use the 'husl' color palette
    unique_categories = sorted(adata_obs_df[category_column].unique())
    plate_palette = sns.color_palette("husl", n_colors=len(unique_categories))

    plt.figure(figsize=(10, 7))

    sns.scatterplot(x=adata_obs_df['total_counts'], y=adata_obs_df['n_genes_by_counts'],
                    hue=adata_obs_df[category_column], alpha=0.5, palette=plate_palette, s=3)

    # Add labels and title
    plt.xlabel('UMI Counts', fontsize=16)
    plt.ylabel('Expressed genes', fontsize=16)
    
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(fontsize=14)

    plt.grid(True, which="both")


def plot_knee(adata, cutoff, category_column='plate', figsize=(10, 7)):
    unique_categories = sorted(adata.obs[category_column].unique())
    color_palette = sns.color_palette("husl", n_colors=len(unique_categories))

    fig, ax = plt.subplots(figsize=figsize)

    for category, color in zip(unique_categories, color_palette):
        # Subset the data for each plate
        plate_subset = adata[adata.obs[category_column] == category]

        # Calculate knee and num_cells for the subset
        knee = np.sort((np.array(plate_subset.layers['raw_counts'].sum(axis=1))).flatten())[::-1]
        cell_set = np.arange(len(knee))
        num_cells = cell_set[::-1][0]

        # Plot the knee for the current plate
        ax.loglog(cell_set, knee, linewidth=5, color=color, label=f"{category} - {num_cells} nuclei")

    ax.axhline(y=cutoff, linewidth=3, color="k", linestyle='--')

    ax.set_ylabel("UMI Counts",fontsize=18)
    ax.set_xlabel("Set of Barcodes",fontsize=18)
    ax.legend(fontsize=14)
    
    ax.tick_params(axis='both', which='major', labelsize=16)


    plt.grid(True, which="both")

def stacked_barplot_proportions(adata, cluster_key, var_key, flip=True, fsize = (12,6), annotations = True):
    m=adata.groupby([cluster_key]).size().to_frame()
    # Group the data by 'cluster_key' and 'var_key', count occurrences, and calculate proportions
    grouped_data = adata.groupby([cluster_key, var_key]).size().unstack().fillna(0)
    proportions = grouped_data.div(grouped_data.sum(axis=1), axis=0)
    
    if var_key == "Genotype":
        color_palette = ['#DA9CC1', '#F4C245', '#C0BFBF', '#55AF5B', '#4F6EAF', '#52A5DB', '#D83026', '#683C91']
        
    elif var_key == "subtype":
        subtype_counts = adata['subtype'].value_counts()
        celltype_sizes = pd.Series(subtype_counts)
        sorted_celltypes = celltype_sizes.sort_values(ascending=False).index
        proportions.set_index('leiden', inplace=True)
        
        # rndom color palette for now
        unique_values = adata[var_key].unique()
        n_unique_values = len(unique_values)
        color_palette = sns.color_palette("husl", n_unique_values)
              
    elif var_key == "Sex" and cluster_key == "leiden":
        color_palette = ['hotpink','dodgerblue']
        
    else:
        unique_values = adata[var_key].unique()
        n_unique_values = len(unique_values)
        color_palette = sns.color_palette("husl", n_unique_values)
    
    # Create the stacked bar plot
    if flip:
        if cluster_key == "leiden":
            proportions = proportions.sort_values(by='leiden', ascending=False)
            m = m.sort_values(by='leiden', ascending=False)
             
        if annotations:
            ax = proportions.plot(kind='barh',color=color_palette, 
                                  stacked=True, figsize=fsize, width = 0.8)
            l = len(adata[var_key].dropna().unique())-1 
            ax.bar_label(ax.containers[l], labels=m[0], label_type='edge', padding = 10, fontsize = 14)
            
            plt.xlim(0,1.15)
            ax.tick_params(axis="x", labelsize = 14)
            ax.tick_params(axis="y", labelsize = 14)
            ax.set_xlabel("Proportion")
            ax.set_ylabel(cluster_key)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            ax.legend(title=var_key, bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(False)
            
            
        else:
            ax = proportions.plot(kind='barh',color=color_palette, 
                                  stacked=True, figsize=fsize, width = 0.8)
            plt.xlim(0,1.15)
            
            ax.tick_params(axis="x", labelsize = 14)
            ax.tick_params(axis="y", labelsize = 14)
            ax.set_xlabel("Proportion")
            ax.set_ylabel(cluster_key)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.get_legend().remove()
            ax.grid(False)

    else:
        if annotations:
            ax = proportions.plot(kind='bar',color=color_palette, 
                                  stacked=True, figsize=fsize, width = 0.8)
            l = len(adata[var_key].dropna().unique())-1 
            ax.bar_label(ax.containers[l], labels=m[0], label_type='edge', padding = 10, fontsize = 14, rotation=90)
            plt.ylim(0,1.15)
            ax.tick_params(axis="x", labelsize = 14, rotation = 40)
            ax.tick_params(axis="y", labelsize = 14)
            ax.set_xlabel(cluster_key)
            ax.set_ylabel("Proportion")
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Show the legend outside the plot
            ax.legend(title=var_key, bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(False)
        
        else:
            ax = proportions.plot(kind='bar',color=color_palette, 
                                  stacked=True, figsize=fsize, width = 0.8)
            plt.ylim(0,1.15)
            ax.tick_params(axis="x", labelsize = 14, rotation = 40)
            ax.tick_params(axis="y", labelsize = 14)
            ax.set_xlabel(cluster_key)
            ax.set_ylabel("Proportion")
            ax.set_title(f'{var_key} by {cluster_key}')
        

            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

            # Show the legend outside the plot
            ax.legend(title=var_key, bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(False)
        

    plt.show()
    
    
def calc_total_counts(adata, obs_col='dataset', layer='counts'):
    # turn into a sparse dataframe
    cols = adata.var.index.tolist()
    inds = adata.obs[obs_col].tolist()
    data = adata.layers[layer].toarray()
    df = pd.DataFrame(data=data, index=inds, columns=cols)
    df.index.name = obs_col

    # add up values on condition (row)
    df = df.groupby(level=0).sum()

    return df  
