
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd


def plot_mappability_distribution(
    data: Dict[str, Dict[str, np.ndarray]], output_dir: str
) -> None:
    """
    Plot the distribution of mappability scores for different k-mer sizes.

    Args:
        data (Dict[str, Dict[str, np.ndarray]]): Mappability data for different k-mer sizes.
        output_dir (str): Directory to save the output plot.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    for kmer, chrom_data in data.items():
        all_values = np.concatenate(list(chrom_data.values()))
        sns.kdeplot(all_values, label=f"{kmer}-mers", ax=ax)

    ax.set_xlabel("Mappability Score")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Mappability Scores")
    ax.legend()
    plt.savefig(os.path.join(output_dir, "mappability_distribution.png"))
    plt.close()


def plot_mappability_changes(
    results: Dict[str, Dict[str, np.ndarray]], output_dir: str
) -> None:
    """
    Plot the changes in mappability scores between different k-mer sizes.

    Args:
        results (Dict[str, Dict[str, np.ndarray]]): Results of mappability comparisons.
        output_dir (str): Directory to save the output plots.
    """
    for comparison, diff_data in results.items():
        fig, ax = plt.subplots(figsize=(10, 6))
        all_diffs = np.concatenate(list(diff_data.values()))
        sns.histplot(all_diffs, kde=True, ax=ax)

        ax.set_xlabel("Mappability Score Difference")
        ax.set_ylabel("Count")
        ax.set_title(f"Mappability Changes: {comparison}")
        plt.savefig(
            os.path.join(output_dir,
                         f"mappability_changes_{comparison}.png")
                         )
        plt.close()

def plot_gene_mappability_changes(gene_mappability: Dict[str, pd.DataFrame], output_dir: str) -> None:
    """
    Plot the changes in gene mappability across different k-mer sizes.

    Args:
        gene_mappability (Dict[str, pd.DataFrame]): Gene mappability data for different k-mer sizes.
        output_dir (str): Directory to save the output plots.
    """
    kmer_sizes = sorted(gene_mappability.keys())
    merged_data = gene_mappability[kmer_sizes[0]].copy()
    merged_data.set_index('gene_id', inplace=True)
    for kmer in kmer_sizes[1:]:
        merged_data[f'mappability_ratio_{kmer}'] = gene_mappability[kmer].set_index('gene_id')['mappability_ratio']
    
    # rename kmer[0] to mappability_ratio_kmer[0]
    merged_data.rename(columns={'mappability_ratio': 'mappability_ratio_{}'.format(kmer_sizes[0])}, inplace=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=merged_data[[f'mappability_ratio_{kmer}' for kmer in kmer_sizes]], ax=ax)
    ax.set_xticklabels(kmer_sizes)
    ax.set_xlabel('K-mer size')
    ax.set_ylabel('Mappability ratio')
    ax.set_title('Distribution of Gene Mappability Ratios Across K-mer Sizes')
    plt.savefig(os.path.join(output_dir, 'gene_mappability_distribution.png'))
    plt.close()

    # Plot heatmap of top 100 genes with highest variance in mappability
    variance = merged_data[[f'mappability_ratio_{kmer}' for kmer in kmer_sizes]].var(axis=1)
    top_100_genes = variance.nlargest(100).index
    top_100_data = merged_data.loc[top_100_genes, [f'mappability_ratio_{kmer}' for kmer in kmer_sizes]]
    
    fig, ax = plt.subplots(figsize=(12, 20))
    sns.heatmap(top_100_data, cmap='YlOrRd', ax=ax)
    ax.set_xlabel('K-mer size')
    ax.set_ylabel('Genes')
    ax.set_title('Top 100 Genes with Highest Variance in Mappability')
    plt.savefig(os.path.join(output_dir, 'top_100_gene_mappability_heatmap.png'))
    plt.close()

    return top_100_data