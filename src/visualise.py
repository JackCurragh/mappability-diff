
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os


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
