"""
Mappability Analysis Tool

This script analyzes and compares mappability across different k-mer sizes
using BigWig files of mappability tracks.

It provides functions to load BigWig data, compare mappability between
different k-mer sizes, and visualize the results.
"""

import os
import argparse
from typing import Dict, List, Tuple
import numpy as np
import pyBigWig

from visualise import plot_mappability_distribution, plot_mappability_changes


def load_bigwig(file_path: str) -> Dict[str, np.ndarray]:
    """
    Load a BigWig file and extract mappability data.

    Args:
        file_path (str): Path to the BigWig file.

    Returns:
        Dict[str, np.ndarray]: A dictionary where keys are chromosome names
        and values are numpy arrays of mappability scores.
    """
    bw = pyBigWig.open(file_path)
    chrom_sizes = dict(bw.chroms().items())
    data = {}
    for chrom, size in chrom_sizes.items():
        data[chrom] = bw.values(chrom, 0, size)
    bw.close()
    return data


def compare_mappability(
    data1: Dict[str, np.ndarray], data2: Dict[str, np.ndarray]
) -> Dict[str, np.ndarray]:
    """
    Compare mappability data between two different k-mer sizes.

    Args:
        data1 (Dict[str, np.ndarray]): Mappability data for the first k-mer size.
        data2 (Dict[str, np.ndarray]): Mappability data for the second k-mer size.

    Returns:
        Dict[str, np.ndarray]: A dictionary of differences in mappability scores.
    """
    diff = {}
    for chrom in data1.keys():
        diff[chrom] = data2[chrom] - data1[chrom]
    return diff


def analyze_mappability_changes(
    bigwig_files: List[str],
) -> Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, Dict[str, np.ndarray]]]:
    """
    Analyze mappability changes across different k-mer sizes.

    Args:
        bigwig_files (List[str]): List of paths to BigWig files.

    Returns:
        Tuple[Dict[str, Dict[str, np.ndarray]], Dict[str, Dict[str, np.ndarray]]]:
        A tuple containing results of mappability comparisons and the original data.
    """
    data = {os.path.basename(f).split("_")[0]: load_bigwig(f) for f in bigwig_files}

    results = {}
    kmer_sizes = sorted(data.keys())

    for i in range(len(kmer_sizes) - 1):
        for j in range(i + 1, len(kmer_sizes)):
            k1, k2 = kmer_sizes[i], kmer_sizes[j]
            diff = compare_mappability(data[k1], data[k2])
            results[f"{k1}_vs_{k2}"] = diff

    return results, data


def main(args: argparse.Namespace) -> None:
    """
    Main function to run the mappability analysis.

    Args:
        args (argparse.Namespace): Command-line arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    bigwig_files = [
        os.path.join(args.bigwig_dir, f)
        for f in os.listdir(args.bigwig_dir)
        if f.endswith(".bw")
    ]

    results, data = analyze_mappability_changes(bigwig_files)

    plot_mappability_distribution(data, args.output_dir)
    plot_mappability_changes(results, args.output_dir)

    if args.verbose:
        print(f"Analysis complete. Output saved to {args.output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze and compare mappability across different k-mer sizes."
    )
    parser.add_argument(
        "bigwig_dir", help="Directory containing BigWig files of mappability tracks"
    )
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Increase output verbosity"
    )
    args = parser.parse_args()

    main(args)
