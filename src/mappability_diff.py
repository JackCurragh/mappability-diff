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
import gffpandas.gffpandas as gffpd
import pandas as pd


from visualise import plot_mappability_distribution, plot_mappability_changes, plot_gene_mappability_changes


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
        diff[chrom] = np.subtract(data2[chrom], data1[chrom])
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

def load_gtf(gtf_file: str):
    """
    Load a GTF file using gffpandas.

    Args:
        gtf_file (str): Path to the GTF file.

    Returns:
        gffpandas.GffPandas: A GffPandas object containing the GTF data.
    """
    return gffpd.read_gff3(gtf_file)

def calculate_gene_mappability(data: Dict[str, np.ndarray], gtf_data) -> pd.DataFrame:
    """
    Calculate mappability for each gene based on exon regions.

    Args:
        data (Dict[str, np.ndarray]): Mappability data for a specific k-mer size.
        gtf_data (gffpandas.GffPandas): GTF data loaded with gffpandas.

    Returns:
        pd.DataFrame: A dataframe with gene mappability information.
    """
    exons = gtf_data.filter_feature_of_type(['exon'])
    gene_mappability = []

    for _, exon in exons.df.iterrows():
        chrom = exon['seq_id']
        start = exon['start'] - 1  # Convert to 0-based
        end = exon['end']
        if 'gene_id' not in exon['attributes']:
            continue

        gene_id = next(i.split(' ')[1] for i in exon['attributes'].split(';') if 'gene_id' in i)

        if chrom in data:
            exon_map = data[chrom][start:end]
            total_bases = len(exon_map)
            mappable_bases = np.sum(np.array(exon_map) > 0.9)

            gene_mappability.append({
                'gene_id': gene_id,
                'total_exon_bases': total_bases,
                'mappable_bases': mappable_bases,
            })

    df = pd.DataFrame(gene_mappability)
    df = df.groupby(['gene_id']).sum().reset_index()
    df['mappability_ratio'] = df['mappable_bases'] / df['total_exon_bases']
    return df

def analyze_gene_mappability_changes(data: Dict[str, Dict[str, np.ndarray]], gtf_data) -> Dict[str, pd.DataFrame]:
    """
    Analyze mappability changes at the gene level across different k-mer sizes.

    Args:
        data (Dict[str, Dict[str, np.ndarray]]): Mappability data for different k-mer sizes.
        gtf_data (gffpandas.GffPandas): GTF data loaded with gffpandas.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary of dataframes with gene mappability information for each k-mer size.
    """
    gene_mappability = {}
    for kmer, mappability_data in data.items():
        gene_mappability[kmer] = calculate_gene_mappability(mappability_data, gtf_data)
    return gene_mappability


def get_gene_specific_mappability(data: Dict[str, np.ndarray], gtf_data) -> pd.DataFrame:
    """
    Return mappability for specific gene list on exon regions.

    Args:
        data (Dict[str, np.ndarray]): Mappability data for a specific k-mer size.
        gtf_data (gffpandas.GffPandas): GTF data loaded with gffpandas.

    Returns:
        pd.DataFrame: A dataframe with gene mappability information.
    """
    exons = gtf_data.filter_feature_of_type(['exon'])
    gene_mappability = []

    for _, exon in exons.df.iterrows():
        chrom = exon['seq_id']
        start = exon['start'] - 1  # Convert to 0-based
        end = exon['end']
        if 'gene_id' not in exon['attributes']:
            continue

        gene_id = next(i.split(' ')[1] for i in exon['attributes'].split(';') if 'gene_id' in i)

        if chrom in data:
            exon_map = data[chrom][start:end]


def main(args: argparse.Namespace) -> None:
    os.makedirs(args.output_dir, exist_ok=True)
    bigwig_files = [os.path.join(args.bigwig_dir, f) for f in os.listdir(args.bigwig_dir) if f.endswith('.bw')]
    
    results, data = analyze_mappability_changes(bigwig_files)
    
    plot_mappability_distribution(data, args.output_dir)
    plot_mappability_changes(results, args.output_dir)
    
    if args.gtf_file:
        gtf_data = load_gtf(args.gtf_file)
        gene_mappability = analyze_gene_mappability_changes(data, gtf_data)
        
        # Save gene mappability data
        for kmer, df in gene_mappability.items():
            df.to_csv(os.path.join(args.output_dir, f'gene_mappability_{kmer}.csv'), index=False)
        most_variable = plot_gene_mappability_changes(gene_mappability, args.output_dir)

    if args.verbose:
        print(f"Analysis complete. Output saved to {args.output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze and compare mappability across different k-mer sizes.")
    parser.add_argument("bigwig_dir", help="Directory containing BigWig files of mappability tracks")
    parser.add_argument("output_dir", help="Directory to save output files")
    parser.add_argument("--gtf_file", help="Path to GTF file for gene-level analysis")
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase output verbosity")
    args = parser.parse_args()

    main(args)