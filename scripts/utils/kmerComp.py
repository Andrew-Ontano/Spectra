#!/usr/bin/env python3
import argparse
import math
import random
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from typing import List, Tuple

# CLI arguments
parser = argparse.ArgumentParser(description="Streamed comparison of k-mer abundance between raw and assembly")
parser.add_argument("-i", "--input_ranked", type=str, help="Ranked TSV file from kmerRank.py")
parser.add_argument("-r", "--raw_dump", type=str, help="Raw jellyfish dump file")
parser.add_argument("-a", "--asm_dump", type=str, help="Assembly jellyfish dump file")
parser.add_argument("-k", "--kmer_size", type=int, default=20, help="K-mer size [default 20]")
parser.add_argument("-o", "--output_prefix", type=str, default="kmerComp_output", help="Output prefix")
parser.add_argument("-f", "--output_format", type=str, default="png", help="Output image format [default png]")
parser.add_argument("-s", "--plot_sample", type=int, default=1000000, help="Number of kmers to sample for plots")
parser.add_argument("-p", "--percentile", type=int, default=1, help="Percentile cutoff for extreme kmers")
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose mode', default=False)

args = parser.parse_args()

# Logging
logging.basicConfig(level=logging.ERROR)
logger = logging.getLogger()
if args.verbose:
    logger.setLevel(logging.INFO)


# Streaming merge
def stream_merge(raw_file: str, asm_file: str, sample_size: int) -> Tuple[List[Tuple[str, int, int]], int]:
    """
    Stream through two sorted dump files, yield matched kmers, and keep a reservoir sample.
    """
    sample = []
    total = 0

    with open(raw_file) as fr, open(asm_file) as fa:
        r = fr.readline().split()
        a = fa.readline().split()

        while r and a:
            if r[0] == a[0]:
                kmer, rc, ac = r[0], int(r[1]), int(a[1])
                total += 1

                # Reservoir sampling
                if len(sample) < sample_size:
                    sample.append((kmer, rc, ac))
                else:
                    j = random.randint(0, total - 1)
                    if j < sample_size:
                        sample[j] = (kmer, rc, ac)

                r = fr.readline().split()
                a = fa.readline().split()

            elif r[0] < a[0]:
                r = fr.readline().split()
            else:
                a = fa.readline().split()

    return sample, total


def stream_ranked(ranked_file: str, sample_size: int) -> Tuple[List[Tuple[str, int, int]], int]:
    """
    Stream through a ranked TSV file and keep a reservoir sample.
    """
    sample = []
    total = 0

    with open(ranked_file) as f:
        header = f.readline()  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                kmer, rc, ac = parts[0], int(parts[1]), int(parts[2])
                total += 1

                # Reservoir sampling
                if len(sample) < sample_size:
                    sample.append((kmer, rc, ac))
                else:
                    j = random.randint(0, total - 1)
                    if j < sample_size:
                        sample[j] = (kmer, rc, ac)

    return sample, total

if args.input_ranked:
    logger.info(f"Reading ranked table from {args.input_ranked} and sampling...")
    sample, total_pairs = stream_ranked(args.input_ranked, args.plot_sample)
elif args.raw_dump and args.asm_dump:
    logger.info("Streaming merge and sampling from jellyfish dumps...")
    sample, total_pairs = stream_merge(args.raw_dump, args.asm_dump, args.plot_sample)
else:
    logger.error("Error: Either -i (ranked table) or both -r and -a (jellyfish dumps) must be provided.")
    parser.print_help()
    exit(1)

logger.info(f"Processed {total_pairs:,} matched kmers; kept {len(sample):,} in sample.")


# Transformations on sample
df = pd.DataFrame(sample, columns=["kmer", "RawCount", "AsmCount"])
df["logRaw"] = np.log10(df["RawCount"] + 1)
df["logAsm"] = np.log10(df["AsmCount"] + 1)
df["reduction"] = df["logAsm"] - df["logRaw"]
df["reductionRank"] = df["reduction"].rank(method="first")

# Percentile filtering on sample
p = args.percentile / 100.0
n = len(df)
low_cut = int(math.ceil(n * p))
high_cut = int(math.floor(n * (1 - p)))
df_sorted = df.sort_values("reductionRank")
df_extreme = pd.concat([df_sorted.iloc[:low_cut], df_sorted.iloc[high_cut:]])
##
xmin = min(df_sorted["RawCount"])
xmax = max(df_sorted["RawCount"])
ymin = min(df_sorted["AsmCount"])
ymax = max(df_sorted["AsmCount"])
##
# Plotting functions
def plot_all(df: pd.DataFrame, df_extreme: pd.DataFrame, args: argparse.Namespace, xmin: float, xmax: float, ymin: float, ymax: float):
    sns.set(style="whitegrid")

    # Scatter
    plt.figure(figsize=(8, 8))
    plt.scatter(df["RawCount"], df["AsmCount"], s=1, alpha=0.2)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("Kmers in raw data")
    plt.ylabel("Kmers in assembly")
    plt.title(f"K={args.kmer_size} coverage (scatter)")
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_scatter.{args.output_format}", dpi=200)
    plt.close()

    # Scatter (extremes only)
    if not df_extreme.empty:
        plt.figure(figsize=(8, 8))
        plt.scatter(df_extreme["RawCount"], df_extreme["AsmCount"],
                    c=df_extreme["reduction"], cmap="coolwarm", s=2, alpha=0.6)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.xlabel("Kmers in raw data")
        plt.ylabel("Kmers in assembly")
        plt.title(f"K={args.kmer_size} extreme kmers (±{args.percentile}%)")
        plt.colorbar(label="Log-fold change")
        plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_scatter_extreme_{args.percentile}pct.{args.output_format}", dpi=200)
        plt.close()

    # ECDF
    plt.figure(figsize=(8, 6))
    x = np.sort(df["reduction"])
    y = np.arange(1, len(x) + 1) / len(x)
    plt.step(x, y, where="post", color="red")
    plt.xlabel("Log-fold kmer change")
    plt.ylabel("Cumulative probability")
    plt.title(f"K={args.kmer_size} empirical cumulative distribution")
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_ecdf.{args.output_format}", dpi=200)
    plt.close()

    # Density
    plt.figure(figsize=(8, 6))
    sns.kdeplot(df["reduction"], fill=True, alpha=0.6, color="#399602")
    plt.xlabel("Log-fold kmer change")
    plt.ylabel("Density")
    plt.title(f"K={args.kmer_size} coverage density")
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_density.{args.output_format}", dpi=200)
    plt.close()

    # Violin
    plt.figure(figsize=(6, 6))
    sns.violinplot(data=[df["logRaw"], df["logAsm"]], palette=["#1f77b4", "#ff7f0e"])
    plt.xticks([0, 1], ["Raw", "Asm"])
    plt.ylabel("Abundance (log10)")
    plt.title(f"K={args.kmer_size} coverage abundance")
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_violin.{args.output_format}", dpi=200)
    plt.close()

    # Back-to-back density plot
    plt.figure(figsize=(6, 6))

    # Compute KDE curves manually once
    y_vals = np.linspace(
        min(df["logRaw"].min(), df["logAsm"].min()),
        max(df["logRaw"].max(), df["logAsm"].max()),
        500
    )

    raw_kde = gaussian_kde(df["logRaw"], bw_method=0.5)(y_vals)
    asm_kde = gaussian_kde(df["logAsm"], bw_method=0.5)(y_vals)

    # Mirror plot: Raw to the left (negative x), Asm to the right
    plt.fill_betweenx(y_vals, -raw_kde, 0, color="#1f77b4", alpha=0.6, label="Raw")
    plt.fill_betweenx(y_vals, 0, asm_kde, color="#ff7f0e", alpha=0.6, label="Asm")

    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Density (mirrored)")
    plt.ylabel("Abundance (log10)")
    plt.title(f"K={args.kmer_size} back-to-back density")
    plt.legend()
    plt.savefig(f"{args.output_prefix}_k{args.kmer_size}_back2back_density.{args.output_format}", dpi=200)
    plt.close()

plot_all(df, df_extreme, args, xmin, xmax, ymin, ymax)

logger.info("Done.")