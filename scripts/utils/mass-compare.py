#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import os
import logging
from matplotlib.colors import LinearSegmentedColormap

# Logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger()

def calculate_jaccard(df, threshold):
    high_present = df['High_Density'] > threshold
    low_present = df['Low_Density'] > threshold
    intersection = np.logical_and(high_present, low_present).sum()
    union = np.logical_or(high_present, low_present).sum()
    return intersection / union if union > 0 else 0.0

def plot_scatter_ends(win_df, output_prefix, end_threshold):
    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(2, 2,  width_ratios=(1, 5), height_ratios=(5, 1),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_hist_x = fig.add_subplot(gs[1, 1], sharex=ax_scatter)
    ax_hist_y = fig.add_subplot(gs[0, 0], sharey=ax_scatter)

    # Background windows (beyond threshold)
    bg = win_df[win_df['MinDist'] > end_threshold]
    ax_scatter.scatter(bg['Low_Density'], bg['High_Density'], alpha=0.3, s=10, color='blue', label='Background')

    # End/N-adjacent regions (within threshold)
    ends = win_df[win_df['MinDist'] <= end_threshold]
    if not ends.empty:
        dists = np.clip(ends['MinDist'], 0, end_threshold)
        norm_dists = dists / end_threshold
        colors = [(.85, 0, 0), (1, 0.75, 0.8)] # Red to Pink
        cm = LinearSegmentedColormap.from_list('outlier_cm', colors, N=100)
        sc = ax_scatter.scatter(ends['Low_Density'], ends['High_Density'], c=norm_dists, cmap=cm, s=10, alpha=0.8, label='End/N-Adjacent', vmin=0, vmax=1)
        cbar_ax = fig.add_axes([0.92, 0.2, 0.02, 0.6])
        cbar = fig.colorbar(sc, cax=cbar_ax)
        cbar.set_label(f'Distance from Feature (0 to {end_threshold})')

    # Add diagonal line
    ax_scatter.plot([0, 1], [0, 1], 'k--', alpha=0.5)

    ax_scatter.set_xlim(-0.01, 1.01)
    ax_scatter.set_ylim(-0.01, 1.01)
    ax_scatter.tick_params(labelleft=False, labelbottom=False)

    # Marginal Histograms
    bins = np.linspace(0, 1, 50)
    ax_hist_x.hist(win_df['Low_Density'], bins=bins, color='blue', alpha=0.7)
    ax_hist_y.hist(win_df['High_Density'], bins=bins, orientation='horizontal', color='blue', alpha=0.7)
    ax_hist_y.invert_xaxis()

    ax_hist_x.set_xlabel('Low Extreme Kmer Density')
    ax_hist_y.set_ylabel('High Extreme Kmer Density')
    ax_scatter.set_title('Scatter Plot: Extreme Kmer Density (Ends Highlighted)', pad=20)

    plt.savefig(f"{output_prefix}_ends_scatter.png", dpi=300)
    plt.close()

def plot_scatter_outliers(win_df, output_prefix, end_threshold):
    fig = plt.figure(figsize=(12, 12))
    gs = fig.add_gridspec(2, 2,  width_ratios=(1, 5), height_ratios=(5, 1),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    ax_scatter = fig.add_subplot(gs[0, 1])
    ax_hist_x = fig.add_subplot(gs[1, 1], sharex=ax_scatter)
    ax_hist_y = fig.add_subplot(gs[0, 0], sharey=ax_scatter)

    # Background windows (non-outliers)
    bg = win_df[~win_df['IsOutlier']]
    ax_scatter.scatter(bg['Low_Density'], bg['High_Density'], alpha=0.3, s=10, color='blue', label='Background')

    # Outliers
    outliers = win_df[win_df['IsOutlier']]
    if not outliers.empty:
        # Near feature outliers (Red-to-Pink)
        near_outliers = outliers[outliers['MinDist'] <= end_threshold]
        if not near_outliers.empty:
            dists = np.clip(near_outliers['MinDist'], 0, end_threshold)
            norm_dists = dists / end_threshold
            colors = [(.85, 0, 0), (1, 0.75, 0.8)] # Red to Pink
            cm = LinearSegmentedColormap.from_list('outlier_cm', colors, N=100)
            sc = ax_scatter.scatter(near_outliers['Low_Density'], near_outliers['High_Density'], c=norm_dists, cmap=cm, s=15, alpha=0.9, label='Near-Feature Outlier', vmin=0, vmax=1)
            cbar_ax = fig.add_axes([0.92, 0.2, 0.02, 0.6])
            cbar = fig.colorbar(sc, cax=cbar_ax)
            cbar.set_label(f'Distance from Feature (0 to {end_threshold})')

        # Far outliers (Green)
        far_outliers = outliers[outliers['MinDist'] > end_threshold]
        if not far_outliers.empty:
            ax_scatter.scatter(far_outliers['Low_Density'], far_outliers['High_Density'], color='green', s=15, alpha=0.9, label='Far Outlier')

    # Add diagonal line
    ax_scatter.plot([0, 1], [0, 1], 'k--', alpha=0.5)

    ax_scatter.set_xlim(-0.01, 1.01)
    ax_scatter.set_ylim(-0.01, 1.01)
    ax_scatter.tick_params(labelleft=False, labelbottom=False)

    # Marginal Histograms
    bins = np.linspace(0, 1, 50)
    ax_hist_x.hist(win_df['Low_Density'], bins=bins, color='blue', alpha=0.7)
    ax_hist_y.hist(win_df['High_Density'], bins=bins, orientation='horizontal', color='blue', alpha=0.7)
    ax_hist_y.invert_xaxis()

    ax_hist_x.set_xlabel('Low Extreme Kmer Density')
    ax_hist_y.set_ylabel('High Extreme Kmer Density')
    ax_scatter.set_title('Scatter Plot: Extreme Kmer Density (Outliers Highlighted)', pad=20)

    plt.savefig(f"{output_prefix}_outliers_scatter.png", dpi=300)
    plt.close()

def plot_end_comparison(stats, output_prefix):
    labels = ['High (Background)', 'High (End)', 'Low (Background)', 'Low (End)']
    values = [
        stats['High_Density_Background'], stats['High_Density_End'],
        stats['Low_Density_Background'], stats['Low_Density_End']
    ]
    colors = ['#ffcccc', 'red', '#ccccff', 'blue']
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(labels, values, color=colors)
    ax.set_ylabel('Mean Density (bp / effective window bp)')
    ax.set_title('Comparison of Extreme Kmer Densities: Sequence Ends vs Background')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_end_comparison.png", dpi=300)
    plt.close()

def parse_ngaps(gff_path):
    gaps = {}
    if not gff_path or not os.path.exists(gff_path):
        return gaps
    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                seqid = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                if seqid not in gaps:
                    gaps[seqid] = []
                gaps[seqid].append((start, end))
        for seqid in gaps:
            gaps[seqid] = sorted(gaps[seqid])
    except Exception as e:
        logger.error(f"Error parsing N-gaps GFF: {e}")
    return gaps

def get_gap_overlap(w_start, w_end, seq_gaps):
    overlap = 0
    for g_start, g_end in seq_gaps:
        o_start = max(w_start, g_start)
        o_end = min(w_end, g_end)
        if o_start <= o_end:
            overlap += (o_end - o_start + 1)
    return overlap

def calculate_nearest_feature(w_start, w_end, length, seq_gaps):
    min_dist = w_start - 1
    f_type = "End-Adjacent"
    f_id = "Left-End"
    r_dist = length - w_end
    if r_dist < min_dist:
        min_dist = r_dist
        f_id = "Right-End"
    for idx, (g_start, g_end) in enumerate(seq_gaps):
        if g_end < w_start:
            d = w_start - g_end - 1
        elif g_start > w_end:
            d = g_start - w_end - 1
        else:
            d = 0
        if d < min_dist:
            min_dist = d
            f_type = "N-Adjacent"
            f_id = f"Gap-{idx}"
            if min_dist == 0:
                break
    return max(0, min_dist), f_type, f_id

def main():
    parser = argparse.ArgumentParser(description="Kmer Mass Compare: Compare extreme kmer accumulations across genome windows")
    parser.add_argument('-i', '--input', required=True, help='Input TSV from mass-query.py')
    parser.add_argument('-o', '--output', default='mass_compare', help='Output prefix for plots and stats')
    parser.add_argument('--jaccard-minimum', type=float, default=0, help='Minimum density threshold for Jaccard coincidence [default 0]')
    parser.add_argument('--ngaps', help='GFF file of N-gap coordinates')
    parser.add_argument('--end-threshold', type=int, default=0, help='Distance threshold for filtering outliers and coloring scatter [default 0]')
    parser.add_argument('--value-column', dest='value_column', default='Basepairs', choices=['Count', 'Basepairs'], help='Column to use for values [default Basepairs]')
    parser.add_argument('--outlier-method', choices=['SD', 'IQR'], default='SD', help='Method for joint outlier detection [default SD]')
    parser.add_argument('--outlier-stat', type=float, default=1.0, help='Multiplier for outlier detection [default 1.0]')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose mode')

    args = parser.parse_args()
    if args.verbose:
        logger.setLevel(logging.INFO)
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found.")
        return

    df = pd.read_csv(args.input, sep='\t')
    if args.value_column not in df.columns:
        logger.error(f"Value column '{args.value_column}' not found in input. Defaulting to 'Count'.")
        args.value_column = 'Count'

    pivot_df = df.pivot_table(index=['Sequence', 'Start', 'End'], columns='Bin', values=args.value_column, fill_value=0).reset_index()
    if 'high' not in pivot_df.columns: pivot_df['high'] = 0
    if 'low' not in pivot_df.columns: pivot_df['low'] = 0

    gaps = parse_ngaps(args.ngaps) if args.ngaps else {}
    seq_lengths = pivot_df.groupby('Sequence')['End'].max().to_dict()

    processed_windows = []
    for idx, row in pivot_df.iterrows():
        seq = row['Sequence']
        w_start = row['Start']
        w_end = row['End']
        seq_gaps = gaps.get(seq, [])
        # We now have EffectiveLen in the TSV from mass-query.py, but for backwards compatibility
        # or if it was modified, we can still use it.
        if 'EffectiveLen' in row:
            effective_len = row['EffectiveLen']
        else:
            gap_len = get_gap_overlap(w_start, w_end, seq_gaps)
            w_len = w_end - w_start + 1
            effective_len = w_len - gap_len

        if effective_len <= 0: continue

        min_dist, f_type, f_id = calculate_nearest_feature(w_start, w_end, seq_lengths[seq], seq_gaps)

        processed_windows.append({
            'Sequence': seq, 'Start': w_start, 'End': w_end,
            'high': row['high'], 'low': row['low'],
            'MinDist': min_dist, 'FeatureType': f_type, 'FeatureID': f"{seq}:{f_id}",
            'EffectiveLen': effective_len,
            'High_Density': row['high'] / effective_len,
            'Low_Density': row['low'] / effective_len
        })

    win_df = pd.DataFrame(processed_windows)
    if len(win_df) < 1:
        logger.error("No valid windows found.")
        return

    # Background Statistics
    if args.outlier_method == 'SD':
        bg_high_stat = win_df['High_Density'].std()
        bg_low_stat = win_df['Low_Density'].std()
        bg_high_center = win_df['High_Density'].mean()
        bg_low_center = win_df['Low_Density'].mean()
    else: # IQR
        bg_high_stat = win_df['High_Density'].quantile(0.75) - win_df['High_Density'].quantile(0.25)
        bg_low_stat = win_df['Low_Density'].quantile(0.75) - win_df['Low_Density'].quantile(0.25)
        bg_high_center = win_df['High_Density'].quantile(0.75)
        bg_low_center = win_df['Low_Density'].quantile(0.75)

    # Per-sequence background stats
    seq_backgrounds = {}
    for seq, group in win_df.groupby('Sequence'):
        if args.outlier_method == 'SD':
            s_high_stat = group['High_Density'].std()
            s_low_stat = group['Low_Density'].std()
            s_high_center = group['High_Density'].mean()
            s_low_center = group['Low_Density'].mean()
        else:
            s_high_stat = group['High_Density'].quantile(0.75) - group['High_Density'].quantile(0.25)
            s_low_stat = group['Low_Density'].quantile(0.75) - group['Low_Density'].quantile(0.25)
            s_high_center = group['High_Density'].quantile(0.75)
            s_low_center = group['Low_Density'].quantile(0.75)

        seq_backgrounds[seq] = {
            'high_stat': s_high_stat if s_high_stat > 0 else 1e-9,
            'low_stat': s_low_stat if s_low_stat > 0 else 1e-9,
            'high_center': s_high_center,
            'low_center': s_low_center
        }

    # Outlier flagging
    outlier_rows = []
    win_df['IsOutlier'] = False

    bg_high_stat = bg_high_stat if bg_high_stat > 0 else 1e-9
    bg_low_stat = bg_low_stat if bg_low_stat > 0 else 1e-9

    for idx, row in win_df.iterrows():
        seq = row['Sequence']
        s_bg = seq_backgrounds[seq]

        # Individual checks (for tagging)
        is_global_high = row['High_Density'] > (bg_high_center + args.outlier_stat * bg_high_stat)
        is_seq_high = row['High_Density'] > (s_bg['high_center'] + args.outlier_stat * s_bg['high_stat'])
        is_global_low = row['Low_Density'] > (bg_low_center + args.outlier_stat * bg_low_stat)
        is_seq_low = row['Low_Density'] > (s_bg['low_center'] + args.outlier_stat * s_bg['low_stat'])

        # Joint checks
        dist_global = np.sqrt(((max(0, row['High_Density'] - bg_high_center) / bg_high_stat)**2) +
                              ((max(0, row['Low_Density'] - bg_low_center) / bg_low_stat)**2))
        dist_seq = np.sqrt(((max(0, row['High_Density'] - s_bg['high_center']) / s_bg['high_stat'])**2) +
                           ((max(0, row['Low_Density'] - s_bg['low_center']) / s_bg['low_stat'])**2))

        if dist_global > args.outlier_stat or dist_seq > args.outlier_stat:
            win_df.at[idx, 'IsOutlier'] = True
            out_types = []
            if dist_global > args.outlier_stat: out_types.append("Global-Joint")
            if dist_seq > args.outlier_stat: out_types.append("Sequence-Joint")
            if is_global_high: out_types.append("Global-High")
            if is_seq_high: out_types.append("Sequence-High")
            if is_global_low: out_types.append("Global-Low")
            if is_seq_low: out_types.append("Sequence-Low")

            outlier_rows.append({
                'Sequence': seq, 'Length': seq_lengths[seq], 'Window': f"{row['Start']}-{row['End']}",
                'High_Density': row['High_Density'], 'Low_Density': row['Low_Density'],
                'MinDist': row['MinDist'], 'FeatureType': row['FeatureType'],
                'FeatureID': row['FeatureID'], 'OutlierType': ",".join(out_types),
                'JointDist_Global': dist_global, 'JointDist_Seq': dist_seq
            })

    df_outliers = pd.DataFrame(outlier_rows)

    # Sequence-level Summary Stats
    seq_stats_list = []
    for seq, group in win_df.groupby('Sequence'):
        s_high = group['high'].sum()
        s_low = group['low'].sum()
        s_asym = (s_high - s_low) / (s_high + s_low) if (s_high + s_low) > 0 else 0

        s_end = group[group['MinDist'] == 0]
        s_bg = group[group['MinDist'] > 0]

        s_high_end = s_end['High_Density'].mean() if not s_end.empty else 0
        s_high_bg = s_bg['High_Density'].mean() if not s_bg.empty else group['High_Density'].mean()
        s_low_end = s_end['Low_Density'].mean() if not s_end.empty else 0
        s_low_bg = s_bg['Low_Density'].mean() if not s_bg.empty else group['Low_Density'].mean()

        seq_stats_list.append({
            'Sequence': seq, 'Length': seq_lengths[seq], 'Total_High': s_high, 'Total_Low': s_low,
            'Asymmetry_Index': s_asym, 'High_Density_End': s_high_end, 'High_Density_Background': s_high_bg,
            'Low_Density_End': s_low_end, 'Low_Density_Background': s_low_bg
        })

    df_seq_stats = pd.DataFrame(seq_stats_list)

    total_high = win_df['high'].sum()
    total_low = win_df['low'].sum()

    # Output files
    stats_file = f"{args.output}.stats"
    with open(stats_file, 'w') as f:
        f.write("# Global Statistics\n")
        f.write(f"Global_Asymmetry_Index\t{(total_high - total_low) / (total_high + total_low) if (total_high + total_low) > 0 else 0:.4f}\n")
        f.write(f"Pearson_Correlation\t{pearsonr(win_df['Low_Density'], win_df['High_Density'])[0] if len(win_df)>1 else 0:.4f}\n")
        f.write(f"Jaccard_Coincidence_Index(>{args.jaccard_minimum})\t{calculate_jaccard(win_df, args.jaccard_minimum):.4f}\n")
        f.write(f"Global_Background_High_Center\t{bg_high_center:.4f}\n")
        f.write(f"Global_Background_High_Stat\t{bg_high_stat:.4f}\n")
        f.write(f"Global_Background_Low_Center\t{bg_low_center:.4f}\n")
        f.write(f"Global_Background_Low_Stat\t{bg_low_stat:.4f}\n")
        f.write("\n# Per-Sequence Statistics\n")
    df_seq_stats.to_csv(stats_file, sep='\t', index=False, mode='a')

    outliers_file = f"{args.output}.outliers.tsv"
    df_outliers.to_csv(outliers_file, sep='\t', index=False)

    logger.info(f"Statistics written to {stats_file}")
    logger.info(f"Outliers written to {outliers_file}")

    # Plotting
    plot_scatter_ends(win_df, args.output, args.end_threshold)
    plot_scatter_outliers(win_df, args.output, args.end_threshold)

    global_end_stats = {
        'High_Density_End': win_df[win_df['MinDist']==0]['High_Density'].mean(),
        'Low_Density_End': win_df[win_df['MinDist']==0]['Low_Density'].mean(),
        'High_Density_Background': win_df[win_df['MinDist']>0]['High_Density'].mean(),
        'Low_Density_Background': win_df[win_df['MinDist']>0]['Low_Density'].mean()
    }
    plot_end_comparison(global_end_stats, args.output)

    logger.info(f"Plots generated with prefix {args.output}")

if __name__ == "__main__":
    main()
