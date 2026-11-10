import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pyBigWig
import warnings

warnings.filterwarnings('ignore')

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'sans-serif']
mpl.rcParams['axes.linewidth'] = 1.5

def main():
    parser = argparse.ArgumentParser(description="Plot biological protein signal enrichment at TAD-X boundaries.")
    parser.add_argument("--bed", type=str, required=True, help="Path to the predicted TADs BED file")
    parser.add_argument("--bw", type=str, required=True, help="Path to the ChIP-seq bigWig file (e.g., CTCF)")
    parser.add_argument("--protein", type=str, default="CTCF", help="Target protein name for the plot title")
    parser.add_argument("--out_img", type=str, default="TAD_X_CTCF_Enrichment.pdf", help="Output image file (.pdf or .png)")
    parser.add_argument("--window", type=int, default=450000, help="Window size (bp) upstream/downstream of boundary")
    parser.add_argument("--bin_size", type=int, default=10000, help="Bin resolution (bp) for extracting signal")

    args = parser.parse_args()
    print(f"\n[*] Extracting {args.protein} signal around TAD-X boundaries...")

    if not os.path.exists(args.bed):
        print(f"[ERROR] BED file not found: {args.bed}")
        return
    if not os.path.exists(args.bw):
        print(f"[ERROR] bigWig file not found: {args.bw}")
        return

    try:
        df = pd.read_csv(args.bed, sep='\t', header=None, usecols=[0, 1, 2], names=['chrom', 'start', 'end'])
        df_melt = pd.concat([
            df[['chrom', 'start']].rename(columns={'start': 'pos'}),
            df[['chrom', 'end']].rename(columns={'end': 'pos'})
        ])
        boundaries = df_melt.drop_duplicates().values.tolist()
    except Exception as e:
        print(f"[ERROR] Failed to parse BED file: {e}")
        return

    n_bins = int((args.window * 2) / args.bin_size)
    x_axis = np.linspace(-args.window / 1000000, args.window / 1000000, n_bins)
    signal_matrix = []

    try:
        bw = pyBigWig.open(args.bw)
        for chrom, pos in boundaries:
            start = int(pos) - args.window
            end = int(pos) + args.window
            if start < 0:
                continue
            try:
                vals = bw.stats(chrom, start, end, type="mean", nBins=n_bins)
                vals = [0.0 if v is None else v for v in vals]
                signal_matrix.append(vals)
            except RuntimeError:
                pass
        bw.close()
    except Exception as e:
        print(f"[ERROR] Failed to process bigWig: {e}")
        return

    if not signal_matrix:
        print("[ERROR] No valid signal extracted. Check chromosome naming compatibility (e.g., 'chr7' vs '7').")
        return

    mean_signal = np.mean(signal_matrix, axis=0)
    print(f"[*] Successfully aggregated signal across {len(signal_matrix)} valid boundaries.")

    print(f"[*] Generating high-resolution {args.protein} enrichment plot...")
    fig, ax = plt.subplots(figsize=(5.5, 4.5))

    color = '#E63946' 

    ax.plot(x_axis, mean_signal, color=color, linewidth=2.5)
    ax.fill_between(x_axis, mean_signal, np.min(mean_signal), color=color, alpha=0.15)

    ax.axvline(x=0, color='black', linestyle='--', linewidth=1.2, alpha=0.6)
    ax.axvspan(-0.05, 0.05, color='gray', alpha=0.08, zorder=0)

    ax.set_title(f"TAD-X: {args.protein} Signal Enrichment", fontsize=16, fontweight='bold', pad=15)
    ax.set_xlabel("Distance from TAD boundary (Mb)", fontsize=14, fontweight='bold')
    ax.set_ylabel(f"Mean {args.protein} Signal", fontsize=14, fontweight='bold')

    ax.set_xlim([-args.window / 1000000, args.window / 1000000])
    ax.set_xticks([-args.window / 1000000, 0, args.window / 1000000])
    ax.tick_params(axis='both', labelsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(args.out_img, dpi=400, bbox_inches='tight', facecolor='white')
    print(f"✅ Plot saved to {args.out_img}\n")

if __name__ == '__main__':
    main()
