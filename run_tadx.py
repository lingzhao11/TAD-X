import argparse
import time
import numpy as np
import scipy.sparse as sp
import pandas as pd
from tad_x import TAD_X 

def main():
    parser = argparse.ArgumentParser(description="TAD-X: FFT-based identification of complex topologically associating domains")
    
    parser.add_argument("--matrix", dest="matrix", type=str, required=True, help="Path to the input Hi-C matrix file (.npz)")
    parser.add_argument("--chrom", dest="chrom", type=str, default="chr7", help="Chromosome name (e.g., chr7)")
    parser.add_argument("--res", dest="resolution", type=int, required=True, help="Resolution of the Hi-C matrix in bp (e.g., 50000)")
    parser.add_argument("--out_bed", dest="out_bed", type=str, default="TAD_X_domains.bed", help="Output BED file for predicted TADs")
    parser.add_argument("--prominence", dest="prominence", type=float, default=0.25, help="Prominence factor for boundary peak detection")
    
    args = parser.parse_args()

    print(f"\n{'='*60}")
    print(f"🚀 TAD-X: Multiscale Topological Domain Inference")
    print(f"{'='*60}")

    try:
        print(f"[*] Loading Hi-C matrix from: {args.matrix} ...")
        mat_sparse = sp.load_npz(args.matrix)
        mat = mat_sparse.toarray().astype(np.float32)
        print(f"    - Matrix shape: {mat.shape}")
    except Exception as e:
        print(f"\n[ERROR] Failed to load .npz matrix: {e}")
        return

    print(f"[*] Running FFT-based spatial scale discovery & DP optimization...")
    start_time = time.perf_counter()
    
    mod = TAD_X(resolution=args.resolution, prominence_factor=args.prominence)
    result = mod.detect(mat)
    
    pred_bin_tads = result["all_tads"]
    
    cost_time = time.perf_counter() - start_time
    print(f"[*] Inference completed in {cost_time:.2f} seconds.")

    pred_tads = []
    for s_bin, e_bin in pred_bin_tads:
        s_bp = int(s_bin * args.resolution)
        e_bp = int((e_bin + 1) * args.resolution)
        pred_tads.append({
            "chr": args.chrom,
            "start": s_bp,
            "end": e_bp
        })
        
    print("-" * 60)
    print(f"✅ Total multiscale domains identified: {len(pred_tads)}")
    print("-" * 60)

    if pred_tads:
        df_out = pd.DataFrame(pred_tads)
        df_out = df_out.sort_values(by=['start', 'end']).drop_duplicates()
        df_out.to_csv(args.out_bed, sep='\t', header=False, index=False)
        print(f"[*] Saved boundaries to {args.out_bed}\n")
    else:
        print("[!] No domains detected for this matrix.\n")

if __name__ == '__main__':
    main()
