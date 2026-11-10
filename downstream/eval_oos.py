import argparse
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.model_selection import KFold
import warnings

warnings.filterwarnings('ignore')

def run_5fold_cv_evaluation(df_eval):

    kf = KFold(n_splits=5, shuffle=True, random_state=31415)
    delta_lls = []
    eps = 1e-10

    for train_idx, test_idx in kf.split(df_eval):
        df_train = df_eval.iloc[train_idx]
        df_test = df_eval.iloc[test_idx].copy()

        # M0 (Distance-only baseline)
        lam0_train = df_train.groupby('dist')['count'].mean().rename('lam0')
        df_test = df_test.merge(lam0_train, on='dist', how='left')
        df_test['lam0'] = df_test['lam0'].fillna(df_train['count'].mean())
        LL0_test = np.sum(df_test['count'] * np.log(df_test['lam0'] + eps) - df_test['lam0'])

        # M1 (Distance + Domain Architecture)
        lam1_train = df_train.groupby(['dist', 'A'])['count'].mean().rename('lam1')
        df_test = df_test.merge(lam1_train, on=['dist', 'A'], how='left')
        
        df_test['lam1'] = df_test['lam1'].fillna(df_test['lam0'])
        LL1_test = np.sum(df_test['count'] * np.log(df_test['lam1'] + eps) - df_test['lam1'])

        delta_lls.append(LL1_test - LL0_test)

    return np.mean(delta_lls), np.std(delta_lls)

def main():
    parser = argparse.ArgumentParser(description="Evaluate out-of-sample information gain for predicted TADs.")
    parser.add_argument("--matrix", type=str, required=True, help="Path to the input Hi-C matrix file (.npz)")
    parser.add_argument("--bed", type=str, required=True, help="Path to the predicted TADs BED file")
    parser.add_argument("--res", type=int, required=True, help="Resolution in bp (e.g., 50000)")
    parser.add_argument("--max_dist", type=int, default=40, help="Maximum distance in bins to evaluate (default: 40)")
    
    args = parser.parse_args()

    print("\n[*] Loading data for Out-of-Sample (OOS) Evaluation...")
    

    try:
        mat_sparse = sp.load_npz(args.matrix)
        matrix = mat_sparse.toarray().astype(np.float64)
    except Exception as e:
        print(f"[ERROR] Failed to load matrix: {e}")
        return

    matrix = np.nan_to_num(matrix, nan=0.0)
    n_bins = matrix.shape[0]

    tads = []
    try:
        df_bed = pd.read_csv(args.bed, sep='\t', header=None, usecols=[0,1,2], names=['chrom', 'start', 'end'])
        for _, row in df_bed.iterrows():
            s = int(row['start']) // args.res
            e = int(row['end']) // args.res - 1
            if s < e:
                tads.append((s, e))
    except Exception as e:
        print(f"[ERROR] Failed to parse BED file: {e}")
        return

    A_mat = np.zeros((n_bins, n_bins), dtype=np.int32)
    for s, e in tads:
        A_mat[s:e+1, s:e+1] += 1


    rows, cols = np.triu_indices(n_bins, k=0)
    dists = cols - rows
    mask = dists <= args.max_dist

    df_base = pd.DataFrame({
        'row': rows[mask],
        'col': cols[mask],
        'dist': dists[mask],
        'count': matrix[rows[mask], cols[mask]]
    })
    df_base['A'] = A_mat[df_base['row'], df_base['col']]

    print(f"[*] Extracted {len(df_base):,} valid pixel pairs (distance <= {args.max_dist} bins).")
    print("[*] Running 5-fold cross-validation. This may take a moment...")
    
    mean_dll, std_dll = run_5fold_cv_evaluation(df_base)

    print("-" * 50)
    print(f"✅ Evaluation Results:")
    print(f"   - Model Complexity (Number of TADs): {len(tads)}")
    print(f"   - OOS Information Gain (\u0394 LL):  {mean_dll:,.2f} \u00B1 {std_dll:,.2f}")
    print("-" * 50)

if __name__ == '__main__':
    main()
