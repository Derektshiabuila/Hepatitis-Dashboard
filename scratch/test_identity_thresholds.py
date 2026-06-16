import os
import pandas as pd

def test_thresholds(virus):
    tsv_path = f"/app/results/{virus}/validated_recombinants.tsv"
    if not os.path.exists(tsv_path):
        print(f"File not found: {tsv_path}")
        return
        
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"\n================================================================================")
    print(f"TESTING IDENTITY THRESHOLDS FOR {virus.upper()} (Total: {len(df)})")
    print(f"================================================================================")
    
    p1_rec = pd.to_numeric(df["p1_identity_recomb"], errors="coerce")
    p2_rec = pd.to_numeric(df["p2_identity_recomb"], errors="coerce")
    p1_non = pd.to_numeric(df["p1_identity_nonrecomb"], errors="coerce")
    p2_non = pd.to_numeric(df["p2_identity_nonrecomb"], errors="coerce")
    
    thresholds = [0.70, 0.75, 0.80, 0.82, 0.85, 0.88, 0.90]
    
    for t in thresholds:
        # A sequence is valid if:
        # In the recomb region, it is close to Parent 2 (minor parent) -> identity >= t
        # In the non-recomb region, it is close to Parent 1 (major parent) -> identity >= t
        passed = (p2_rec >= t) & (p1_non >= t)
        count = passed.sum()
        pct = count / len(df) * 100
        print(f"  Threshold >= {t:.2f} (80%): {count:5d} ({pct:.1f}%) remaining")

if __name__ == "__main__":
    test_thresholds("hbv")
    test_thresholds("hcv")
    test_thresholds("hev")
