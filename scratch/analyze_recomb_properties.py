import os
import pandas as pd
import numpy as np

def analyze_recomb(tsv_path, virus_name):
    if not os.path.exists(tsv_path):
        print(f"File not found: {tsv_path}")
        return
        
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"\n================================================================================")
    print(f"ANALYZING {virus_name.upper()} VALIDATED RECOMBINANTS ({len(df)} rows)")
    print(f"================================================================================")
    
    if len(df) == 0:
        print("No rows to analyze.")
        return
        
    # Check column names
    print("Columns:", list(df.columns))
    
    # 1. Inspect P-values
    p_vals = pd.to_numeric(df["p_value"], errors="coerce")
    print("\nP-value Distribution:")
    print(f"  Min:  {p_vals.min():.2E}")
    print(f"  Max:  {p_vals.max():.2E}")
    print(f"  Mean: {p_vals.mean():.2E}")
    print(f"  Median: {p_vals.median():.2E}")
    print(f"  Number of p-values > 0.05: {(p_vals > 0.05).sum()} / {len(df)} ({(p_vals > 0.05).sum()/len(df)*100:.1f}%)")
    
    # 2. Inspect Parent Identities and Switching Deltas
    p1_rec = pd.to_numeric(df["p1_identity_recomb"], errors="coerce")
    p2_rec = pd.to_numeric(df["p2_identity_recomb"], errors="coerce")
    p1_non = pd.to_numeric(df["p1_identity_nonrecomb"], errors="coerce")
    p2_non = pd.to_numeric(df["p2_identity_nonrecomb"], errors="coerce")
    
    # Delta identities
    rec_delta = p2_rec - p1_rec
    nonrec_delta = p1_non - p2_non
    
    print("\nIdentity Differences (Support for Parent Switching):")
    print("  Recombinant Region (Parent 2 - Parent 1 Identity):")
    print(f"    Min:  {rec_delta.min():.4f}")
    print(f"    Max:  {rec_delta.max():.4f}")
    print(f"    Mean: {rec_delta.mean():.4f}")
    print(f"    Median: {rec_delta.median():.4f}")
    print(f"    Delta < 0.01 (trivially small support): {(rec_delta < 0.01).sum()} / {len(df)} ({(rec_delta < 0.01).sum()/len(df)*100:.1f}%)")
    print(f"    Delta < 0.005: {(rec_delta < 0.005).sum()} / {len(df)} ({(rec_delta < 0.005).sum()/len(df)*100:.1f}%)")
    
    print("  Non-Recombinant Region (Parent 1 - Parent 2 Identity):")
    print(f"    Min:  {nonrec_delta.min():.4f}")
    print(f"    Max:  {nonrec_delta.max():.4f}")
    print(f"    Mean: {nonrec_delta.mean():.4f}")
    print(f"    Median: {nonrec_delta.median():.4f}")
    print(f"    Delta < 0.01: {(nonrec_delta < 0.01).sum()} / {len(df)} ({(nonrec_delta < 0.01).sum()/len(df)*100:.1f}%)")
    
    # 3. Analyze Parent Genotypes
    # We can check if parents start with "ref_" or are query sequences
    p1_is_ref = df["parent_1"].astype(str).str.startswith("ref_")
    p2_is_ref = df["parent_2"].astype(str).str.startswith("ref_")
    print("\nParent type classification:")
    print(f"  Parent 1 is reference: {p1_is_ref.sum()} / {len(df)} ({p1_is_ref.sum()/len(df)*100:.1f}%)")
    print(f"  Parent 2 is reference: {p2_is_ref.sum()} / {len(df)} ({p2_is_ref.sum()/len(df)*100:.1f}%)")
    
    # Let's inspect a few sample rows
    print("\nSample recombinants (top 5 rows):")
    cols_to_show = ["sequence_id", "genotype", "parent_1", "parent_2", "p1_identity_recomb", "p2_identity_recomb", "p1_identity_nonrecomb", "p2_identity_nonrecomb", "p_value"]
    print(df[cols_to_show].head(5).to_string(index=False))

if __name__ == "__main__":
    analyze_recomb("/app/results/hbv/validated_recombinants.tsv", "hbv")
    analyze_recomb("/app/results/hcv/validated_recombinants.tsv", "hcv")
    analyze_recomb("/app/results/hev/validated_recombinants.tsv", "hev")
