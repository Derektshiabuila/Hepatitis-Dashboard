import os
import sys
import pickle
import pandas as pd

# Ensure the root project directory is on python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data_loader import get_data_path, normalize_accession_id

print("================================================================================")
print("DASHBOARD RECOMBINATION DIAGNOSTIC CHECK")
print("================================================================================")

viruses = ["hbv", "hcv", "hev"]
for virus in viruses:
    recomb_file = f"results/{virus}/validated_recombinants.tsv"
    recomb_path = get_data_path(recomb_file)
    print(f"\n[{virus.upper()}] Recombination File Check:")
    print(f"  Path: {recomb_path}")
    if os.path.exists(recomb_path):
        size = os.path.getsize(recomb_path)
        print(f"  Status: EXISTS (Size: {size:,} bytes)")
        try:
            df = pd.read_csv(recomb_path, sep="\t", dtype=str)
            print(f"  Rows: {len(df):,}")
            print(f"  Columns: {list(df.columns)}")
            if len(df) > 0:
                id_col = None
                for col in ["sequence_id", "accession", "ID", "sample", "Sample"]:
                    if col in df.columns:
                        id_col = col
                        break
                print(f"  Matched ID column: {id_col}")
                print(f"  Sample IDs: {list(df[id_col].head(3)) if id_col else 'None'}")
                print(f"  Sample is_recombinant values: {list(df['is_recombinant'].head(3)) if 'is_recombinant' in df.columns else 'None'}")
        except Exception as e:
            print(f"  Error reading file: {e}")
    else:
        print("  Status: NOT FOUND")

cache_file = get_data_path("results/preprocessed_data_store.pkl")
print(f"\n[CACHE] Cache File Check:")
print(f"  Path: {cache_file}")
if os.path.exists(cache_file):
    size = os.path.getsize(cache_file)
    print(f"  Status: EXISTS (Size: {size:,} bytes)")
    try:
        with open(cache_file, "rb") as f:
            data = pickle.load(f)
        print("  Status: Loaded successfully")
        for key in ["hbv_data", "hcv_data", "hev_data"]:
            if key in data:
                df = data[key]
                print(f"  DataFrame '{key}':")
                print(f"    Total Rows: {len(df):,}")
                if "is_recombinant" in df.columns:
                    recomb_vals = df["is_recombinant"].value_counts().to_dict()
                    print(f"    is_recombinant value counts: {recomb_vals}")
                else:
                    print("    is_recombinant column: NOT FOUND")
                if "genotype" in df.columns:
                    recomb_genotypes = (df["genotype"] == "Recombinant").sum()
                    print(f"    Sequences with genotype 'Recombinant': {recomb_genotypes}")
                    genotype_vals = df["genotype"].value_counts().head(5).to_dict()
                    print(f"    Top 5 genotypes: {genotype_vals}")
                else:
                    print("    genotype column: NOT FOUND")
            else:
                print(f"  DataFrame '{key}': NOT FOUND in cache")
    except Exception as e:
        print(f"  Error loading cache: {e}")
else:
    print("  Status: NOT FOUND")
print("================================================================================")
