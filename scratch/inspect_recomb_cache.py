import pickle
import os

def inspect_cache():
    cache_path = "/app/results/preprocessed_data_store.pkl"
    if not os.path.exists(cache_path):
        # Fallback to local
        cache_path = "results/preprocessed_data_store.pkl"
        if not os.path.exists(cache_path):
            print("Cache file not found.")
            return
            
    print(f"Loading cache file: {cache_path}")
    with open(cache_path, "rb") as f:
        data = pickle.load(f)
        
    for key in ["hbv_data", "hcv_data", "hev_data"]:
        if key in data:
            df = data[key]
            print(f"\n[{key.upper()}] Shape: {df.shape}")
            print("Columns:", list(df.columns))
            
            # Check is_recombinant
            if "is_recombinant" in df.columns:
                print("is_recombinant value counts:")
                print(df["is_recombinant"].value_counts().to_dict())
            else:
                print("is_recombinant column: NOT FOUND")
                
            # Check genotype Recombinant
            if "genotype" in df.columns:
                recomb_count = (df["genotype"] == "Recombinant").sum()
                print(f"Number of rows with genotype 'Recombinant': {recomb_count}")
            
            # Check recombination_class
            if "recombination_class" in df.columns:
                print("recombination_class value counts:")
                print(df["recombination_class"].value_counts().to_dict())
                
                # Check for Recombinant rows specifically
                if "genotype" in df.columns:
                    recomb_rows = df[df["genotype"] == "Recombinant"]
                    if len(recomb_rows) > 0:
                        print("recombination_class value counts for Recombinant sequences:")
                        print(recomb_rows["recombination_class"].value_counts().to_dict())
            else:
                print("recombination_class column: NOT FOUND")
        else:
            print(f"{key} not found in cache.")

if __name__ == "__main__":
    inspect_cache()
