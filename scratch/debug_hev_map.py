import os
import pandas as pd

def debug():
    recomb_file = "/app/results/hev/validated_recombinants.tsv"
    anno_file = "/app/refs/hev/ref_annotations.tsv"
    
    if not os.path.exists(recomb_file) or not os.path.exists(anno_file):
        print("Files not found.")
        return
        
    anno_df = pd.read_csv(anno_file, sep="\t")
    print("Annotations columns:", list(anno_df.columns))
    print("First 3 annotations:")
    print(anno_df.head(3).to_dict(orient="records"))
    
    rec_df = pd.read_csv(recomb_file, sep="\t")
    print("\nRecombinants columns:", list(rec_df.columns))
    print("First 3 recombinants:")
    print(rec_df[["sequence_id", "parent_1", "parent_2"]].head(3).to_dict(orient="records"))
    
    p1 = rec_df.loc[0, "parent_1"]
    print(f"\nParent 1: '{p1}' (type: {type(p1)})")
    
    # Check if parent_1 exists in anno_df sequenceID
    id_col = anno_df.columns[0]
    matched = anno_df[anno_df[id_col].astype(str).str.strip() == p1.replace("ref_", "")]
    print("Match by strip ref_:", matched.to_dict(orient="records"))

if __name__ == "__main__":
    debug()
