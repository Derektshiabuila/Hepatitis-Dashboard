import os
import pandas as pd
import re

def clean_genotype(g, virus):
    if pd.isna(g):
        return None
    g = str(g).strip().upper()
    if not g:
        return None
    if virus == "hcv":
        match = re.match(r"([1-7])([A-Z]*)", g, re.IGNORECASE)
        if match:
            major = match.group(1)
            subtype = match.group(0).lower()
            return {"major": major, "subtype": subtype}
        return {"major": g, "subtype": g.lower()}
    return {"major": g, "subtype": g}

def inspect_recombs(virus, thresholds=[0.80, 0.85, 0.88]):
    recomb_file = f"/app/results/{virus}/validated_recombinants.tsv"
    anno_file = f"/app/refs/{virus}/ref_annotations.tsv"
    
    if not os.path.exists(recomb_file):
        print(f"File not found: {recomb_file}")
        return
        
    df = pd.read_csv(recomb_file, sep="\t")
    
    # Load annotations to resolve parent genotypes
    anno_map = {}
    if os.path.exists(anno_file):
        anno_df = pd.read_csv(anno_file, sep="\t")
        anno_df.columns = [c.lower() for c in anno_df.columns]
        id_col = "sequenceid" if "sequenceid" in anno_df.columns else anno_df.columns[0]
        geno_col = "genotype" if "genotype" in anno_df.columns else anno_df.columns[1]
        for _, row in anno_df.iterrows():
            seq_id = str(row[id_col]).strip()
            anno_map[seq_id] = row[geno_col]
            anno_map[seq_id.replace("ref_", "")] = row[geno_col]

    p1_rec = pd.to_numeric(df["p1_identity_recomb"], errors="coerce")
    p2_rec = pd.to_numeric(df["p2_identity_recomb"], errors="coerce")
    p1_non = pd.to_numeric(df["p1_identity_nonrecomb"], errors="coerce")
    p2_non = pd.to_numeric(df["p2_identity_nonrecomb"], errors="coerce")

    for t in thresholds:
        passed_idx = (p2_rec >= t) & (p1_non >= t)
        passed_df = df[passed_idx].copy()
        
        print(f"\n================================================================================")
        print(f"{virus.upper()} RECOMBINANTS WITH IDENTITY >= {t:.2f} (Count: {len(passed_df)})")
        print(f"================================================================================")
        
        if len(passed_df) == 0:
            print("None found.")
            continue
            
        # Classify and map parent genotypes
        classified = []
        for _, row in passed_df.iterrows():
            p1 = row["parent_1"]
            p2 = row["parent_2"]
            p1_geno = anno_map.get(p1, "Unknown")
            p2_geno = anno_map.get(p2, "Unknown")
            
            p1_info = clean_genotype(p1_geno, virus)
            p2_info = clean_genotype(p2_geno, virus)
            
            if p1_info and p2_info and p1_info["major"] != "UNKNOWN" and p2_info["major"] != "UNKNOWN":
                if p1_info["major"] != p2_info["major"]:
                    cls = "Inter-Genotypic"
                elif p1_info["subtype"] != p2_info["subtype"]:
                    cls = "Intra-Genotypic (Inter-Subtype)"
                else:
                    cls = "Intra-Subtype"
            else:
                cls = "Unknown"
                
            classified.append({
                "sequence_id": row["sequence_id"],
                "genotype": row["genotype"],
                "parent_1": p1,
                "p1_geno": p1_geno,
                "parent_2": p2,
                "p2_geno": p2_geno,
                "p2_id_recomb": row["p2_identity_recomb"],
                "p1_id_nonrecomb": row["p1_identity_nonrecomb"],
                "type": cls
            })
            
        class_df = pd.DataFrame(classified)
        print("Class counts:")
        print(class_df["type"].value_counts().to_string())
        
        print("\nSample sequences passing threshold:")
        print(class_df.head(10).to_string(index=False))

if __name__ == "__main__":
    inspect_recombs("hbv", [0.85])
    inspect_recombs("hcv", [0.80, 0.85, 0.88])
    inspect_recombs("hev", [0.80, 0.85])
