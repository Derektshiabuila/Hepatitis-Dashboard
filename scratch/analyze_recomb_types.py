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
        # HCV genotype annotation like "1a", "1b", "2a"
        match = re.match(r"([1-7])([A-Z]*)", g, re.IGNORECASE)
        if match:
            major = match.group(1)
            subtype = match.group(0).lower()
            return {"major": major, "subtype": subtype}
        return {"major": g, "subtype": g.lower()}
        
    elif virus == "hbv":
        # HBV genotype annotation like "A", "B", "C"
        # Just use genotype letter as major and subtype
        return {"major": g[0], "subtype": g}
        
    elif virus == "hev":
        # HEV genotype annotation like "4c", "3e", "3"
        match = re.match(r"([1-8])([A-Z]*)", g, re.IGNORECASE)
        if match:
            major = match.group(1)
            subtype = match.group(0).lower()
            return {"major": major, "subtype": subtype}
        return {"major": g, "subtype": g.lower()}
        
    return {"major": g, "subtype": g}

def analyze_recomb_types(virus):
    recomb_file = f"/app/results/{virus}/validated_recombinants.tsv"
    anno_file = f"/app/refs/{virus}/ref_annotations.tsv"
    
    if not os.path.exists(recomb_file):
        print(f"\n[{virus.upper()}] Recombination file not found: {recomb_file}")
        return
        
    if not os.path.exists(anno_file):
        # Fallback to local path just in case
        anno_file = f"refs/{virus}/ref_annotations.tsv"
        if not os.path.exists(anno_file):
            # Try to build mapping from ref_msa headers
            print(f"\n[{virus.upper()}] Annotations file not found: {anno_file}")
            return

    # Load annotations
    # Read columns: SequenceID (or sequenceID) and Genotype (or genotype)
    anno_df = pd.read_csv(anno_file, sep="\t")
    # standardize columns
    anno_df.columns = [c.lower() for c in anno_df.columns]
    id_col = "sequenceid" if "sequenceid" in anno_df.columns else anno_df.columns[0]
    geno_col = "genotype" if "genotype" in anno_df.columns else anno_df.columns[1]
    
    # Create map from reference accession to genotype
    anno_map = {}
    for _, row in anno_df.iterrows():
        seq_id = str(row[id_col]).strip()
        # Clean prefix ref_ if present in map lookup
        clean_id = seq_id.replace("ref_", "")
        anno_map[seq_id] = row[geno_col]
        anno_map[clean_id] = row[geno_col]

    # Load validated recombinants
    rec_df = pd.read_csv(recomb_file, sep="\t")
    
    print(f"\n================================================================================")
    print(f"CLASSIFYING {virus.upper()} RECOMBINATION EVENTS ({len(rec_df)} validated recombinants)")
    print(f"================================================================================")
    
    if len(rec_df) == 0:
        print("No recombinants found.")
        return
        
    counts = {
        "inter_genotypic": 0,
        "intra_genotypic_inter_subtype": 0,
        "intra_subtype": 0,
        "unknown": 0
    }
    
    examples = {
        "inter_genotypic": [],
        "intra_genotypic_inter_subtype": [],
        "intra_subtype": []
    }
    
    for _, row in rec_df.iterrows():
        cand_id = row["sequence_id"]
        p1 = row["parent_1"]
        p2 = row["parent_2"]
        cand_geno_raw = row["genotype"]
        
        # Get genotypes of parents
        p1_geno_raw = anno_map.get(p1, None)
        p2_geno_raw = anno_map.get(p2, None)
        
        if p1_geno_raw is None or p2_geno_raw is None:
            counts["unknown"] += 1
            continue
            
        p1_info = clean_genotype(p1_geno_raw, virus)
        p2_info = clean_genotype(p2_geno_raw, virus)
        
        if p1_info is None or p2_info is None:
            counts["unknown"] += 1
            continue
            
        classification = None
        if p1_info["major"] != p2_info["major"]:
            classification = "inter_genotypic"
        elif p1_info["subtype"] != p2_info["subtype"]:
            classification = "intra_genotypic_inter_subtype"
        else:
            classification = "intra_subtype"
            
        counts[classification] += 1
        
        if len(examples[classification]) < 5:
            examples[classification].append({
                "child": cand_id,
                "child_geno": cand_geno_raw,
                "p1": p1,
                "p1_geno": p1_geno_raw,
                "p2": p2,
                "p2_geno": p2_geno_raw,
                "p_value": row["p_value"]
            })
            
    print("\nRecombination Class Counts:")
    for cls, count in counts.items():
        pct = count / len(rec_df) * 100
        print(f"  {cls:30s}: {count:5d} ({pct:.1f}%)")
        
    for cls in ["inter_genotypic", "intra_genotypic_inter_subtype", "intra_subtype"]:
        if examples[cls]:
            print(f"\nExamples of {cls.upper()}:")
            print(pd.DataFrame(examples[cls]).to_string(index=False))

if __name__ == "__main__":
    analyze_recomb_types("hbv")
    analyze_recomb_types("hcv")
    analyze_recomb_types("hev")
