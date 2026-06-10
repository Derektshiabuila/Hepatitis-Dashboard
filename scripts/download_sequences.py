# download_sequences.py

import os
import pandas as pd
import re
from Bio import Entrez
from hepatitis_data_collection import (
    VirusFilter,
    fetch_all_virus_metadata,
    parse_response_to_dictionaries,
    fetch_sequences_batch
)

# -----------------------------
# Read Snakemake inputs
# -----------------------------
virus = snakemake.params.virus.lower()
email = snakemake.params.email
minlength = snakemake.params.minlength
maxlength = snakemake.params.maxlength

fasta_out = snakemake.output.fasta
meta_out  = snakemake.output.meta

Entrez.email = email

# -----------------------------
# Virus-specific configuration
# -----------------------------
if virus == "hbv":
    taxid = "10407"
    minlength = 2900
    maxlength = 3500  # HBV is typically ~3,200 bp
    exclude_terms = [
        "bat", "chiroptera", "rodent", "primate", "avian", "bird",
        "chimpanzee", "gorilla", "orangutan", "monkey", "macaca",
        "callithrix", "marmoset", "baboon", "gibbon", "lemur",
        "tupaia", "tree shrew", "squirrel", "woodchuck", "duck",
        "heron", "crane", "goose", "parrot", "finch"
    ]
elif virus == "hcv":
    taxid = "3052230"
    minlength = 8000
    maxlength = 11000  # HCV is typically ~9,600 bp
    exclude_terms = []
elif virus == "hev":
    taxid = "291484"
    minlength = 7000
    maxlength = 7500  # HeV is typically ~7,200 bp
    exclude_terms = []
else:
    raise ValueError(f"Unknown virus: {virus}")

print(f"Fetching NCBI metadata for {virus.upper()} (taxid={taxid})")

# -----------------------------
# Fetch metadata
# -----------------------------
filters = VirusFilter(taxon=taxid)
records = fetch_all_virus_metadata(filters=filters)

df = pd.DataFrame(parse_response_to_dictionaries(records))

# -----------------------------
# Filtering
# -----------------------------
df = df.dropna(subset=["location", "isolate_collection_date"])
df["length"] = df["length"].astype(float)
df = df[(df["length"] >= minlength) & (df["length"] <= maxlength)]

print(f"Filtered to {len(df)} sequences with length between {minlength} and {maxlength} bp")
print(f"Length distribution:\n{df['length'].describe()}")

print(f"Initial count after length filter: {len(df)}")

# -----------------------------
# DEBUG: Check DataFrame structure
# -----------------------------
print("\n=== DEBUG INFORMATION ===")
print(f"DataFrame shape: {df.shape}")
print(f"Columns available: {list(df.columns)}")

# Check what columns actually contain data
print("\n=== CHECKING COLUMN CONTENT ===")
for col in df.columns:
    if df[col].dtype == 'object':  # Text columns
        non_empty = df[col].notna().sum()
        sample_vals = df[col].dropna().unique()[:3]
        print(f"{col}: {non_empty} non-empty values. Samples: {sample_vals}")

# -----------------------------
# Filter for human sequences only (HBV-specific)
# -----------------------------
if virus == "hbv":
    print(f"\n=== APPLYING HUMAN FILTER FOR HBV ===")
    initial_count = len(df)
    
    # Reset index
    df = df.reset_index(drop=True)
    
    # Convert all text columns to lowercase
    text_columns = [col for col in df.columns if df[col].dtype == 'object']
    for col in text_columns:
        df[col] = df[col].astype(str).str.lower()
    
    # SIMPLE AND EFFECTIVE FILTERING
    # Since we don't have reliable "human" indicators, we'll just exclude non-human sequences
    
    # Track which sequences we're removing
    removed_counts = {}
    
    for term in exclude_terms:
        mask = pd.Series([False] * len(df))
        for col in text_columns:
            # Simple string contains without regex groups
            col_mask = df[col].str.contains(term, na=False)
            mask = mask | col_mask
        
        if mask.any():
            removed_counts[term] = mask.sum()
            df = df[~mask].copy()
    
    # Print summary
    print(f"\n=== FILTERING SUMMARY ===")
    print(f"Initial sequences: {initial_count}")
    print(f"Final sequences after filtering: {len(df)}")
    print(f"Total filtered out: {initial_count - len(df)}")
    
    if removed_counts:
        print("\nSequences removed by term:")
        for term, count in removed_counts.items():
            if count > 0:
                print(f"  '{term}': {count} sequences")
    
    # Quick check for any remaining bat sequences
    if len(df) > 0:
        bat_check = pd.Series([False] * len(df))
        for col in text_columns:
            bat_check = bat_check | df[col].str.contains("bat", na=False)
        
        if bat_check.any():
            print(f"\nWARNING: Still found {bat_check.sum()} sequences with 'bat' after filtering!")
            # Show what columns contain 'bat'
            for col in text_columns:
                if df[col].str.contains("bat", na=False).any():
                    bat_examples = df[df[col].str.contains("bat", na=False)]["accession_id"].head(3).tolist()
                    print(f"  Found in column '{col}': examples {bat_examples}")

# -----------------------------
# Write outputs SAFELY
# -----------------------------
tmp_fasta = fasta_out + ".tmp"
tmp_meta  = meta_out + ".tmp"

# Create output directory
os.makedirs(os.path.dirname(fasta_out), exist_ok=True)

# Metadata first
df.to_csv(tmp_meta, sep="\t", index=False)
print(f"\nMetadata saved to temporary file: {tmp_meta}")

# Sequence download
if len(df) > 0:
    accessions = df["accession_id"].tolist()
    print(f"Fetching {len(accessions)} sequences from NCBI...")
    
    try:
        fetch_sequences_batch(
            accession_list=accessions,
            batch_size=200,
            output_file=tmp_fasta,
            sleep_time=0.1
        )
        print("Sequences downloaded successfully.")
    except Exception as e:
        print(f"Error downloading sequences: {e}")
        # Clean up and exit with error
        if os.path.exists(tmp_fasta):
            os.remove(tmp_fasta)
        if os.path.exists(tmp_meta):
            os.remove(tmp_meta)
        raise
else:
    print("WARNING: No sequences to download after filtering!")
    # Create empty files
    with open(tmp_fasta, 'w') as f:
        pass
    with open(tmp_meta, 'w') as f:
        f.write("\t".join(df.columns) + "\n" if len(df.columns) > 0 else "")

# -----------------------------
# Atomic move to final outputs
# -----------------------------
os.replace(tmp_fasta, fasta_out)
os.replace(tmp_meta, meta_out)

print(f"\n=== DOWNLOAD COMPLETE ===")
print(f"FASTA file: {fasta_out}")
print(f"Metadata file: {meta_out}")
print(f"Total sequences: {len(df)}")

# Quick verification
if os.path.exists(fasta_out):
    # Count sequences in FASTA
    with open(fasta_out, 'r') as f:
        seq_count = sum(1 for line in f if line.startswith('>'))
    print(f"Sequences in FASTA file: {seq_count}")