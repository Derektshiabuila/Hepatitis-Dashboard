import pandas as pd
import sys

csv_file = snakemake.input.rdp
meta_file = snakemake.input.meta
out = snakemake.output[0]

df = pd.read_csv(csv_file)
recombinants = set(df["Recombinant"].dropna())

meta = pd.read_csv(meta_file, sep="\t")
meta["is_recombinant"] = meta["accession"].isin(recombinants)

meta.to_csv(out, sep="\t", index=False)
