import pandas as pd

meta = pd.read_csv(snakemake.input.recomb, sep="\t")
g2p = pd.read_csv(snakemake.input.g2p, sep="\t")

df = meta.merge(g2p, left_on="accession", right_on="id", how="left")

df.to_csv(snakemake.output[0], sep="\t", index=False)

