import sys, pandas as pd

EXPECTED = [
    "sample","virus","genotype","subgenotype",
    "gene","drug","mutation","detected","type"
]

df = pd.read_csv(sys.argv[1], sep="\t")
missing = set(EXPECTED) - set(df.columns)

if missing:
    raise SystemExit(f"Missing columns: {missing}")

