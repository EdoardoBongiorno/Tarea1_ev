import sys
import pandas as pd

with open(snakemake.log[0], "w") as log_f:
    sys.stderr = log_f
    sys.stdout = log_f

    counts = pd.read_csv(snakemake.input.counts, index_col=0)
    meta = pd.read_csv(snakemake.input.metadata, index_col=0)

    common = [s for s in counts.columns if s in meta.index]
    meta_full = meta.loc[common].copy()
    meta_full.index.name = "sample_id"

    meta_full.to_csv(snakemake.output.samples, sep="\t")
    print(f"Coorte completa salvata: {len(meta_full)} campioni.")