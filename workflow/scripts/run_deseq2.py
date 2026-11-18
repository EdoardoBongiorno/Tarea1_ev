import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

with open(snakemake.log[0], "w") as log_f:
    sys.stderr = log_f
    sys.stdout = log_f

    counts = pd.read_csv(snakemake.input.counts, index_col=0)
    meta = pd.read_csv(snakemake.input.samples, sep="\t", index_col=0)

    counts = counts[meta.index]

    min_count = int(snakemake.params.min_count)
    k = int(meta["site"].value_counts().min())
    keep = (counts >= min_count).sum(axis=1) >= k
    counts_filt = counts.loc[keep]

    print(f"Subset: {snakemake.wildcards.subset} | k={k} | Geni mantenuti: {len(counts_filt)}/{len(counts)}")

    design_formula = snakemake.params.design
    if "plate" in design_formula and meta["plate"].nunique() < 2:
        design_formula = "~ site"

    dds = DeseqDataSet(
        counts=counts_filt.T.astype(int),
        metadata=meta,
        design=design_formula,
        ref_level=["site", "Solid_Tissue_Normal"],
        n_cpus=snakemake.threads,
        quiet=True,
    )
    dds.deseq2()

    stat_res = DeseqStats(
        dds,
        contrast=["site", "Primary_Solid_Tumor", "Solid_Tissue_Normal"],
        alpha=float(snakemake.params.alpha),
        n_cpus=snakemake.threads,
        quiet=True,
    )
    stat_res.summary()

    res_df = stat_res.results_df.copy()
    res_df.index.name = "gene"
    res_df.to_csv(snakemake.output.all_genes, sep="\t")

    de_df = res_df[
        (res_df["padj"] < float(snakemake.params.alpha))
        & (res_df["log2FoldChange"].abs() > float(snakemake.params.lfc_threshold))
    ].copy()
    de_df.to_csv(snakemake.output.de_genes, sep="\t")
    print(f"Geni DE significativi trovati: {len(de_df)}")