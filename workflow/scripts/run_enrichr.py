import sys
import pandas as pd
import gseapy as gp

with open(snakemake.log[0], "w") as log_f:
    sys.stderr = log_f
    sys.stdout = log_f

    all_df = pd.read_csv(snakemake.input.all_genes, sep="\t").dropna(subset=["padj", "log2FoldChange"])
    alpha = float(snakemake.params.alpha)
    lfc_thresh = float(snakemake.params.lfc_threshold)
    direction = snakemake.wildcards.direction
    library_name = snakemake.params.library_name

    # Background esplicito: tutti i geni che hanno superato il filtro di espressione
    background_genes = all_df["gene"].astype(str).tolist()

    # Selezione geni up o down regolati
    if direction == "up":
        gene_list = all_df.loc[
            (all_df["padj"] < alpha) & (all_df["log2FoldChange"] > lfc_thresh), "gene"
        ].astype(str).tolist()
    elif direction == "down":
        gene_list = all_df.loc[
            (all_df["padj"] < alpha) & (all_df["log2FoldChange"] < -lfc_thresh), "gene"
        ].astype(str).tolist()
    else:
        raise ValueError(f"Direzione non valida: {direction}")

    print(f"Subset: {snakemake.wildcards.subset} | {direction} | {library_name} | Geni query: {len(gene_list)} | Background: {len(background_genes)}")

    if len(gene_list) == 0:
        pd.DataFrame().to_csv(snakemake.output.table, sep="\t", index=False)
        sys.exit(0)

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=[library_name],
        organism="human",
        background=background_genes,
        outdir=None,
    )

    res = enr.results
    res.to_csv(snakemake.output.table, sep="\t", index=False)
    print(f"Termini trovati: {len(res)}")