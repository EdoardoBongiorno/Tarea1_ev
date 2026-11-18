import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

with open(snakemake.log[0], "w") as log_f:
    sys.stderr = log_f
    sys.stdout = log_f

    df = pd.read_csv(snakemake.input.all_genes, sep="\t").dropna(subset=["padj", "log2FoldChange"])
    alpha = float(snakemake.params.alpha)
    lfc_thresh = float(snakemake.params.lfc_threshold)

    df["status"] = "NS"
    df.loc[(df["padj"] < alpha) & (df["log2FoldChange"] > lfc_thresh), "status"] = "Up"
    df.loc[(df["padj"] < alpha) & (df["log2FoldChange"] < -lfc_thresh), "status"] = "Down"

    colors = {"Up": "#d62728", "Down": "#1f77b4", "NS": "#bdbdbd"}

    fig, ax = plt.subplots(figsize=(7, 5))
    for status, color in colors.items():
        sub = df[df["status"] == status]
        ax.scatter(
            sub["log2FoldChange"],
            -np.log10(sub["padj"].clip(lower=1e-300)),
            c=color,
            label=f"{status} ({len(sub)})",
            alpha=0.6,
            s=12,
            edgecolors="none",
        )
    ax.axvline(lfc_thresh, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_thresh, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(alpha), color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change (Tumor vs Normal)")
    ax.set_ylabel("-log10(padj)")
    ax.set_title(f"Volcano Plot ({snakemake.wildcards.subset})")
    ax.legend()
    plt.tight_layout()
    plt.savefig(snakemake.output.volcano, dpi=300)
    plt.close()

    fig, ax = plt.subplots(figsize=(7, 5))
    for status, color in colors.items():
        sub = df[df["status"] == status]
        ax.scatter(
            np.log10(sub["baseMean"] + 1),
            sub["log2FoldChange"],
            c=color,
            label=f"{status} ({len(sub)})",
            alpha=0.6,
            s=12,
            edgecolors="none",
        )
    ax.axhline(lfc_thresh, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-lfc_thresh, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xlabel("log10(baseMean + 1)")
    ax.set_ylabel("log2 Fold Change (Tumor vs Normal)")
    ax.set_title(f"MA Plot ({snakemake.wildcards.subset})")
    ax.legend()
    plt.tight_layout()
    plt.savefig(snakemake.output.ma, dpi=300)
    plt.close()