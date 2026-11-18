rule prepare_full:
    input:
        counts=config["raw_data"]["counts"],
        metadata=config["raw_data"]["metadata"]
    output:
        samples="results/full/samples.tsv"
    log:
        "logs/prepare_full.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/prepare_full.py"


rule run_deseq2:
    input:
        counts=config["raw_data"]["counts"],
        samples="results/{subset}/samples.tsv"
    output:
        all_genes="results/{subset}/deseq2_results.tsv",
        de_genes="results/{subset}/de_genes.tsv"
    params:
        min_count=config["filtering"]["min_count"],
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
        design=config["deseq2"]["design"]
    wildcard_constraints:
        subset=r"full|n[0-9]+_s[0-9]+"
    threads: 4
    log:
        "logs/deseq2_{subset}.log"
    benchmark:
        "benchmarks/deseq2_{subset}.tsv"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/run_deseq2.py"


rule plot_de:
    input:
        all_genes="results/{subset}/deseq2_results.tsv"
    output:
        volcano="results/{subset}/volcano_plot.png",
        ma="results/{subset}/ma_plot.png"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"]
    wildcard_constraints:
        subset=r"full|n[0-9]+_s[0-9]+"
    log:
        "logs/plot_de_{subset}.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/plot_de.py"