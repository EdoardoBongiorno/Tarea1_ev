rule run_enrichr:
    input:
        all_genes="results/{subset}/deseq2_results.tsv"
    output:
        table="results/{subset}/enrichr_{direction}_{library}.tsv"
    params:
        alpha=config["deseq2"]["alpha"],
        lfc_threshold=config["deseq2"]["lfc_threshold"],
        library_name=lambda wildcards: config["enrichment"]["libraries"][wildcards.library]
    wildcard_constraints:
        subset=r"full|n[0-9]+_s[0-9]+",
        direction=r"up|down",
        library=r"go_bp|reactome"
    retries: 2
    log:
        "logs/enrichr_{subset}_{direction}_{library}.log"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/run_enrichr.py"