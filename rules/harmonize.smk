rule build_annotation:
    output:
        tsv = f"{PROCESSED}/{{project}}/annotation.tsv"
    log:
        "logs/build_annotation/{project}.log"
    conda:
        "envs/r-tcgabiolinks.yaml"
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30
    params:
        sample_type = config["sample_type"]
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/build_annotation.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            > {log} 2>&1
        """
