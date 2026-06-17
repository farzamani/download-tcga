# Download rules — one per omics modality.
# Each rule calls a dedicated R script that handles GDCquery/GDCdownload/GDCprepare
# and writes a TSV to results/raw/{project}/.
# On any error the R script writes an empty (header-only) TSV so that downstream
# rules are not blocked and the summary step reports 0 samples.

rule download_rna:
    output:
        tsv = f"{RAW}/{{project}}/rna.tsv"
    log:
        "logs/download_rna/{project}.log"
    conda:
        "envs/r-tcgabiolinks.yaml"
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 120      # minutes
    params:
        sample_type = config["sample_type"]
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/download_rna.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            > {log} 2>&1
        """


rule download_mirna:
    output:
        tsv = f"{RAW}/{{project}}/mirna.tsv"
    log:
        "logs/download_mirna/{project}.log"
    conda:
        "envs/r-tcgabiolinks.yaml"
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 60
    params:
        sample_type = config["sample_type"]
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/download_mirna.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            > {log} 2>&1
        """


rule download_methylation:
    output:
        tsv = f"{RAW}/{{project}}/methylation.tsv"
    log:
        "logs/download_methylation/{project}.log"
    conda:
        "envs/r-tcgabiolinks.yaml"
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 240      # methylation files are large
    params:
        sample_type = config["sample_type"],
        max_cpgs    = config.get("max_cpgs", 50000)
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/download_methylation.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.max_cpgs} \
            > {log} 2>&1
        """


rule download_cnv:
    output:
        tsv = f"{RAW}/{{project}}/cnv.tsv"
    log:
        "logs/download_cnv/{project}.log"
    conda:
        "envs/r-tcgabiolinks.yaml"
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 60
    params:
        sample_type = config["sample_type"]
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/download_cnv.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            > {log} 2>&1
        """
