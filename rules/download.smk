rule download_rna:
    output:
        tsv = f"{RAW}/{{project}}/rna.tsv"
    log:
        "logs/download_rna/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 120
    params:
        sample_type = config["sample_type"],
        gdc_cache   = config["dirs"]["gdc_cache"]
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {params.gdc_cache}
        Rscript scripts/download_rna.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.gdc_cache} \
            > {log} 2>&1
        """


rule download_mirna:
    output:
        tsv = f"{RAW}/{{project}}/mirna.tsv"
    log:
        "logs/download_mirna/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 60
    params:
        sample_type = config["sample_type"],
        gdc_cache   = config["dirs"]["gdc_cache"]
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {params.gdc_cache}
        Rscript scripts/download_mirna.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.gdc_cache} \
            > {log} 2>&1
        """


rule download_methylation:
    output:
        tsv = f"{RAW}/{{project}}/methylation.tsv"
    log:
        "logs/download_methylation/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 240
    params:
        sample_type = config["sample_type"],
        gdc_cache   = config["dirs"]["gdc_cache"],
        max_cpgs    = config.get("max_cpgs", 50000)
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {params.gdc_cache}
        Rscript scripts/download_methylation.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.gdc_cache} \
            {params.max_cpgs} \
            > {log} 2>&1
        """


rule download_cnv:
    output:
        tsv = f"{RAW}/{{project}}/cnv.tsv"
    log:
        "logs/download_cnv/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 60
    params:
        sample_type = config["sample_type"],
        gdc_cache   = config["dirs"]["gdc_cache"]
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {params.gdc_cache}
        Rscript scripts/download_cnv.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.gdc_cache} \
            > {log} 2>&1
        """
