rule download_mrna:
    output:
        tsv = f"{RAW}/{{project}}/mrna.tsv"
    log:
        "logs/download_mrna/{project}.log"
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
        Rscript scripts/download_mrna.R \
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
        runtime = 480      # methylation chunks are large; allow extra time for retries
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


# GDC's own gene-summarized copy number product, mapped/aligned to the same
# protein-coding Ensembl gene set as mrna.tsv (via gene_annotation.tsv), so
# CNV and expression describe identical genes for the shared latent space.
rule download_cnv_gene:
    input:
        gene_annotation = f"{ANNOT}/gene_annotation.tsv"
    output:
        tsv = f"{RAW}/{{project}}/cnv_gene.tsv"
    log:
        "logs/download_cnv_gene/{project}.log"
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
        Rscript scripts/download_cnv_gene.R \
            {wildcards.project} \
            {params.sample_type} \
            {output.tsv} \
            {params.gdc_cache} \
            {input.gene_annotation} \
            > {log} 2>&1
        """
