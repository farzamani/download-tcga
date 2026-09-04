# Merge rules — combine per-project TSVs into one pan-TCGA file per modality.
# Wide-format modalities (rna, mirna, methylation): features are intersected so
# the merged matrix contains only probes/genes present in every project.
# Long-format modalities (cnv, annotation): simple row-bind, no feature alignment needed.

rule merge_mrna:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{RAW}/{{project}}/mrna.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/mrna.tsv"
    log:
        "logs/merge/mrna.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 60
    shell:
        # --max-ppsize raises R's PROTECT stack (default 50000); rbindlist()
        # on ~60k-column wide tables (full transcriptome) overflows it otherwise.
        """
        mkdir -p $(dirname {output.tsv})
        Rscript --max-ppsize=500000 {input.script} wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_mirna:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{RAW}/{{project}}/mirna.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/mirna.tsv"
    log:
        "logs/merge/mirna.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript --max-ppsize=500000 {input.script} wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_methylation:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{RAW}/{{project}}/methylation.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/methylation.tsv"
    log:
        "logs/merge/methylation.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    params:
        top_var_flag = (f"--top-var={config['merged_max_cpgs']}"
                         if config.get("merged_max_cpgs") else "")
    threads: 2
    resources:
        mem_mb  = 64000,   # 33 projects × up to 20k CpGs is large
        runtime = 120
    shell:
        # --top-var keeps only the N common CpGs with highest pooled variance,
        # since the per-project cap (max_cpgs) alone leaves each project with
        # a different cancer-type-specific top set.
        """
        mkdir -p $(dirname {output.tsv})
        Rscript --max-ppsize=500000 {input.script} wide {output.tsv} \
            {params.top_var_flag} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_cnv_gene:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{RAW}/{{project}}/cnv_gene.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/cnv_gene.tsv"
    log:
        "logs/merge/cnv_gene.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 32000,
        runtime = 60
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript --max-ppsize=500000 {input.script} wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_cnv:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{RAW}/{{project}}/cnv.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/cnv.tsv"
    log:
        "logs/merge/cnv.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 16000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript {input.script} long {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_annotation:
    input:
        script = "scripts/merge_modality.R",
        tsvs   = expand(f"{ANNOT}/{{project}}/annotation.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/annotation.tsv"
    log:
        "logs/merge/annotation.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 15
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript {input.script} long {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """
