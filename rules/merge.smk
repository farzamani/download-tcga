# Merge rules — combine per-project TSVs into one pan-TCGA file per modality.
# Wide-format modalities (rna, mirna, methylation): features are intersected so
# the merged matrix contains only probes/genes present in every project.
# Long-format modalities (cnv, annotation): simple row-bind, no feature alignment needed.

rule merge_mrna:
    input:
        tsvs = expand(f"{RAW}/{{project}}/mrna.tsv", project=PROJECTS)
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
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/merge_modality.R wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_mirna:
    input:
        tsvs = expand(f"{RAW}/{{project}}/mirna.tsv", project=PROJECTS)
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
        Rscript scripts/merge_modality.R wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_methylation:
    input:
        tsvs = expand(f"{RAW}/{{project}}/methylation.tsv", project=PROJECTS)
    output:
        tsv = f"{MERGED}/methylation.tsv"
    log:
        "logs/merge/methylation.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 2
    resources:
        mem_mb  = 64000,   # 33 projects × up to 50k CpGs is large
        runtime = 120
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/merge_modality.R wide {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_cnv:
    input:
        tsvs = expand(f"{RAW}/{{project}}/cnv.tsv", project=PROJECTS)
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
        Rscript scripts/merge_modality.R long {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """


rule merge_annotation:
    input:
        tsvs = expand(f"{ANNOT}/{{project}}/annotation.tsv", project=PROJECTS)
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
        Rscript scripts/merge_modality.R long {output.tsv} {input.tsvs} \
            > {log} 2>&1
        """
