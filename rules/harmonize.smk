# Per-project annotation: clinical, subtype, tumor purity
rule build_annotation:
    output:
        tsv = f"{ANNOT}/{{project}}/annotation.tsv"
    log:
        "logs/build_annotation/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/build_annotation.R \
            {wildcards.project} \
            {output.tsv} \
            > {log} 2>&1
        """


# Global gene annotation (runs once, not per project)
# Keyed on Ensembl gene ID — same IDs used as columns in mrna.tsv
rule annotate_genes:
    output:
        tsv = f"{ANNOT}/gene_annotation.tsv"
    log:
        "logs/annotate_genes.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/annotate_genes.R {output.tsv} > {log} 2>&1
        """


# Global miRNA annotation (runs once, not per project)
# Keyed on miRNA_ID — same IDs used as columns in mirna.tsv
rule annotate_mirna:
    output:
        tsv = f"{ANNOT}/mirna_annotation.tsv"
    log:
        "logs/annotate_mirna.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/annotate_mirna.R {output.tsv} > {log} 2>&1
        """
