def _summary_inputs(wildcards):
    inputs = {}
    if config["modalities"]["run_mrna"]:
        inputs["mrna"] = f"{RAW}/{wildcards.project}/mrna.tsv"
    if config["modalities"]["run_mirna"]:
        inputs["mirna"] = f"{RAW}/{wildcards.project}/mirna.tsv"
    if config["modalities"]["run_methylation"]:
        inputs["methylation"] = f"{RAW}/{wildcards.project}/methylation.tsv"
    if config["modalities"]["run_cnv"]:
        inputs["cnv"] = f"{RAW}/{wildcards.project}/cnv.tsv"
    return inputs


rule summarize_project:
    input:
        unpack(_summary_inputs)
    output:
        tsv = f"{SUMMARY}/{{project}}.modality_summary.tsv"
    log:
        "logs/summarize/{project}.log"
    conda:
        os.path.join(workflow.basedir, "envs/r-tcgabiolinks.yaml")
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 15
    params:
        raw_dir = lambda wildcards: f"{RAW}/{wildcards.project}",
        project = lambda wildcards: wildcards.project
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        Rscript scripts/summarize_modalities.R \
            {params.project} \
            {params.raw_dir} \
            {output.tsv} \
            > {log} 2>&1
        """
