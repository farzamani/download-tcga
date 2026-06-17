import csv

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Load project list (stdlib only — no numpy/pandas required)
# ---------------------------------------------------------------------------
PROJECTS = []
with open(config["projects_file"]) as _fh:
    for row in csv.DictReader(_fh, delimiter="\t"):
        PROJECTS.append(row["project_id"])

RAW       = config["dirs"]["raw"]
PROCESSED = config["dirs"]["processed"]
SUMMARY   = config["dirs"]["summary"]
MERGED    = config["dirs"]["merged"]

# ---------------------------------------------------------------------------
# Rule modules
# ---------------------------------------------------------------------------
include: "rules/download.smk"
include: "rules/harmonize.smk"
include: "rules/summary.smk"
include: "rules/merge.smk"

# ---------------------------------------------------------------------------
# Helper: collect all expected outputs based on config toggles
# ---------------------------------------------------------------------------
def all_targets():
    targets = []

    if config["modalities"]["run_rna"]:
        targets += expand(f"{RAW}/{{project}}/rna.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/rna.tsv")

    if config["modalities"]["run_mirna"]:
        targets += expand(f"{RAW}/{{project}}/mirna.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/mirna.tsv")

    if config["modalities"]["run_methylation"]:
        targets += expand(f"{RAW}/{{project}}/methylation.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/methylation.tsv")

    if config["modalities"]["run_cnv"]:
        targets += expand(f"{RAW}/{{project}}/cnv.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/cnv.tsv")

    if config["modalities"]["run_annotation"]:
        targets += expand(f"{PROCESSED}/{{project}}/annotation.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/annotation.tsv")

    # Per-project summary (always produced)
    targets += expand(f"{SUMMARY}/{{project}}.modality_summary.tsv", project=PROJECTS)

    return targets


# ---------------------------------------------------------------------------
# Default target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        all_targets()
