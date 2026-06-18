import csv
import os

configfile: "config/config.yaml"

# ---------------------------------------------------------------------------
# Load project list (stdlib only — no numpy/pandas required)
# config["projects"] overrides projects_file when set to a non-null list
# ---------------------------------------------------------------------------
_all_projects = []
with open(config["projects_file"]) as _fh:
    for row in csv.DictReader(_fh, delimiter="\t"):
        _all_projects.append(row["project_id"])

_subset = config.get("projects") or []
PROJECTS = [p for p in _all_projects if p in _subset] if _subset else _all_projects

RAW       = config["dirs"]["raw"]
PROCESSED = config["dirs"]["processed"]
SUMMARY   = config["dirs"]["summary"]
MERGED    = config["dirs"]["merged"]
ANNOT     = config["dirs"]["annotation"]

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
        targets.append(f"{ANNOT}/gene_annotation.tsv")

    if config["modalities"]["run_mirna"]:
        targets += expand(f"{RAW}/{{project}}/mirna.tsv", project=PROJECTS)
        targets.append(f"{MERGED}/mirna.tsv")
        targets.append(f"{ANNOT}/mirna_annotation.tsv")

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
