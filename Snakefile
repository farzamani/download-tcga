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

RAW     = config["dirs"]["raw"]
SUMMARY = config["dirs"]["summary"]
MERGED  = config["dirs"]["merged"]
ANNOT   = config["dirs"]["annotation"]

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

    # Note: per-project results/raw/{project}/*.tsv are NOT listed as targets
    # here even though the rules exist — they're marked temp() in
    # download.smk and consumed by the merge/summarize rules below, which is
    # enough to build them via the DAG. Listing them as targets would force
    # Snakemake to always (re-)materialize them (defeating temp() cleanup)
    # since `rule all` would otherwise demand their continued existence.
    if config["modalities"]["run_mrna"]:
        targets.append(f"{MERGED}/mrna.tsv")
        targets.append(f"{ANNOT}/gene_annotation.tsv")

    if config["modalities"]["run_mirna"]:
        targets.append(f"{MERGED}/mirna.tsv")
        targets.append(f"{ANNOT}/mirna_annotation.tsv")

    if config["modalities"]["run_methylation"]:
        targets.append(f"{MERGED}/methylation.tsv")

    if config["modalities"].get("run_cnv_gene"):
        targets.append(f"{MERGED}/cnv_gene.tsv")

    if config["modalities"]["run_annotation"]:
        targets += expand(f"{ANNOT}/{{project}}/annotation.tsv", project=PROJECTS)
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


# ---------------------------------------------------------------------------
# Remove the GDC download cache and the per-project raw/ tree on successful
# completion. Final results live in merged/ and annotation/ — raw/ is a
# working area (per-project TSVs already marked temp() and consumed by the
# merge/summarize rules; this sweep also clears now-empty per-project
# subdirectories temp() cleanup leaves behind). On failure both are
# preserved so runs can resume without re-downloading.
# ---------------------------------------------------------------------------
onsuccess:
    import shutil
    cache_dir = config["dirs"]["gdc_cache"]
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
        print(f"Cleaned up GDC cache: {cache_dir}")
    if os.path.exists(RAW):
        shutil.rmtree(RAW)
        print(f"Cleaned up raw data directory: {RAW}")
