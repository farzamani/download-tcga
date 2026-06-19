# download-tcga

Snakemake workflow to download and harmonise TCGA multi-omics data using [TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks/). All data are open-access GDC files.

## Requirements

- Snakemake ≥ 7.0
- Conda / Mamba

## Quick start

```bash
# dry run
snakemake -n --use-conda

# local run
snakemake --use-conda --cores 4

# SLURM
snakemake --executor slurm --use-conda --cores 4 --jobs 100 \
    --default-resources slurm_account=<YOUR_ACCOUNT>
```

## Configuration

All settings live in `config/config.yaml`.

| Key | Default | Description |
|-----|---------|-------------|
| `projects` | `null` (all 33) | List of TCGA project IDs to download, or `null` to run all projects in `projects_file` |
| `projects_file` | `config/projects.tsv` | Full list of 33 TCGA projects |
| `sample_type` | `all` | `all` downloads every sample type; restrict with a code e.g. `TP` (Primary Tumor), `NT` (Solid Tissue Normal) |
| `max_cpgs` | `5000` | Top-N most-variable CpG probes to retain per project. Set to `null` to keep all (very large) |
| `modalities.run_mrna` | `true` | Toggle each modality independently |
| `modalities.run_mirna` | `true` | |
| `modalities.run_methylation` | `true` | |
| `modalities.run_cnv` | `true` | |
| `modalities.run_annotation` | `true` | |

## Outputs

The GDC download cache (`results/cache/`) is removed automatically on successful completion. All other outputs are retained.

### Per-project files

| Path | Description |
|------|-------------|
| `results/raw/{project}/mrna.tsv` | Samples × protein-coding genes, STAR unstranded counts (ENSG IDs, no version suffix) |
| `results/raw/{project}/mirna.tsv` | Samples × mature miRNA arms (hsa-miR-9-5p / hsa-miR-9-3p), raw counts |
| `results/raw/{project}/methylation.tsv` | Samples × top-N CpG probes, Illumina 450k/27k beta values |
| `results/raw/{project}/cnv.tsv` | Copy-number segments, long format |
| `results/annotation/{project}/annotation.tsv` | One row per sample: clinical, Pan-Cancer Atlas subtype, tumor purity (4 methods + consensus) |
| `results/summary/{project}.modality_summary.tsv` | Per-modality sample counts and patient-level intersections (2–4 way) |

### Pan-TCGA merged files

| Path | Description |
|------|-------------|
| `results/merged/mrna.tsv` | All projects stacked; feature columns are the intersection across projects |
| `results/merged/mirna.tsv` | All projects stacked |
| `results/merged/methylation.tsv` | All projects stacked; feature columns are the intersection across projects |
| `results/merged/cnv.tsv` | All projects stacked (long format) |
| `results/merged/annotation.tsv` | All projects stacked |

### Annotation dictionaries

| Path | Description |
|------|-------------|
| `results/annotation/gene_annotation.tsv` | Ensembl gene ID → HGNC symbol, coordinates (protein-coding genes only) |
| `results/annotation/mirna_annotation.tsv` | Precursor → 5p/3p mature arm names and miRBase v22 accessions |

### File format

All data files use `barcode` (full TCGA aliquot barcode) as the sole identifier column — all other metadata (patient ID, sample type, clinical, purity) is in `annotation.tsv` and joined on `barcode`:

```r
library(data.table)
mrna <- fread("results/merged/mrna.tsv")
ann  <- fread("results/merged/annotation.tsv")
data <- ann[mrna, on = "barcode"]
```

See `docs/columns.md` for a full description of every column in every output file.

## License

MIT
