# download-tcga

Snakemake workflow to download and harmonise TCGA multi-omics data (mRNA, miRNA, methylation, CNV, annotation) using [TCGAbiolinks](https://bioconductor.org/packages/TCGAbiolinks/). All data are open-access GDC files.

This repo was built using Claude Code.

## Requirements

- Snakemake ≥ 7.0
- Conda / Mamba

## Usage

```bash
# dry run
snakemake -n

# full run
snakemake --use-conda --cores 4
```

For SLURM, pass resources via `--cluster` or a profile:

```bash
snakemake --use-conda --cores 4 \
  --cluster "sbatch --mem={resources.mem_mb}M --time={resources.runtime} --cpus-per-task={threads}"
```

## Configuration

Edit `config/config.yaml` to toggle modalities, set sample type, or change output directories. Edit `config/projects.tsv` to restrict to specific TCGA projects (all 33 are included by default).

## Outputs

Per project in `results/raw/{project}/` and `results/processed/{project}/`, plus pan-TCGA merged files in `results/merged/`:

| File | Description |
|---|---|
| `rna.tsv` | Samples × genes STAR counts |
| `mirna.tsv` | Samples × miRNAs read counts |
| `methylation.tsv` | Samples × CpGs beta values (top 50k by variance) |
| `cnv.tsv` | Copy number segments (long format) |
| `annotation.tsv` | Clinical covariates, subtype, tumor purity |
| `results/summary/{project}.modality_summary.tsv` | Sample counts and cross-modality overlap |

All TSVs share a common first five columns: `barcode`, `patient_id`, `sample_id`, `project`, `sample_type`.

See `docs/design.md` for data model details and known caveats per project.

## License

MIT
