# Column Reference

## annotation.tsv

Paths: `results/annotation/{project}/annotation.tsv` (per-project) and `results/merged/annotation.tsv` (pan-TCGA).

One row per GDC sample aliquot.  Clinical columns are patient-level and are repeated across all aliquots from the same patient.

| Column | Type | Description |
|--------|------|-------------|
| `barcode` | string | Full TCGA aliquot barcode (e.g. `TCGA-A1-A0SO-01A-11R-A084-07`). Unique per aliquot. Join key for all data files. |
| `patient_id` | string | First 12 characters of the barcode — uniquely identifies a patient across all TCGA projects. |
| `sample_id` | string | First 16 characters of the barcode — identifies a tissue sample (one patient may have multiple: tumor, normal, metastatic). |
| `project` | string | TCGA project identifier, e.g. `TCGA-BRCA`. |
| `sample_type` | string | Human-readable GDC sample type label, e.g. `Primary Tumor`, `Solid Tissue Normal`, `Blood Derived Normal`. |
| `sample_type_code` | string | Two-digit numeric code from barcode positions 14–15. Common values: `01` = Primary Solid Tumor, `06` = Metastatic, `10` = Blood Derived Normal, `11` = Solid Tissue Normal. Full list at https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes |
| `subtype` | string | Pan-Cancer Atlas consensus subtype (`Subtype_Selected` from `PanCancerAtlas_subtypes()`). `NA` for projects not covered by the Atlas. Example values for BRCA: `BRCA.Basal`, `BRCA.LumA`, `BRCA.LumB`, `BRCA.Her2`, `BRCA.Normal`. |

### Tumor purity columns

Tumor purity is the fraction of a bulk tumor sample that consists of cancer cells (as opposed to stromal, immune, or other normal cells).  From the TCGAbiolinks `Tumor.purity` dataset (Aran et al. 2015).

| Column | Source | Description |
|--------|--------|-------------|
| `tumor_purity` | Consensus (CPE) | Consensus Purity Estimate — median of available method-specific estimates. Use this as the primary purity column. Range 0–1. |
| `ESTIMATE_purity` | ESTIMATE | Derived from gene-expression profiles of immune and stromal cell signatures. A high immune or stromal infiltration score pushes this estimate lower. |
| `ABSOLUTE_purity` | ABSOLUTE | Estimates purity and ploidy jointly from somatic copy-number alterations and loss-of-heterozygosity. |
| `LUMP_purity` | LUMP | Leukocyte Unmethylation for Purity — uses DNA methylation at leukocyte-specific CpG loci as a proxy for immune contamination. |
| `IHC_purity` | IHC | Pathologist estimate from immunohistochemistry slides. Most direct measurement; not available for all samples. |

> **Note**: The ESTIMATE algorithm also computes `StromalScore` and `ImmuneScore` (higher = more infiltration) before deriving purity. These intermediate scores are not in the `Tumor.purity` bundled dataset. If you need them, run `TCGAanalyze_estimate()` on the mRNA expression matrix.

### Clinical columns

Patient-level data from GDC clinical XML (`GDCquery_clinic`). Repeated for every sample from the same patient.

| Column | Type | Description |
|--------|------|-------------|
| `gender` | string | Patient biological sex: `male` or `female`. |
| `age_at_diagnosis` | integer | Age in years at the time of initial pathologic diagnosis. |
| `tumor_stage` | string | Pathologic stage (AJCC), e.g. `Stage I`, `Stage IIA`, `Stage IV`. Conventions vary by project. |
| `vital_status` | string | `Alive` or `Dead` at last follow-up. |
| `days_to_death` | integer | Days from diagnosis to death. `NA` for patients still alive. |
| `days_to_last_follow_up` | integer | Days from diagnosis to last recorded contact. Use with `vital_status` for survival analysis. |

---

## Joining data files to annotation

All data files use `barcode` as the sole identifier column.  Attach clinical or purity metadata by joining on `barcode`:

```r
library(data.table)
mrna <- fread("results/merged/mrna.tsv")
ann  <- fread("results/merged/annotation.tsv")
data <- ann[mrna, on = "barcode"]   # left join: annotation columns prepended
```

---

## mrna.tsv (per-project and merged)

One row per sample, one column per protein-coding gene.  Gene IDs match `results/annotation/gene_annotation.tsv`.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `ENSG...` | Raw unstranded read counts (STAR aligner, GDC harmonised pipeline). Ensembl version suffix stripped (e.g. `ENSG00000000003`, not `ENSG00000000003.15`). |

---

## mirna.tsv (per-project and merged)

One row per sample, one column per **mature miRNA arm**.  Where multiple precursor loci produce the same arm (e.g. `hsa-mir-9-1/2/3` all produce `hsa-miR-9-5p`), their counts are summed.  Precursors with no arm annotation in miRBase v22 retain their original precursor ID.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `hsa-miR-...-5p` | Raw read count for the 5p arm (summed across contributing loci) |
| `hsa-miR-...-3p` | Raw read count for the 3p arm (summed across contributing loci) |
| `hsa-mir-...` | Precursor ID retained when no arm annotation exists in miRBase v22 |

---

## methylation.tsv (per-project and merged)

One row per sample, one column per CpG probe.  Only the top `max_cpgs` most-variable probes are retained (default 5,000).

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `cgXXXXXXXX` | Illumina 450k/27k probe ID. Beta value = methylated / (methylated + unmethylated + 100). Range 0–1; 0 = fully unmethylated, 1 = fully methylated. |

---

## cnv.tsv (per-project and merged)

Long format — one row per copy-number segment per sample.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `chromosome` | Chromosome of the segment |
| `start` | Genomic start position (1-based) |
| `end` | Genomic end position (inclusive) |
| `num_probes` | Number of SNP array probes in the segment |
| `segment_mean` | Log2 copy-number ratio (tumor / normal). 0 ≈ diploid; positive = amplification; negative = deletion. |

---

## gene_annotation.tsv

Path: `results/annotation/gene_annotation.tsv`

One row per protein-coding Ensembl gene (GRCh38, queried from Ensembl at run time via biomaRt).

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID (`ENSG...`, no version suffix) — matches columns in `mrna.tsv` |
| `gene_name` | HGNC symbol (e.g. `TP53`, `BRCA1`) |
| `gene_type` | Biotype — always `protein_coding` in this table |
| `chromosome` | Chromosome (1–22, X, Y, MT) |
| `start` | Transcription start position (GRCh38, 1-based) |
| `end` | Transcription end position |
| `strand` | `+` or `-` |

---

## mirna_annotation.tsv

Path: `results/annotation/mirna_annotation.tsv`

One row per precursor miRNA (human, miRBase v22).  The `mirna_id` column links to `mirna.tsv` — either directly (for unannotated precursors) or via the mature arm names (for annotated precursors).

| Column | Description |
|--------|-------------|
| `mirna_id` | Precursor miRNA name (e.g. `hsa-mir-21`) |
| `accession` | miRBase accession for the precursor (prefix `MI0`) |
| `mature_id_5p` | Name of the 5p mature arm (e.g. `hsa-miR-21-5p`) — matches columns in `mirna.tsv` |
| `mature_acc_5p` | miRBase accession for the 5p arm (prefix `MIMAT`) |
| `mature_id_3p` | Name of the 3p mature arm (e.g. `hsa-miR-21-3p`) — matches columns in `mirna.tsv` |
| `mature_acc_3p` | miRBase accession for the 3p arm (prefix `MIMAT`) |
