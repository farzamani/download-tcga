# Column Reference

## annotation.tsv (per-project and merged)

One row per GDC sample aliquot.  Clinical columns are patient-level and are repeated across all aliquots from the same patient.

| Column | Type | Description |
|--------|------|-------------|
| `barcode` | string | Full TCGA aliquot barcode (e.g. `TCGA-A1-A0SO-01A-11R-A084-07`). Unique per aliquot. |
| `patient_id` | string | First 12 characters of the barcode — uniquely identifies a patient across all TCGA projects. |
| `sample_id` | string | First 16 characters of the barcode — identifies a tissue sample (one patient may have multiple: tumor, normal, metastatic). |
| `project` | string | TCGA project identifier, e.g. `TCGA-BRCA`. |
| `sample_type` | string | Human-readable GDC sample type label, e.g. `Primary Tumor`, `Solid Tissue Normal`, `Blood Derived Normal`. |
| `sample_type_code` | string | Two-digit numeric code from barcode positions 14–15. Common values: `01` = Primary Solid Tumor, `06` = Metastatic, `10` = Blood Derived Normal, `11` = Solid Tissue Normal. Full list at https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes |
| `subtype` | string | Molecular subtype from the TCGA Pan-Cancer Atlas (`TCGAquery_subtype`). Available for most major cancer types; `NA` for projects with no curated subtype. Example values for BRCA: `BRCA_Basal`, `BRCA_LumA`, `BRCA_LumB`, `BRCA_Her2`, `BRCA_Normal`. |

### Tumor purity columns

Tumor purity is the fraction of a bulk tumor sample that consists of cancer cells (as opposed to stromal, immune, or other normal cells).  Five independent estimation methods are included; they often disagree because they use different data modalities.

| Column | Source | Description |
|--------|--------|-------------|
| `tumor_purity` | Consensus (CPE) | Consensus Purity Estimate — median of all available method-specific estimates for that sample. Use this as the primary purity column. Range 0–1. |
| `ABSOLUTE_purity` | ABSOLUTE | Estimates purity and ploidy jointly from somatic copy-number alterations and loss-of-heterozygosity. Tends to be the most accurate when sufficient copy-number signal is present. |
| `ESTIMATE_purity` | ESTIMATE | Derived from gene-expression profiles of immune and stromal cell signatures. A high immune or stromal infiltration score pushes this estimate lower. Does not require DNA data. |
| `LUMP_purity` | LUMP | Leukocyte Unmethylation for Purity — uses DNA methylation at leukocyte-specific CpG loci as a proxy for immune contamination. |
| `IHC_purity` | IHC | Pathologist estimate obtained from immunohistochemistry slides. The most direct measurement but the least available (not all samples have slides). |
| `TIMER_purity` | TIMER | Tumor Immune Estimation Resource — deconvolves bulk gene expression into immune-cell fractions and derives purity as the residual. |

> **Note on "stromal / immune percentage"**: the ESTIMATE algorithm internally computes a `StromalScore` and `ImmuneScore` (higher = more infiltration) before deriving purity.  These intermediate scores are not stored in the TCGA tumor-purity consensus dataset used here; `ESTIMATE_purity` is the final purity value from that algorithm.  If you need the raw stromal/immune scores, run `TCGAanalyze_estimate()` on your RNA expression matrix.

### Clinical columns

All clinical data is patient-level (from GDC clinical XML / `GDCquery_clinic`).

| Column | Type | Description |
|--------|------|-------------|
| `gender` | string | Patient biological sex: `male` or `female`. |
| `age_at_diagnosis` | integer | Age in years at the time of initial pathologic diagnosis. |
| `tumor_stage` | string | Pathologic stage at diagnosis (AJCC system), e.g. `Stage I`, `Stage IIA`, `Stage IV`. Not all projects use identical staging conventions. |
| `vital_status` | string | `Alive` or `Dead` at last follow-up. |
| `days_to_death` | integer | Days from initial diagnosis to death. `NA` for patients who were alive at last follow-up. |
| `days_to_last_follow_up` | integer | Days from initial diagnosis to the last recorded follow-up contact. Use this together with `vital_status` for survival analysis. |

---

## Joining data files to annotation

All data files use `barcode` as the sole identifier.  To attach clinical or purity metadata, join on `barcode` to `annotation.tsv`:

```r
library(data.table)
rna  <- fread("results/merged/rna.tsv")
ann  <- fread("results/merged/annotation.tsv")
data <- ann[rna, on = "barcode"]   # left join: annotation columns prepended
```

---

## rna.tsv (per-project and merged)

One row per sample, one column per gene.  Join to `annotation.tsv` on `barcode` for metadata.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `ENSG...` | Raw unstranded read counts (STAR aligner, GDC harmonised pipeline) |

---

## mirna.tsv (per-project and merged)

One row per sample, one column per precursor miRNA.  Join to `annotation.tsv` on `barcode` for metadata.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `hsa-mir-...` | Raw read counts (BCGSC miRNA Profiling pipeline) |

---

## methylation.tsv (per-project and merged)

One row per sample, one column per CpG probe.  Only the top `max_cpgs` most-variable probes are retained (configurable; default 5,000).  Join to `annotation.tsv` on `barcode` for metadata.

| Column | Description |
|--------|-------------|
| `barcode` | Full TCGA aliquot barcode (join key) |
| `cgXXXXXXXX` | Illumina 450k/27k probe ID. Beta value = methylated / (methylated + unmethylated + 100). Ranges 0–1; 0 = fully unmethylated, 1 = fully methylated. |

---

## cnv.tsv (per-project and merged)

Long format — one row per copy-number segment per sample.  Join to `annotation.tsv` on `barcode` for metadata.

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

One row per protein-coding Ensembl gene (GRCh38, Ensembl release queried at run time via biomaRt).

| Column | Description |
|--------|-------------|
| `gene_id` | Ensembl gene ID (`ENSG...`) — matches columns in `rna.tsv` |
| `gene_name` | HGNC symbol (e.g. `TP53`, `BRCA1`) |
| `gene_type` | Biotype — always `protein_coding` in this table |
| `chromosome` | Chromosome (1–22, X, Y, MT) |
| `start` | Transcription start position (GRCh38, 1-based) |
| `end` | Transcription end position |
| `strand` | `+` or `-` |

---

## mirna_annotation.tsv

One row per precursor miRNA (hsa, miRBase v22).

| Column | Description |
|--------|-------------|
| `mirna_id` | Precursor miRNA name (e.g. `hsa-mir-21`) — matches columns in `mirna.tsv` |
| `accession` | miRBase accession for the precursor (prefix `MI0`) |
| `mature_id_5p` | Name of the 5p mature product (e.g. `hsa-miR-21-5p`), if annotated |
| `mature_acc_5p` | miRBase accession for the 5p mature product (prefix `MIMAT`) |
| `mature_id_3p` | Name of the 3p mature product (e.g. `hsa-miR-21-3p`), if annotated |
| `mature_acc_3p` | miRBase accession for the 3p mature product (prefix `MIMAT`) |
