# Design Notes: download-tcga

## Data model

### Sample identifier hierarchy

TCGA barcodes encode a hierarchy of identifiers within a single string:

```
TCGA-02-0001-01C-01D-0182-06
│    │    │   │   │   │    └─ plate, portion, analyte, ...
│    │    │   │   └───────── vial
│    │    │   └───────────── sample type code (01 = Primary Solid Tumor)
│    │    └───────────────── patient number
│    └────────────────────── tissue source site
└─────────────────────────── program
```

The workflow uses three levels:

| Field | Characters | Example | Notes |
|---|---|---|---|
| `patient_id` | 1–12 | `TCGA-02-0001` | One patient; may have multiple samples |
| `sample_id` | 1–16 | `TCGA-02-0001-01` | One vial/sample type per patient |
| `barcode` | full | `TCGA-02-0001-01C-01D-0182-06` | Unique per aliquot |

**When to join on which field:**

- `patient_id`: use when matching across modalities that may have been assayed
  from different aliquots of the same tumor (e.g., RNA-seq from one aliquot,
  methylation from another). This is the broadest join and loses sample-type
  specificity.
- `sample_id`: use when you want to verify that the same physical sample
  (same tumor, same vial) was assayed across modalities. Tighter, may drop
  patients assayed on only one platform.
- `barcode`: use only when aliquot-level traceability is required (e.g., batch
  correction using plate information).

The annotation table is keyed on `patient_id` (clinical data is patient-level).
Expression, methylation, and CNV TSVs carry all three fields so analysts can
choose the appropriate join key.

### CNV format

CNV output is **long-format** (one row per genomic segment), not wide. This is
intentional: segments vary in number across samples, making a wide matrix
impractical. Typical downstream steps (GISTIC2, CBS-based tools) expect
segment-level input.

Columns: `barcode`, `patient_id`, `sample_id`, `project`, `sample_type`,
`chromosome`, `start`, `end`, `num_probes`, `segment_mean`.

### Methylation format

Methylation output is **wide-format** (samples × CpGs) with beta values [0, 1].
The default cap of 50,000 CpGs (top variance) keeps file sizes manageable
while retaining the most informative sites. Set `max_cpgs: null` in
`config.yaml` to retain all probes (~450,000 for the 450k array).

---

## Matching strategy across modalities

Recommended approach for multi-omics integration:

1. Start with the annotation table (`annotation.tsv`) as the reference patient list.
2. Join each modality on `patient_id` (the 12-character prefix).
3. For samples where multiple aliquots exist (duplicates at `sample_id` level),
   keep the aliquot with the highest RNA-seq alignment rate or the lowest
   proportion of missing methylation values — this is project-specific and
   not automated in this workflow.
4. The summary file (`{project}.modality_summary.tsv`) reports pairwise
   patient-level overlap so you can see how many patients have complete data
   before running any integration step.

### Known ambiguity: TCGA-OV

TCGA-OV was assayed on both Illumina 450k and 27k arrays across different
sample batches. The methylation script queries 450k first and falls back to 27k.
If both are available, only 450k samples are returned (to keep probe sets
consistent). Mixing 450k and 27k beta values without renormalisation is not
recommended.

---

## Why Snakemake

- **Rule-based resumption**: Snakemake tracks outputs and only re-runs rules
  whose outputs are missing or stale. Long GDC downloads (methylation can
  take hours) survive interrupted runs without restarting from scratch.
- **Config-driven DAG**: enabling or disabling a modality via a boolean in
  `config.yaml` changes which rules are included in the DAG at parse time,
  with no code changes needed.
- **Conda integration**: each rule uses the same `r-tcgabiolinks.yaml`
  environment, pinning R and Bioconductor package versions without requiring
  a system-wide R installation.
- **Parallelism**: with `--cores N`, independent project × modality
  combinations run in parallel automatically.
- **Dry-run support**: `snakemake -n` lists all planned jobs without touching
  GDC, making it easy to validate config before submitting to a cluster.

---

## Known caveats

### Methylation availability

| Project | Array | Notes |
|---|---|---|
| TCGA-OV | 27k + 450k | Two batches; only 450k returned by default |
| TCGA-LAML | 450k | Blood tumour; "TP" sample type may return 0 samples — use "TB" |
| TCGA-DLBC | 450k | Small cohort (~48 samples) |
| TCGA-MESO | 450k | Very small cohort (~87 samples); methylation may be incomplete |

### CNV availability

- TCGA-LAML (acute myeloid leukaemia) uses `"TB"` (Blood Derived Normal) and
  `"TBM"` (Bone Marrow) sample types rather than `"TP"`. Setting `sample_type: "TP"`
  will return 0 CNV segments for LAML. Override per-project if needed.
- Some projects were processed with different SNP array platforms (Affymetrix SNP 6.0
  vs. Illumina Infinium), which affects segment resolution. The workflow does not
  normalise across platforms.

### miRNA availability

TCGA-GBM has a subset of samples with miRNA data; earlier collection periods used
a different sequencing protocol. Expect fewer miRNA samples than RNA-seq samples.

### Clinical data completeness

Not all clinical fields are available for all projects. Fields such as
`tumor_stage` may be missing for haematological malignancies. The annotation
script uses `NA` for missing values and does not impute.

### Sample type codes

The default `sample_type: "TP"` (Primary Solid Tumor) excludes:
- Blood-derived tumours (LAML, DLBC) — use `"TB"` or `"TBM"`
- Metastatic samples — use `"TM"`
- Recurrent tumours — use `"TR"`

To process multiple sample types, run separate workflow instances with different
`config.yaml` files or extend the Snakefile to iterate over a list of types.
