suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  outfile
# Produces one row per GDC sample (not per patient).  Actual sample barcodes
# and sample-type labels come from the GDC file manifest; clinical covariates
# from GDCquery_clinic are patient-level and are joined across all samples
# belonging to the same patient.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: build_annotation.R <project> <outfile>")

project <- args[1]
outfile <- normalizePath(args[2], mustWork = FALSE)

write_empty <- function(path) {
  fwrite(data.frame(
    barcode = character(), patient_id = character(), sample_id = character(),
    project = character(), sample_type = character(), sample_type_code = character(),
    subtype = character(),
    tumor_purity    = numeric(), ABSOLUTE_purity = numeric(),
    ESTIMATE_purity = numeric(), LUMP_purity     = numeric(),
    IHC_purity      = numeric(), TIMER_purity    = numeric(),
    gender = character(), age_at_diagnosis = numeric(), tumor_stage = character(),
    vital_status = character(), days_to_death = numeric(),
    days_to_last_follow_up = numeric()
  ), path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty annotation.tsv for ", project)
}

tryCatch({

  # ------------------------------------------------------------------
  # 1. Sample manifest — actual barcodes + sample-type labels from GDC
  # Try data types in preference order; first successful query is used.
  # ------------------------------------------------------------------
  manifest <- NULL
  for (spec in list(
    list(data.category = "Transcriptome Profiling",
         data.type     = "Gene Expression Quantification",
         workflow.type = "STAR - Counts"),
    list(data.category = "Copy Number Variation",
         data.type     = "Masked Copy Number Segment"),
    list(data.category = "Transcriptome Profiling",
         data.type     = "miRNA Expression Quantification",
         workflow.type = "BCGSC miRNA Profiling")
  )) {
    q <- tryCatch(do.call(GDCquery, c(list(project = project), spec)), error = function(e) NULL)
    if (is.null(q)) next
    r <- as.data.table(getResults(q))
    bc_col <- intersect(c("cases", "submitter_id"), colnames(r))[1]
    st_col <- intersect(c("sample_type", "cases.sample_type"), colnames(r))[1]
    if (is.na(bc_col)) next
    manifest <- r[, .(
      barcode     = get(bc_col),
      sample_type = if (is.na(st_col)) NA_character_ else get(st_col)
    )]
    manifest[, patient_id       := substr(barcode, 1, 12)]
    manifest[, sample_id        := substr(barcode, 1, 16)]
    manifest[, sample_type_code := substr(barcode, 14, 15)]
    manifest <- unique(manifest, by = "barcode")
    break
  }

  if (is.null(manifest) || nrow(manifest) == 0) {
    message("Could not retrieve sample manifest for ", project, "; writing empty annotation")
    write_empty(outfile)
    quit(save = "no", status = 0)
  }

  manifest[, project := project]

  # ------------------------------------------------------------------
  # 2. Clinical data (patient-level — joined to every sample via patient_id)
  # ------------------------------------------------------------------
  safe_col <- function(dt, candidates, default = NA_character_) {
    m <- intersect(candidates, colnames(dt))
    if (length(m)) dt[[m[1]]] else rep(default, nrow(dt))
  }

  clinical_raw <- tryCatch(
    GDCquery_clinic(project = project, type = "clinical"),
    error = function(e) { message("Clinical query failed (skipping): ", conditionMessage(e)); NULL }
  )

  if (!is.null(clinical_raw) && nrow(clinical_raw) > 0) {
    clin <- as.data.table(clinical_raw)
    clin_ann <- data.table(
      patient_id             = safe_col(clin, c("submitter_id", "bcr_patient_barcode")),
      gender                 = safe_col(clin, c("gender", "sex")),
      age_at_diagnosis       = suppressWarnings(as.numeric(safe_col(clin, c("age_at_index", "age_at_diagnosis")))),
      tumor_stage            = safe_col(clin, c("ajcc_pathologic_stage", "tumor_stage", "clinical_stage")),
      vital_status           = safe_col(clin, c("vital_status")),
      days_to_death          = suppressWarnings(as.numeric(safe_col(clin, c("days_to_death")))),
      days_to_last_follow_up = suppressWarnings(as.numeric(safe_col(clin, c("days_to_last_follow_up"))))
    )
    ann <- merge(manifest, clin_ann, by = "patient_id", all.x = TRUE)
  } else {
    ann <- copy(manifest)
    for (col in c("gender", "tumor_stage", "vital_status")) ann[, (col) := NA_character_]
    for (col in c("age_at_diagnosis", "days_to_death", "days_to_last_follow_up")) ann[, (col) := NA_real_]
  }

  # ------------------------------------------------------------------
  # 3. Molecular subtype (not available for all projects; NA if missing)
  # ------------------------------------------------------------------
  ann[, subtype := NA_character_]
  tryCatch({
    sub_df <- TCGAquery_subtype(tumor = sub("TCGA-", "", project))
    if (!is.null(sub_df) && nrow(sub_df) > 0) {
      sub_dt  <- as.data.table(sub_df)
      sub_col <- intersect(c("Subtype_Selected", "Subtype_mRNA", "Subtype", "subtype"), colnames(sub_dt))
      pid_col <- intersect(c("patient", "bcr_patient_barcode", "submitter_id"), colnames(sub_dt))
      if (length(sub_col) && length(pid_col)) {
        sub_map        <- as.character(sub_dt[[sub_col[1]]])
        names(sub_map) <- sub_dt[[pid_col[1]]]
        ann[, subtype := unname(sub_map[patient_id])]
      }
    }
  }, error = function(e) {
    message("Subtype query failed for ", project, " (skipping): ", conditionMessage(e))
  })

  # ------------------------------------------------------------------
  # 4. Tumor purity — 5 independent algorithms + consensus estimate
  # TCGAtumor_purity() matches on the first 16 chars of the barcode.
  # Algorithms:
  #   CPE      = consensus (median of available methods) → renamed tumor_purity
  #   ABSOLUTE = based on somatic copy-number + LOH
  #   ESTIMATE = based on gene-expression (high immune/stromal → lower purity)
  #   LUMP     = based on DNA methylation at leukocyte-specific loci
  #   IHC      = pathologist estimate from immunohistochemistry
  #   TIMER    = based on immune-cell deconvolution
  # ------------------------------------------------------------------
  purity_dt <- tryCatch({
    pur <- TCGAtumor_purity(samples = unique(ann$sample_id))
    if (is.null(pur) || nrow(pur) == 0) return(NULL)
    setDT(pur)
    id_col <- intersect(c("TCGA_sample", "sample", "barcode"), colnames(pur))[1]
    if (is.na(id_col)) return(NULL)
    setnames(pur, id_col, "sample_id")
    pur
  }, error = function(e) {
    message("Tumor purity query failed for ", project, " (skipping): ", conditionMessage(e))
    NULL
  })

  if (!is.null(purity_dt)) {
    rename_if <- function(dt, old, new) {
      if (old %in% colnames(dt) && !(new %in% colnames(dt))) setnames(dt, old, new)
    }
    rename_if(purity_dt, "CPE",      "tumor_purity")
    rename_if(purity_dt, "ABSOLUTE", "ABSOLUTE_purity")
    rename_if(purity_dt, "ESTIMATE", "ESTIMATE_purity")
    rename_if(purity_dt, "LUMP",     "LUMP_purity")
    rename_if(purity_dt, "IHC",      "IHC_purity")
    rename_if(purity_dt, "TIMER",    "TIMER_purity")

    keep_cols <- c("sample_id",
                   intersect(c("tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                               "LUMP_purity",  "IHC_purity",      "TIMER_purity"),
                             colnames(purity_dt)))
    ann <- merge(ann, purity_dt[, ..keep_cols], by = "sample_id", all.x = TRUE)
  } else {
    for (col in c("tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                  "LUMP_purity", "IHC_purity", "TIMER_purity"))
      ann[, (col) := NA_real_]
  }

  # ------------------------------------------------------------------
  # 5. Final column order
  # ------------------------------------------------------------------
  preferred <- c("barcode", "patient_id", "sample_id", "project",
                 "sample_type", "sample_type_code", "subtype",
                 "tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                 "LUMP_purity",  "IHC_purity",      "TIMER_purity",
                 "gender", "age_at_diagnosis", "tumor_stage",
                 "vital_status", "days_to_death", "days_to_last_follow_up")
  have  <- intersect(preferred, colnames(ann))
  extra <- setdiff(colnames(ann), have)
  setcolorder(ann, c(have, extra))

  fwrite(ann, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(ann), " sample records for ", project, " to ", outfile)

}, error = function(e) {
  message("ERROR in build_annotation for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
