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
abs_path <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile  <- abs_path(args[2])

write_empty <- function(path) {
  fwrite(data.frame(
    barcode = character(), patient_id = character(), sample_id = character(),
    project = character(), sample_type = character(), sample_type_code = character(),
    subtype = character(),
    tumor_purity    = numeric(), ABSOLUTE_purity = numeric(),
    ESTIMATE_purity = numeric(), LUMP_purity     = numeric(),
    IHC_purity      = numeric(),
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
  # 3. Molecular subtype — Pan-Cancer Atlas subtypes (all projects)
  # PanCancerAtlas_subtypes() returns a bundled tibble keyed on 12-char
  # patient barcode (pan.samplesID).  Subtype_Selected is the recommended
  # consensus subtype for each cancer type.
  # ------------------------------------------------------------------
  ann[, subtype := NA_character_]
  tryCatch({
    all_sub     <- as.data.table(PanCancerAtlas_subtypes())
    cancer_type <- sub("TCGA-", "", project)
    proj_sub    <- all_sub[cancer.type == cancer_type]
    if (nrow(proj_sub) > 0) {
      # pan.samplesID holds the full aliquot barcode (e.g. TCGA-3C-AAAU-01A-11R-A41B-07);
      # truncate to 12 chars to get the patient-level barcode for joining.
      proj_sub[, pid := substr(pan.samplesID, 1, 12)]
      sub_map        <- as.character(proj_sub$Subtype_Selected)
      names(sub_map) <- proj_sub$pid
      ann[, subtype := unname(sub_map[patient_id])]
    }
  }, error = function(e) {
    message("Subtype query failed for ", project, " (skipping): ", conditionMessage(e))
  })

  # ------------------------------------------------------------------
  # 4. Tumor purity — 4 algorithms + consensus (from bundled TCGAbiolinks dataset)
  # Tumor.purity is keyed on Sample.ID (16-char barcode, e.g. TCGA-A1-A0SO-01A).
  # Values may be stored with comma-decimal notation depending on R locale;
  # convert to numeric defensively.
  # Columns:
  #   CPE      = consensus purity estimate (median of available methods)
  #   ESTIMATE = expression-based (high immune/stromal → lower purity)
  #   ABSOLUTE = copy-number / LOH based
  #   LUMP     = DNA methylation at leukocyte-specific loci
  #   IHC      = pathologist estimate from immunohistochemistry
  # ------------------------------------------------------------------
  purity_dt <- tryCatch({
    data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
    pur <- as.data.table(Tumor.purity)
    to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
    for (col in intersect(c("CPE", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC"), colnames(pur)))
      pur[, (col) := to_num(get(col))]
    id_col <- intersect(c("Sample.ID", "sample_id"), colnames(pur))[1]
    if (is.na(id_col)) return(NULL)
    setnames(pur, id_col, "sample_id")
    setnames(pur, c("CPE", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC"),
             c("tumor_purity", "ESTIMATE_purity", "ABSOLUTE_purity", "LUMP_purity", "IHC_purity"),
             skip_absent = TRUE)
    pur[sample_id %in% ann$sample_id]
  }, error = function(e) {
    message("Tumor purity load failed (skipping): ", conditionMessage(e))
    NULL
  })

  if (!is.null(purity_dt) && nrow(purity_dt) > 0) {
    keep_cols <- c("sample_id",
                   intersect(c("tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                               "LUMP_purity",  "IHC_purity"),
                             colnames(purity_dt)))
    ann <- merge(ann, purity_dt[, ..keep_cols], by = "sample_id", all.x = TRUE)
  } else {
    for (col in c("tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                  "LUMP_purity", "IHC_purity"))
      ann[, (col) := NA_real_]
  }

  # ------------------------------------------------------------------
  # 5. Final column order
  # ------------------------------------------------------------------
  preferred <- c("barcode", "patient_id", "sample_id", "project",
                 "sample_type", "sample_type_code", "subtype",
                 "tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity",
                 "LUMP_purity",  "IHC_purity",
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
