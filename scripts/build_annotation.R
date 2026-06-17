suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile
# sample_type is informational here; clinical data is patient-level (no filter)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: build_annotation.R <project> <sample_type> <outfile>")

project     <- args[1]
sample_type <- args[2]   # recorded in output but clinical query is always all patients
outfile     <- args[3]

write_empty <- function(path) {
  fwrite(data.frame(
    barcode = character(), patient_id = character(), sample_id = character(),
    project = character(), sample_type = character(), subtype = character(),
    tumor_purity = numeric(), gender = character(),
    age_at_diagnosis = numeric(), tumor_stage = character(),
    vital_status = character(), days_to_death = numeric(),
    days_to_last_follow_up = numeric()
  ), path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty annotation.tsv for ", project)
}

tryCatch({
  # ------------------------------------------------------------------
  # 1. Clinical data (always patient-level, no sample_type filter)
  # ------------------------------------------------------------------
  clinical <- tryCatch(
    GDCquery_clinic(project = project, type = "clinical"),
    error = function(e) {
      message("Clinical query failed for ", project, ": ", conditionMessage(e)); NULL
    }
  )

  if (is.null(clinical) || nrow(clinical) == 0) {
    message("No clinical data for ", project, "; writing empty annotation")
    write_empty(outfile)
    quit(save = "no", status = 0)
  }

  clinical <- as.data.table(clinical)

  safe_col <- function(dt, candidates, default = NA_character_) {
    m <- intersect(candidates, colnames(dt))
    if (length(m)) dt[[m[1]]] else rep(default, nrow(dt))
  }

  ann <- data.table(
    patient_id             = safe_col(clinical, c("submitter_id", "bcr_patient_barcode")),
    gender                 = safe_col(clinical, c("gender", "sex")),
    age_at_diagnosis       = suppressWarnings(as.numeric(safe_col(clinical, c("age_at_index", "age_at_diagnosis")))),
    tumor_stage            = safe_col(clinical, c("ajcc_pathologic_stage", "tumor_stage", "clinical_stage")),
    vital_status           = safe_col(clinical, c("vital_status")),
    days_to_death          = suppressWarnings(as.numeric(safe_col(clinical, c("days_to_death")))),
    days_to_last_follow_up = suppressWarnings(as.numeric(safe_col(clinical, c("days_to_last_follow_up"))))
  )

  # ------------------------------------------------------------------
  # 2. Subtype (not available for all projects)
  # ------------------------------------------------------------------
  subtype_col <- rep(NA_character_, nrow(ann))
  tryCatch({
    sub_df <- TCGAquery_subtype(tumor = sub("TCGA-", "", project))
    if (!is.null(sub_df) && nrow(sub_df) > 0) {
      sub_df  <- as.data.table(sub_df)
      sub_col <- intersect(c("Subtype_Selected", "Subtype", "subtype"), colnames(sub_df))
      pid_col <- intersect(c("patient", "bcr_patient_barcode", "submitter_id"), colnames(sub_df))
      if (length(sub_col) && length(pid_col)) {
        sub_map           <- sub_df[[sub_col[1]]]
        names(sub_map)    <- sub_df[[pid_col[1]]]
        subtype_col       <- unname(sub_map[ann$patient_id])
      }
    }
  }, error = function(e) {
    message("Subtype query failed for ", project, " (skipping): ", conditionMessage(e))
  })
  ann[, subtype := subtype_col]

  # ------------------------------------------------------------------
  # 3. Tumor purity
  # ------------------------------------------------------------------
  purity_col <- rep(NA_real_, nrow(ann))
  tryCatch({
    pur_df <- TCGAtumor_purity(barcodes = ann$patient_id, db = "ESTIMATE")
    if (!is.null(pur_df) && nrow(pur_df) > 0) {
      pur_dt  <- as.data.table(pur_df)
      pid_col <- intersect(c("patient_id", "submitter_id", "sample"), colnames(pur_dt))
      pur_col <- intersect(c("purity", "TumorPurity", "Purity"), colnames(pur_dt))
      if (length(pid_col) && length(pur_col)) {
        pur_map        <- as.numeric(pur_dt[[pur_col[1]]])
        names(pur_map) <- pur_dt[[pid_col[1]]]
        purity_col     <- unname(pur_map[ann$patient_id])
      }
    }
  }, error = function(e) {
    message("Tumor purity query failed for ", project, " (skipping): ", conditionMessage(e))
  })
  ann[, tumor_purity := purity_col]

  # ------------------------------------------------------------------
  # 4. Finalise — clinical data is patient-level; record config sample_type
  # ------------------------------------------------------------------
  ann[, barcode     := patient_id]
  ann[, sample_id   := patient_id]
  ann[, project     := project]
  ann[, sample_type := sample_type]

  setcolorder(ann, c("barcode", "patient_id", "sample_id", "project", "sample_type",
                     "subtype", "tumor_purity", "gender", "age_at_diagnosis",
                     "tumor_stage", "vital_status", "days_to_death",
                     "days_to_last_follow_up"))

  fwrite(ann, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(ann), " patient records to ", outfile)

}, error = function(e) {
  message("ERROR in build_annotation for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
