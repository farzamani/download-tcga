suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(jsonlite)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  outfile
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: build_annotation.R <project> <outfile>")

project  <- args[1]
abs_path <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile  <- abs_path(args[2])

write_empty <- function(path) {
  fwrite(data.frame(
    barcode = character(), patient_id = character(), sample_id = character(),
    project = character(), project_name = character(), primary_site = character(),
    sample_type = character(), sample_type_code = character(), subtype = character(),
    tumor_purity = numeric(), ABSOLUTE_purity = numeric(),
    ESTIMATE_purity = numeric(), LUMP_purity = numeric(), IHC_purity = numeric(),
    gender = character(), race = character(), ethnicity = character(),
    age_at_diagnosis = integer(), tumor_stage = character(),
    primary_diagnosis = character(), tissue_or_organ_of_origin = character(),
    vital_status = character(), days_to_death = numeric(),
    days_to_last_follow_up = numeric()
  ), path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty annotation.tsv for ", project)
}

tryCatch({

  # ------------------------------------------------------------------
  # 1. Sample manifest — barcodes + sample-type labels from GDC
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
  # 2. Project-level metadata: full name and primary anatomical site
  # ------------------------------------------------------------------
  tryCatch({
    all_proj <- as.data.table(getGDCprojects())
    row      <- all_proj[id == project]
    if (nrow(row) > 0) {
      manifest[, project_name := as.character(row$name[1])]
      manifest[, primary_site := paste(unlist(row$primary_site[[1]]), collapse = ", ")]
    } else {
      manifest[, c("project_name", "primary_site") := NA_character_]
    }
  }, error = function(e) {
    message("getGDCprojects failed (skipping): ", conditionMessage(e))
    manifest[, c("project_name", "primary_site") := NA_character_]
  })

  # ------------------------------------------------------------------
  # 3. Clinical data — direct GDC REST API query
  #
  # GDCquery_clinic() has a data.table version conflict that causes it to
  # crash whenever patients have multiple follow-up or diagnosis records
  # (e.g. TCGA-LGG, TCGA-BRCA).  Querying the GDC API directly avoids
  # the bug entirely and also gives us additional fields (race, ethnicity,
  # primary_diagnosis, tissue_or_organ_of_origin).
  #
  # age_at_diagnosis is returned in DAYS by the GDC API; converted to years.
  # ------------------------------------------------------------------
  gdc_clinical <- function(proj) {
    fields <- paste(c(
      "submitter_id",
      "demographic.gender", "demographic.race", "demographic.ethnicity",
      "demographic.vital_status", "demographic.days_to_death",
      "diagnoses.age_at_diagnosis", "diagnoses.ajcc_pathologic_stage",
      "diagnoses.primary_diagnosis", "diagnoses.tissue_or_organ_of_origin",
      "diagnoses.days_to_last_follow_up"
    ), collapse = ",")
    filt <- sprintf(
      '{"op":"=","content":{"field":"project.project_id","value":"%s"}}', proj
    )
    url <- paste0(
      "https://api.gdc.cancer.gov/cases?",
      "filters=", URLencode(filt, reserved = TRUE),
      "&fields=", fields,
      "&size=10000&format=JSON"
    )
    tryCatch(fromJSON(url)$data$hits, error = function(e) {
      message("GDC clinical API failed (skipping): ", conditionMessage(e))
      NULL
    })
  }

  # Extract the first non-NA value from a list-of-data.frames column
  first_val <- function(lst, col) {
    vapply(lst, function(x) {
      if (is.null(x) || !is.data.frame(x) || !col %in% names(x)) return(NA_character_)
      v <- na.omit(as.character(x[[col]]))
      if (length(v) == 0) NA_character_ else v[1]
    }, character(1))
  }

  hits <- gdc_clinical(project)
  if (!is.null(hits) && nrow(hits) > 0) {
    clin <- data.table(
      patient_id                = hits$submitter_id,
      gender                    = hits$demographic$gender,
      race                      = hits$demographic$race,
      ethnicity                 = hits$demographic$ethnicity,
      vital_status              = hits$demographic$vital_status,
      days_to_death             = suppressWarnings(as.numeric(hits$demographic$days_to_death)),
      age_at_diagnosis          = as.integer(
        suppressWarnings(as.numeric(first_val(hits$diagnoses, "age_at_diagnosis"))) / 365.25
      ),
      tumor_stage               = first_val(hits$diagnoses, "ajcc_pathologic_stage"),
      primary_diagnosis         = first_val(hits$diagnoses, "primary_diagnosis"),
      tissue_or_organ_of_origin = first_val(hits$diagnoses, "tissue_or_organ_of_origin"),
      days_to_last_follow_up    = suppressWarnings(
        as.numeric(first_val(hits$diagnoses, "days_to_last_follow_up"))
      )
    )
    ann <- merge(manifest, clin, by = "patient_id", all.x = TRUE)
  } else {
    ann <- copy(manifest)
    for (col in c("gender", "race", "ethnicity", "tumor_stage", "vital_status",
                  "primary_diagnosis", "tissue_or_organ_of_origin"))
      ann[, (col) := NA_character_]
    for (col in c("age_at_diagnosis", "days_to_death", "days_to_last_follow_up"))
      ann[, (col) := NA_integer_]
  }

  # ------------------------------------------------------------------
  # 4. Molecular subtype — Pan-Cancer Atlas
  # pan.samplesID holds the full aliquot barcode; truncate to 12 chars
  # to get the patient-level barcode for joining.
  # ------------------------------------------------------------------
  ann[, subtype := NA_character_]
  tryCatch({
    all_sub     <- as.data.table(PanCancerAtlas_subtypes())
    cancer_type <- sub("TCGA-", "", project)
    proj_sub    <- all_sub[cancer.type == cancer_type]
    if (nrow(proj_sub) > 0) {
      proj_sub[, pid := substr(pan.samplesID, 1, 12)]
      sub_map        <- as.character(proj_sub$Subtype_Selected)
      names(sub_map) <- proj_sub$pid
      ann[, subtype := unname(sub_map[patient_id])]
    }
  }, error = function(e) {
    message("Subtype query failed for ", project, " (skipping): ", conditionMessage(e))
  })

  # ------------------------------------------------------------------
  # 5. Tumor purity — bundled TCGAbiolinks dataset
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
  # 6. Final column order
  # ------------------------------------------------------------------
  preferred <- c(
    "barcode", "patient_id", "sample_id",
    "project", "project_name", "primary_site",
    "sample_type", "sample_type_code", "subtype",
    "tumor_purity", "ABSOLUTE_purity", "ESTIMATE_purity", "LUMP_purity", "IHC_purity",
    "gender", "race", "ethnicity",
    "age_at_diagnosis", "tumor_stage",
    "primary_diagnosis", "tissue_or_organ_of_origin",
    "vital_status", "days_to_death", "days_to_last_follow_up"
  )
  have  <- intersect(preferred, colnames(ann))
  extra <- setdiff(colnames(ann), have)
  setcolorder(ann, c(have, extra))

  fwrite(ann, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(ann), " sample records for ", project, " to ", outfile)

}, error = function(e) {
  message("ERROR in build_annotation for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
