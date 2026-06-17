suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(matrixStats)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  [max_cpgs]
# sample_type: TCGA code such as "TP", or "all" to download every sample type
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: download_methylation.R <project> <sample_type> <outfile> [max_cpgs]")

project     <- args[1]
sample_type <- args[2]
outfile     <- args[3]
max_cpgs    <- if (length(args) >= 4 && args[4] != "null") as.integer(args[4]) else NULL

barcode_to_patient <- function(bc) substr(bc, 1, 12)
barcode_to_sample  <- function(bc) substr(bc, 1, 16)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character(), patient_id = character(),
                    sample_id = character(), project = character(),
                    sample_type = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty methylation.tsv for ", project)
}

tryCatch({
  # Try 450k first; fall back to 27k if not available for this project
  query <- tryCatch({
    q_args <- list(
      project       = project,
      data.category = "DNA Methylation",
      data.type     = "Methylation Beta Value",
      platform      = "Illumina Human Methylation 450"
    )
    if (sample_type != "all") q_args$sample.type <- sample_type
    do.call(GDCquery, q_args)
  }, error = function(e) {
    message("450k not found, trying 27k for ", project)
    q_args <- list(
      project       = project,
      data.category = "DNA Methylation",
      data.type     = "Methylation Beta Value",
      platform      = "Illumina Human Methylation 27"
    )
    if (sample_type != "all") q_args$sample.type <- sample_type
    do.call(GDCquery, q_args)
  })

  GDCdownload(query, method = "api", files.per.chunk = 50)
  se <- GDCprepare(query)

  beta_mat <- assay(se)   # CpGs × samples

  if (!is.null(max_cpgs) && nrow(beta_mat) > max_cpgs) {
    row_vars <- rowVars(beta_mat, na.rm = TRUE)
    keep     <- order(row_vars, decreasing = TRUE)[seq_len(max_cpgs)]
    beta_mat <- beta_mat[keep, , drop = FALSE]
    message("Kept top ", max_cpgs, " variable CpGs out of ", length(row_vars))
  }

  barcodes <- colnames(beta_mat)

  out <- as.data.frame(t(beta_mat))
  out <- cbind(
    data.frame(barcode     = barcodes,
               patient_id  = barcode_to_patient(barcodes),
               sample_id   = barcode_to_sample(barcodes),
               project     = project,
               sample_type = substr(barcodes, 14, 15),
               stringsAsFactors = FALSE),
    out
  )

  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(out), " samples × ", ncol(out) - 5, " CpGs to ", outfile)

}, error = function(e) {
  message("ERROR in download_methylation for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
