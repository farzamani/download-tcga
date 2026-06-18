suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(matrixStats)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache  [max_cpgs]
# sample_type: TCGA code such as "TP", or "all" to download every sample type
# gdc_cache:   directory where GDCdownload stores raw files (e.g. results/cache)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_methylation.R <project> <sample_type> <outfile> <gdc_cache> [max_cpgs]")

project     <- args[1]
sample_type <- args[2]
abs_path  <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile   <- abs_path(args[3])
gdc_cache <- abs_path(args[4])
max_cpgs  <- if (length(args) >= 5 && args[5] != "null") as.integer(args[5]) else NULL

dir.create(gdc_cache, recursive = TRUE, showWarnings = FALSE)
setwd(gdc_cache)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character()),
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

  # Methylation chunks are large (~130 MB each at 10 files/chunk).
  # Retry up to 3 times on transient GDC network errors.
  max_attempts <- 3L
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch({
      GDCdownload(query, method = "api", files.per.chunk = 10, directory = gdc_cache)
      "ok"
    }, error = function(e) conditionMessage(e))
    if (identical(result, "ok")) break
    if (attempt < max_attempts) {
      message("Download attempt ", attempt, " failed — retrying: ", result)
      Sys.sleep(30)
    } else {
      stop("GDCdownload failed after ", max_attempts, " attempts: ", result)
    }
  }

  # summarizedExperiment=FALSE returns a plain matrix (CpGs × samples) and
  # avoids the sesameData dependency introduced in TCGAbiolinks >= 2.29
  beta_mat <- GDCprepare(query, directory = gdc_cache, summarizedExperiment = FALSE)

  if (!is.null(max_cpgs) && nrow(beta_mat) > max_cpgs) {
    row_vars <- rowVars(beta_mat, na.rm = TRUE)
    keep     <- order(row_vars, decreasing = TRUE)[seq_len(max_cpgs)]
    beta_mat <- beta_mat[keep, , drop = FALSE]
    message("Kept top ", max_cpgs, " variable CpGs out of ", length(row_vars))
  }

  barcodes <- colnames(beta_mat)

  out <- as.data.frame(t(beta_mat))
  out <- cbind(data.frame(barcode = barcodes, stringsAsFactors = FALSE), out)

  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(out), " samples × ", ncol(out) - 1L, " CpGs to ", outfile)

}, error = function(e) {
  message("ERROR in download_methylation for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
