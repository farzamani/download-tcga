suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache
# sample_type: TCGA code such as "TP", or "all" to download every sample type
# gdc_cache:   directory where GDCdownload stores raw files (e.g. results/cache)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_mirna.R <project> <sample_type> <outfile> <gdc_cache>")

project     <- args[1]
sample_type <- args[2]
outfile     <- args[3]
gdc_cache   <- args[4]

barcode_to_patient <- function(bc) substr(bc, 1, 12)
barcode_to_sample  <- function(bc) substr(bc, 1, 16)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character(), patient_id = character(),
                    sample_id = character(), project = character(),
                    sample_type = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty mirna.tsv for ", project)
}

tryCatch({
  # "miRNA Expression Quantification" returns miRNA-level data directly
  # (one row per miRNA, already aggregated from isoforms by GDC).
  query_args <- list(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "miRNA Expression Quantification"
  )
  if (sample_type != "all") query_args$sample.type <- sample_type

  query <- do.call(GDCquery, query_args)
  GDCdownload(query, method = "api", files.per.chunk = 100, directory = gdc_cache)
  se <- GDCprepare(query, directory = gdc_cache)

  # Prefer reads_per_million_miRNA_mapped; fall back to the first available assay
  assay_name <- if ("reads_per_million_miRNA_mapped" %in% assayNames(se)) {
    "reads_per_million_miRNA_mapped"
  } else {
    assayNames(se)[1]
  }
  expr_mat <- assay(se, assay_name)   # miRNAs × samples

  barcodes <- colnames(expr_mat)

  out <- as.data.frame(t(expr_mat))
  out <- cbind(
    data.frame(barcode     = barcodes,
               patient_id  = barcode_to_patient(barcodes),
               sample_id   = barcode_to_sample(barcodes),
               project     = project,
               sample_type = substr(barcodes, 14, 15),
               stringsAsFactors = FALSE),
    out
  )

  fwrite(out, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(out), " samples to ", outfile)

}, error = function(e) {
  message("ERROR in download_mirna for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
