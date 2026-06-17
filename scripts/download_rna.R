suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile
# sample_type: TCGA code such as "TP", or "all" to download every sample type
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: download_rna.R <project> <sample_type> <outfile>")

project     <- args[1]
sample_type <- args[2]
outfile     <- args[3]

barcode_to_patient <- function(bc) substr(bc, 1, 12)
barcode_to_sample  <- function(bc) substr(bc, 1, 16)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character(), patient_id = character(),
                    sample_id = character(), project = character(),
                    sample_type = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty rna.tsv for ", project)
}

tryCatch({
  query_args <- list(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  if (sample_type != "all") query_args$sample.type <- sample_type

  query <- do.call(GDCquery, query_args)
  GDCdownload(query, method = "api", files.per.chunk = 100)
  se <- GDCprepare(query)

  counts_mat <- assay(se, "unstranded")   # genes × samples
  barcodes   <- colnames(counts_mat)

  out <- as.data.frame(t(counts_mat))
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
  message("ERROR in download_rna for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
