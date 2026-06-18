suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_mirna.R <project> <sample_type> <outfile> <gdc_cache>")

project     <- args[1]
sample_type <- args[2]
outfile     <- normalizePath(args[3], mustWork = FALSE)
gdc_cache   <- normalizePath(args[4], mustWork = FALSE)

dir.create(gdc_cache, recursive = TRUE, showWarnings = FALSE)
setwd(gdc_cache)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty mirna.tsv for ", project)
}

tryCatch({
  query_args <- list(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "miRNA Expression Quantification",
    workflow.type = "BCGSC miRNA Profiling"
  )
  if (sample_type != "all") query_args$sample.type <- sample_type

  query <- do.call(GDCquery, query_args)
  GDCdownload(query, method = "api", files.per.chunk = 100, directory = gdc_cache)

  # GDCprepare returns a wide data.frame for miRNA (not SummarizedExperiment):
  # columns: miRNA_ID | read_count_<barcode> | reads_per_million_..._<barcode> | cross-mapped_<barcode>
  df <- GDCprepare(query, directory = gdc_cache)

  # Extract raw read count columns
  rc_cols  <- grep("^read_count_", colnames(df), value = TRUE)
  if (length(rc_cols) == 0) stop("No read_count_ columns found in GDCprepare output")

  barcodes <- sub("^read_count_", "", rc_cols)

  # Build miRNAs × samples count matrix, then transpose
  count_mat <- as.matrix(df[, rc_cols, drop = FALSE])
  rownames(count_mat) <- df$miRNA_ID
  colnames(count_mat) <- barcodes

  out <- as.data.frame(t(count_mat))
  out <- cbind(data.frame(barcode = barcodes, stringsAsFactors = FALSE), out)

  fwrite(out, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(out), " samples x ",
          length(rc_cols), " miRNAs to ", outfile)

}, error = function(e) {
  message("ERROR in download_mirna for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
