suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache
# sample_type: full tissue.definition text such as "Primary Tumor"
#              (NOT the shortLetterCode "TP" — see TCGAbiolinks::getBarcodeDefinition()),
#              or "all" to download every sample type
# gdc_cache:   directory where GDCdownload stores raw files (e.g. results/cache)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_mrna.R <project> <sample_type> <outfile> <gdc_cache>")

project     <- args[1]
sample_type <- args[2]
# Resolve to absolute paths before setwd() changes the working directory.
# normalizePath(mustWork=FALSE) is unreliable for non-existent paths on macOS;
# file.path(getwd(), ...) always works.
abs_path  <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile   <- abs_path(args[3])
gdc_cache <- abs_path(args[4])

# Redirect all TCGAbiolinks temp files (tar.gz chunks, MANIFEST.txt, query
# cache) to gdc_cache by making it the working directory for this script.
dir.create(gdc_cache, recursive = TRUE, showWarnings = FALSE)
setwd(gdc_cache)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty mrna.tsv for ", project)
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
  GDCdownload(query, method = "api", files.per.chunk = 100, directory = gdc_cache)

  # summarizedExperiment=FALSE skips TCGAbiolinks' colData/clinical merge
  # (makeSEfromTranscriptomeProfilingSTAR -> colDataPrepare), which errors
  # on projects like TCGA-LAML whose indexed clinical data is missing
  # columns (e.g. disease_response) that TCGAbiolinks assumes exist. We
  # only need raw counts, so build the matrix ourselves from the wide
  # data.table GDCprepare returns instead.
  df <- GDCprepare(query, directory = gdc_cache, summarizedExperiment = FALSE)
  if (!"gene_id" %in% colnames(df)) {
    stop("GDCprepare(summarizedExperiment=FALSE) returned an unexpected shape ",
         "(no gene_id column) - TCGAbiolinks internals may have changed")
  }
  df <- df[grepl("^ENSG", df$gene_id), ]
  if (nrow(df) == 0) stop("no ENSG gene rows found after filtering GDCprepare output")

  count_cols <- grep("^unstranded_", colnames(df), value = TRUE)
  if (length(count_cols) == 0) {
    stop("no unstranded_* count columns found - TCGAbiolinks internals may have changed")
  }
  counts_mat <- as.matrix(df[, ..count_cols])
  rownames(counts_mat) <- df$gene_id
  colnames(counts_mat) <- sub("^unstranded_", "", count_cols)
  barcodes <- colnames(counts_mat)

  # Strip Ensembl version suffix so IDs match gene_annotation.tsv
  # ENSG00000000003.15 → ENSG00000000003
  rownames(counts_mat) <- sub("\\.[0-9]+$", "", rownames(counts_mat))

  out <- as.data.frame(t(counts_mat))
  out <- cbind(data.frame(barcode = barcodes, stringsAsFactors = FALSE), out)

  fwrite(out, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(out), " samples x ",
          ncol(out) - 1L, " genes to ", outfile)

}, error = function(e) {
  message("ERROR in download_mrna for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
