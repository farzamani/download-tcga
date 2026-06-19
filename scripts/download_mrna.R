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
  se <- GDCprepare(query, directory = gdc_cache)

  counts_mat <- assay(se, "unstranded")   # genes × samples
  barcodes   <- colnames(counts_mat)

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
