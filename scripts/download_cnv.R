suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache
# sample_type: TCGA code such as "TP", or "all" to download every sample type
# gdc_cache:   directory where GDCdownload stores raw files (e.g. results/cache)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_cnv.R <project> <sample_type> <outfile> <gdc_cache>")

project     <- args[1]
sample_type <- args[2]
outfile     <- normalizePath(args[3], mustWork = FALSE)
gdc_cache   <- normalizePath(args[4], mustWork = FALSE)

dir.create(gdc_cache, recursive = TRUE, showWarnings = FALSE)
setwd(gdc_cache)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character(), chromosome = character(),
                    start = integer(), end = integer(),
                    num_probes = integer(), segment_mean = numeric()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty cnv.tsv for ", project)
}

tryCatch({
  query_args <- list(
    project       = project,
    data.category = "Copy Number Variation",
    data.type     = "Masked Copy Number Segment"
  )
  if (sample_type != "all") query_args$sample.type <- sample_type

  query <- do.call(GDCquery, query_args)
  GDCdownload(query, method = "api", files.per.chunk = 100, directory = gdc_cache)
  cnv_df <- GDCprepare(query, directory = gdc_cache)

  # Standardise column names (may vary between TCGAbiolinks versions)
  rename_col <- function(df, candidates) {
    m <- intersect(candidates, colnames(df))
    if (length(m)) df[[m[1]]] else rep(NA, nrow(df))
  }

  barcodes <- rename_col(cnv_df, c("Sample", "barcode", "Tumor_Sample_Barcode"))

  out <- data.table(
    barcode      = barcodes,
    chromosome   = rename_col(cnv_df, c("Chromosome", "chrom", "Chr")),
    start        = as.integer(rename_col(cnv_df, c("Start", "loc.start", "Start_Position"))),
    end          = as.integer(rename_col(cnv_df, c("End", "loc.end", "End_Position"))),
    num_probes   = as.integer(rename_col(cnv_df, c("Num_Probes", "num.mark"))),
    segment_mean = as.numeric(rename_col(cnv_df, c("Segment_Mean", "seg.mean")))
  )

  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(out), " segments from ",
          uniqueN(out$barcode), " samples to ", outfile)

}, error = function(e) {
  message("ERROR in download_cnv for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
