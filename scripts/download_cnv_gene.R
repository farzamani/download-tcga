suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache  gene_annotation_file
# sample_type: full tissue.definition text such as "Primary Tumor"
#              (NOT the shortLetterCode "TP" — see TCGAbiolinks::getBarcodeDefinition()),
#              or "all" to download every sample type
# gdc_cache:   directory where GDCdownload stores raw files (e.g. results/cache)
# gene_annotation_file: annotation/gene_annotation.tsv (from annotate_genes.R)
#   — the mapping step: restricts/aligns CNV genes to the same Ensembl gene
#   IDs used in mrna.tsv, so both modalities describe the same gene set.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: download_cnv_gene.R <project> <sample_type> <outfile> <gdc_cache> <gene_annotation_file>")
}

project      <- args[1]
sample_type  <- args[2]
abs_path     <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile      <- abs_path(args[3])
gdc_cache    <- abs_path(args[4])
gene_ann_file <- abs_path(args[5])

dir.create(gdc_cache, recursive = TRUE, showWarnings = FALSE)
setwd(gdc_cache)

write_empty <- function(path) {
  fwrite(data.frame(barcode = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty cnv_gene.tsv for ", project)
}

tryCatch({
  # Genes to align to: the same protein-coding Ensembl gene set used for
  # mrna.tsv. GDC's gene-level copy number covers a broader gene model
  # (not protein-coding only), so this is the actual "mapping" step —
  # a plain ID intersection, no genomic-coordinate overlap needed since
  # GDC already gene-summarized the segments on their end.
  gene_ann <- fread(gene_ann_file, showProgress = FALSE)
  if (!"gene_id" %in% colnames(gene_ann)) {
    stop("gene_annotation_file has no gene_id column: ", gene_ann_file)
  }
  if (nrow(gene_ann) == 0) {
    stop("gene_annotation_file is empty (0 genes): ", gene_ann_file,
         " - annotate_genes.R likely failed; ",
         "check logs/annotate_genes.log and rerun that rule before retrying this one")
  }
  target_genes <- unique(gene_ann$gene_id)

  query_args <- list(
    project       = project,
    data.category = "Copy Number Variation",
    data.type     = "Gene Level Copy Number"
  )
  if (sample_type != "all") query_args$sample.type <- sample_type

  query <- do.call(GDCquery, query_args)
  GDCdownload(query, method = "api", files.per.chunk = 100, directory = gdc_cache)

  # summarizedExperiment=FALSE skips TCGAbiolinks' colData/clinical merge
  # (make_se_from_gene_level_copy_number -> colDataPrepare), which hits the
  # same per-project clinical-column bug as mRNA/methylation (see
  # download_mrna.R). We only need the copy_number values, so build the
  # matrix ourselves from the wide data.frame GDCprepare returns instead.
  # GDCprepare returns a tibble here (not a data.table like the mRNA path),
  # so force a plain data.frame to guarantee base `[` column-selection
  # semantics regardless of what class future TCGAbiolinks versions return.
  df <- as.data.frame(GDCprepare(query, directory = gdc_cache, summarizedExperiment = FALSE))
  if (!"gene_id" %in% colnames(df)) {
    stop("GDCprepare(summarizedExperiment=FALSE) returned an unexpected shape ",
         "(no gene_id column) - TCGAbiolinks internals may have changed")
  }

  # Per-sample columns are named "<barcode>_copy_number" /
  # "<barcode>_min_copy_number" / "<barcode>_max_copy_number"; keep only the
  # plain (non-min/max) copy_number value per sample.
  cn_cols <- grep("_copy_number$", colnames(df), value = TRUE)
  cn_cols <- cn_cols[!grepl("_(min|max)_copy_number$", cn_cols)]
  if (length(cn_cols) == 0) {
    stop("no <barcode>_copy_number columns found - TCGAbiolinks internals may have changed")
  }

  gene_ids <- sub("\\.[0-9]+$", "", df$gene_id)  # strip Ensembl version suffix
  keep     <- gene_ids %in% target_genes
  if (!any(keep)) {
    stop("none of this project's CNV gene_ids matched gene_annotation_file — ",
         "check gene ID formats agree")
  }

  cn_mat <- as.matrix(df[keep, cn_cols, drop = FALSE])
  mode(cn_mat) <- "numeric"
  rownames(cn_mat) <- gene_ids[keep]
  colnames(cn_mat) <- sub("_copy_number$", "", cn_cols)

  # A gene can appear more than once in GDC's file (rare, alt scaffolds);
  # keep the first occurrence to guarantee one column per gene downstream.
  cn_mat <- cn_mat[!duplicated(rownames(cn_mat)), , drop = FALSE]

  barcodes <- colnames(cn_mat)
  out <- as.data.frame(t(cn_mat))
  out <- cbind(data.frame(barcode = barcodes, stringsAsFactors = FALSE), out)

  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(out), " samples x ", ncol(out) - 1L,
          " genes (of ", length(target_genes), " target genes) to ", outfile)

}, error = function(e) {
  message("ERROR in download_cnv_gene for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
