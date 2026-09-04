suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
  library(AnnotationFilter)
  library(GenomicRanges)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: outfile
# Produces a gene annotation table for human protein-coding genes only,
# keyed on Ensembl gene ID (the row IDs in rna.tsv).
#
# Uses the offline Bioconductor package EnsDb.Hsapiens.v86 (Ensembl release
# 86, GRCh38) rather than a live biomaRt/Ensembl query — biomaRt needs
# outbound HTTPS to ensembl.org and its mirrors, which many HPC compute
# nodes block/proxy (seen here as HTTP 500/502/403 from all three mirrors).
# EnsDb.Hsapiens.v86 is installed via conda and needs no network access.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: annotate_genes.R <outfile>")
outfile <- args[1]

STANDARD_CHROMS <- c(as.character(1:22), "X", "Y", "MT")

write_empty <- function(path) {
  fwrite(data.frame(gene_id = character(), gene_name = character(),
                    gene_type = character(), chromosome = character(),
                    start = integer(), end = integer(), strand = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty gene_annotation.tsv")
}

tryCatch({
  edb <- EnsDb.Hsapiens.v86
  g   <- genes(edb, filter = GeneBiotypeFilter("protein_coding"))

  ann <- data.table(
    gene_id    = sub("\\.[0-9]+$", "", g$gene_id),  # strip version, if any
    gene_name  = g$gene_name,
    gene_type  = g$gene_biotype,
    chromosome = as.character(seqnames(g)),
    start      = start(g),
    end        = end(g),
    strand     = as.character(strand(g))
  )

  # Restrict to standard chromosomes
  ann <- ann[chromosome %in% STANDARD_CHROMS]
  setorder(ann, chromosome, start)

  fwrite(ann, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(ann), " gene records to ", outfile)

}, error = function(e) {
  message("ERROR in annotate_genes: ", conditionMessage(e))
  write_empty(outfile)
})
