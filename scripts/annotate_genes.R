suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: outfile
# Produces a gene annotation table for human protein-coding genes only,
# keyed on Ensembl gene ID (the row IDs in rna.tsv).
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
  # ensembl.org occasionally redirects/rate-limits or is unreachable from
  # some networks (e.g. HPC compute nodes with restricted outbound access);
  # try a couple of Ensembl's other mirror hosts before giving up.
  hosts <- c("https://ensembl.org", "https://useast.ensembl.org", "https://asia.ensembl.org")
  mart  <- NULL
  last_error <- NULL
  for (h in hosts) {
    mart <- tryCatch(
      useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = h),
      error = function(e) { last_error <<- e; NULL }
    )
    if (!is.null(mart)) {
      message("Connected to Ensembl BioMart via ", h)
      break
    }
    message("Could not reach Ensembl BioMart via ", h, " — ", conditionMessage(last_error))
  }
  if (is.null(mart)) stop("Could not reach any Ensembl BioMart mirror: ", conditionMessage(last_error))

  ann <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                   "chromosome_name", "start_position", "end_position", "strand"),
    filters    = "biotype",
    values     = "protein_coding",
    mart       = mart
  )
  setDT(ann)
  setnames(ann, c("gene_id", "gene_name", "gene_type",
                  "chromosome", "start", "end", "strand"))

  # Restrict to standard chromosomes and convert strand to +/-
  ann <- ann[chromosome %in% STANDARD_CHROMS]
  ann[, strand := ifelse(strand == 1, "+", "-")]
  setorder(ann, chromosome, start)

  fwrite(ann, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(ann), " gene records to ", outfile)

}, error = function(e) {
  message("ERROR in annotate_genes: ", conditionMessage(e))
  write_empty(outfile)
})
