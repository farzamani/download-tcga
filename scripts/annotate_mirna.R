suppressPackageStartupMessages({
  library(miRBaseConverter)
  library(biomaRt)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: outfile
# Produces a miRNA annotation table keyed on the miRNA_ID used in mirna.tsv
# (e.g. hsa-mir-21). Columns: mirna_id, accession, mature_id_5p, mature_id_3p,
# mature_acc_5p, mature_acc_3p, chromosome, start, end, strand.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: annotate_mirna.R <outfile>")
outfile <- args[1]

STANDARD_CHROMS <- c(as.character(1:22), "X", "Y", "MT")

write_empty <- function(path) {
  fwrite(data.frame(mirna_id = character(), accession = character(),
                    mature_id_5p = character(), mature_id_3p = character(),
                    mature_acc_5p = character(), mature_acc_3p = character(),
                    chromosome = character(), start = integer(),
                    end = integer(), strand = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty mirna_annotation.tsv")
}

tryCatch({
  # ------------------------------------------------------------------
  # 1. Get all human miRNA IDs and accessions from miRBaseConverter
  # ------------------------------------------------------------------
  all_ids <- getAllMiRNAs(version = "v22", type = "all", species = "hsa")
  setDT(all_ids)

  # Precursor miRNA table
  pre <- all_ids[Type == "Pre-miRNA", .(mirna_id = Accession,      # MI accession
                                        mirna_name = Name)]         # hsa-mir-21

  # Mature miRNA table
  mat <- all_ids[Type == "Mature miRNA", .(mature_acc = Accession,
                                           mature_id  = Name,
                                           pre_name   = gsub("-[35]p$", "", Name))]

  # Split into 5p and 3p
  mat_5p <- mat[grepl("-5p$", mature_id),
                .(pre_name, mature_id_5p = mature_id, mature_acc_5p = mature_acc)]
  mat_3p <- mat[grepl("-3p$", mature_id),
                .(pre_name, mature_id_3p = mature_id, mature_acc_3p = mature_acc)]

  # Join on precursor name
  ann <- pre[mat_5p, on = .(mirna_name = pre_name), nomatch = NA]
  ann <- ann[mat_3p, on = .(mirna_name = pre_name), nomatch = NA]
  setnames(ann, "mirna_id", "accession")
  setnames(ann, "mirna_name", "mirna_id")

  # ------------------------------------------------------------------
  # 2. Genomic coordinates via Ensembl biomaRt
  # ------------------------------------------------------------------
  mart <- tryCatch(
    useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://ensembl.org"),
    error = function(e) { message("biomaRt unavailable, skipping coordinates"); NULL }
  )

  if (!is.null(mart)) {
    coords <- getBM(
      attributes = c("mirbase_id", "chromosome_name",
                     "start_position", "end_position", "strand"),
      filters    = "biotype",
      values     = "miRNA",
      mart       = mart
    )
    setDT(coords)
    setnames(coords, c("mirna_id", "chromosome", "start", "end", "strand"))
    coords <- coords[chromosome %in% STANDARD_CHROMS]
    coords[, strand := ifelse(strand == 1, "+", "-")]
    # mirbase_id in Ensembl uses names like "MI0000077"; map via accession
    ann <- coords[ann, on = .(mirna_id = accession), nomatch = NA]
    setnames(ann, "mirna_id", "chromosome_mirna_id")   # Ensembl mirbase_id
    setnames(ann, "i.mirna_id", "mirna_id")
  } else {
    ann[, c("chromosome", "start", "end", "strand") := NA]
  }

  keep_cols <- intersect(c("mirna_id", "accession", "mature_id_5p", "mature_id_3p",
                            "mature_acc_5p", "mature_acc_3p",
                            "chromosome", "start", "end", "strand"),
                         colnames(ann))
  ann <- ann[, ..keep_cols]

  fwrite(ann, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(ann), " miRNA records to ", outfile)

}, error = function(e) {
  message("ERROR in annotate_mirna: ", conditionMessage(e))
  write_empty(outfile)
})
