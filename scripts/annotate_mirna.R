suppressPackageStartupMessages({
  library(miRBaseConverter)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: outfile
# Produces a miRNA annotation table keyed on precursor miRNA ID (e.g. hsa-mir-21),
# which matches the column names used in mirna.tsv.
# Source: miRBaseConverter v22 (bundled database, no internet required).
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: annotate_mirna.R <outfile>")
outfile <- args[1]

write_empty <- function(path) {
  fwrite(data.frame(mirna_id = character(), accession = character(),
                    mature_id_5p = character(), mature_acc_5p = character(),
                    mature_id_3p = character(), mature_acc_3p = character()),
         path, sep = "\t", quote = FALSE)
  message("WARNING: wrote empty mirna_annotation.tsv")
}

tryCatch({
  # getAllMiRNAs returns columns: Accession, Name, Sequence (no Type column).
  # Distinguish precursor vs mature by Accession prefix:
  #   MI0...   = precursor
  #   MIMAT... = mature
  all_mirnas <- getAllMiRNAs(version = "v22", type = "all", species = "hsa")
  setDT(all_mirnas)

  pre <- all_mirnas[startsWith(Accession, "MI0"),
                    .(mirna_id = Name, accession = Accession)]

  mat <- all_mirnas[startsWith(Accession, "MIMAT"),
                    .(mature_id = Name, mature_acc = Accession)]

  # Derive expected precursor name from mature name:
  #   hsa-miR-21-5p   → hsa-mir-21
  #   hsa-let-7a-2-3p → hsa-let-7a-2
  # Step 1: strip -5p / -3p / * suffix
  # Step 2: mature IDs use capital miR; precursor names use lowercase mir
  mat[, pre_name := gsub("-[35]p$|\\*$", "", mature_id)]
  mat[, pre_name := gsub("-miR-", "-mir-", pre_name)]
  mat[, pre_name := gsub("-miR$",  "-mir",  pre_name)]

  # Keep one 5p and one 3p per precursor (first match wins)
  mat5 <- mat[grepl("-5p$", mature_id)][!duplicated(pre_name),
              .(pre_name, mature_id_5p = mature_id, mature_acc_5p = mature_acc)]
  mat3 <- mat[grepl("-3p$", mature_id)][!duplicated(pre_name),
              .(pre_name, mature_id_3p = mature_id, mature_acc_3p = mature_acc)]

  ann <- merge(pre, mat5, by.x = "mirna_id", by.y = "pre_name", all.x = TRUE)
  ann <- merge(ann, mat3, by.x = "mirna_id", by.y = "pre_name", all.x = TRUE)

  setcolorder(ann, c("mirna_id", "accession",
                     "mature_id_5p", "mature_acc_5p",
                     "mature_id_3p", "mature_acc_3p"))
  setorder(ann, mirna_id)

  fwrite(ann, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(ann), " precursor miRNA records to ", outfile)

}, error = function(e) {
  message("ERROR in annotate_mirna: ", conditionMessage(e))
  write_empty(outfile)
})
