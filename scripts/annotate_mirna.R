suppressPackageStartupMessages({
  library(miRBaseConverter)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: outfile
# Produces a miRNA annotation table keyed on precursor miRNA ID (e.g. hsa-mir-21),
# which matches the column names used in mirna.tsv.
# Source: miRBaseConverter v22 (bundled database, no internet required).
#
# Uses miRNA_PrecursorToMature() for reliable precursor→mature mapping.
# Name-pattern stripping fails for multi-locus miRNAs (e.g. hsa-mir-9-1/2/3)
# where the locus suffix prevents a clean name match to mature arms.
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
  # All human miRNAs from bundled miRBase v22 database
  all_mirnas <- getAllMiRNAs(version = "v22", type = "all", species = "hsa")
  setDT(all_mirnas)

  # Separate precursors (MI0...) and matures (MIMAT...)
  pre <- all_mirnas[startsWith(Accession, "MI0"),
                    .(mirna_id = Name, accession = Accession)]
  mat <- all_mirnas[startsWith(Accession, "MIMAT"),
                    .(mature_id = Name, mature_acc = Accession)]

  # miRNA_PrecursorToMature() uses the miRBase database directly, so it
  # correctly maps hsa-mir-9-1 → hsa-miR-9-5p + hsa-miR-9-3p without
  # relying on name-pattern stripping that breaks for locus-suffixed names.
  # Returns columns: OriginalName, Mature1, Mature2  (names only, no accessions)
  mature_raw <- as.data.table(miRNA_PrecursorToMature(pre$mirna_id, version = "v22"))
  setnames(mature_raw, "OriginalName", "mirna_id")

  ann <- merge(pre, mature_raw, by = "mirna_id", all.x = TRUE)

  # Mature1/Mature2 order is not guaranteed to be 5p-first.
  # Assign arms by inspecting the -5p/-3p suffix in the mature name.
  col_or_na <- function(nm) {
    if (nm %in% colnames(ann)) ann[[nm]] else rep(NA_character_, nrow(ann))
  }
  m1 <- col_or_na("Mature1")
  m2 <- col_or_na("Mature2")

  pick_id <- function(id1, id2, pattern) {
    ifelse(!is.na(id1) & grepl(pattern, id1), id1,
    ifelse(!is.na(id2) & grepl(pattern, id2), id2, NA_character_))
  }

  ann[, mature_id_5p := pick_id(m1, m2, "-5p$")]
  ann[, mature_id_3p := pick_id(m1, m2, "-3p$")]

  # Look up accessions for each mature from the matures table
  acc_map        <- mat$mature_acc
  names(acc_map) <- mat$mature_id
  ann[, mature_acc_5p := unname(acc_map[mature_id_5p])]
  ann[, mature_acc_3p := unname(acc_map[mature_id_3p])]

  out <- ann[, .(mirna_id, accession,
                 mature_id_5p, mature_acc_5p,
                 mature_id_3p, mature_acc_3p)]
  setorder(out, mirna_id)

  fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
  message("SUCCESS: wrote ", nrow(out), " precursor miRNA records to ", outfile)

}, error = function(e) {
  message("ERROR in annotate_mirna: ", conditionMessage(e))
  write_empty(outfile)
})
