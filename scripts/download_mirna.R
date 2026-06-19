suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(miRBaseConverter)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  sample_type  outfile  gdc_cache
# Output columns are mature arm IDs where available (hsa-miR-9-5p, hsa-miR-9-3p)
# rather than precursor IDs (hsa-mir-9-1, hsa-mir-9-2, hsa-mir-9-3).
# When multiple precursors produce the same mature arm, their counts are summed.
# Precursors with no arm annotation in miRBase v22 retain their original ID.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: download_mirna.R <project> <sample_type> <outfile> <gdc_cache>")

project     <- args[1]
sample_type <- args[2]
abs_path  <- function(p) if (startsWith(p, "/")) p else file.path(getwd(), p)
outfile   <- abs_path(args[3])
gdc_cache <- abs_path(args[4])

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

  # GDCprepare returns a wide data.frame:
  # columns: miRNA_ID | read_count_<barcode> | reads_per_million_..._<barcode> | cross-mapped_<barcode>
  df <- GDCprepare(query, directory = gdc_cache)

  rc_cols  <- grep("^read_count_", colnames(df), value = TRUE)
  if (length(rc_cols) == 0) stop("No read_count_ columns found in GDCprepare output")

  barcodes  <- sub("^read_count_", "", rc_cols)
  count_mat <- as.matrix(df[, rc_cols, drop = FALSE])
  rownames(count_mat) <- df$miRNA_ID   # precursor IDs, e.g. hsa-mir-9-1
  colnames(count_mat) <- barcodes
  storage.mode(count_mat) <- "numeric"

  # ------------------------------------------------------------------
  # Map precursor IDs → mature arm names using miRBase v22 database.
  # Uses miRNA_PrecursorToMature() which handles multi-locus families
  # correctly (hsa-mir-9-1/2/3 → hsa-miR-9-5p + hsa-miR-9-3p).
  # ------------------------------------------------------------------
  all_pre <- rownames(count_mat)
  mature_map <- as.data.table(miRNA_PrecursorToMature(all_pre, version = "v22"))
  setnames(mature_map, "OriginalName", "mirna_id")

  col_or_na <- function(dt, nm) {
    if (nm %in% colnames(dt)) dt[[nm]] else rep(NA_character_, nrow(dt))
  }
  m1 <- col_or_na(mature_map, "Mature1")
  m2 <- col_or_na(mature_map, "Mature2")

  pick_arm <- function(id1, id2, pattern) {
    ifelse(!is.na(id1) & grepl(pattern, id1), id1,
    ifelse(!is.na(id2) & grepl(pattern, id2), id2, NA_character_))
  }
  mature_map[, arm_5p := pick_arm(m1, m2, "-5p$")]
  mature_map[, arm_3p := pick_arm(m1, m2, "-3p$")]

  # Build feature map: each precursor → one or more output column names
  feat_5p   <- mature_map[!is.na(arm_5p), .(mirna_id, feature = arm_5p)]
  feat_3p   <- mature_map[!is.na(arm_3p), .(mirna_id, feature = arm_3p)]
  feat_none <- mature_map[is.na(arm_5p) & is.na(arm_3p), .(mirna_id, feature = mirna_id)]
  unmapped  <- setdiff(all_pre, mature_map$mirna_id)
  if (length(unmapped) > 0)
    feat_none <- rbind(feat_none, data.table(mirna_id = unmapped, feature = unmapped))
  feat_map <- rbindlist(list(feat_5p, feat_3p, feat_none), fill = TRUE)

  # Sum counts from all precursors contributing to the same mature arm
  all_feats <- unique(feat_map$feature)
  out_mat <- matrix(0, nrow = ncol(count_mat), ncol = length(all_feats),
                    dimnames = list(barcodes, all_feats))
  for (feat in all_feats) {
    pres <- intersect(feat_map[feature == feat, mirna_id], rownames(count_mat))
    if (length(pres) > 0)
      out_mat[, feat] <- colSums(count_mat[pres, , drop = FALSE], na.rm = TRUE)
  }

  out <- as.data.frame(out_mat)
  out <- cbind(data.frame(barcode = rownames(out_mat), stringsAsFactors = FALSE), out)

  fwrite(out, outfile, sep = "\t", quote = FALSE)
  message("SUCCESS: wrote ", nrow(out), " samples x ", length(all_feats),
          " miRNA features (arm-level) to ", outfile)

}, error = function(e) {
  message("ERROR in download_mirna for ", project, ": ", conditionMessage(e))
  write_empty(outfile)
})
