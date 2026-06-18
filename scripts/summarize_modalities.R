suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  raw_dir  processed_dir  outfile
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: summarize_modalities.R <project> <raw_dir> <processed_dir> <outfile>")
}

project       <- args[1]
raw_dir       <- args[2]
processed_dir <- args[3]
outfile       <- args[4]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
read_col <- function(path, col) {
  if (!file.exists(path)) return(character(0))
  dt <- tryCatch(fread(path, select = col, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(character(0))
  unique(dt[[col]])
}

to_patient <- function(ids) unique(substr(ids, 1, 12))

# ---------------------------------------------------------------------------
# Collect unique sample identifiers per modality
# ---------------------------------------------------------------------------
modalities <- list(
  rna         = file.path(raw_dir, "rna.tsv"),
  mirna       = file.path(raw_dir, "mirna.tsv"),
  methylation = file.path(raw_dir, "methylation.tsv"),
  cnv         = file.path(raw_dir, "cnv.tsv"),
  annotation  = file.path(processed_dir, "annotation.tsv")
)

id_col <- c(rna = "barcode", mirna = "barcode", methylation = "barcode",
            cnv = "barcode", annotation = "barcode")

sample_sets <- setNames(lapply(names(modalities), function(mod) {
  read_col(modalities[[mod]], id_col[[mod]])
}), names(modalities))

# ---------------------------------------------------------------------------
# Per-modality sample counts
# ---------------------------------------------------------------------------
counts <- data.table(
  project  = project,
  modality = names(sample_sets),
  n_samples = vapply(sample_sets, length, integer(1))
)

# ---------------------------------------------------------------------------
# N-way patient-level overlap (2-way through 5-way, for modalities with data)
# ---------------------------------------------------------------------------
active <- names(which(vapply(sample_sets, length, integer(1)) > 0))

if (length(active) >= 2) {
  overlap_rows <- list()
  for (k in 2:min(length(active), 5)) {
    combos <- combn(active, k, simplify = FALSE)
    for (combo in combos) {
      patient_sets <- lapply(combo, function(m) to_patient(sample_sets[[m]]))
      overlap_rows <- c(overlap_rows, list(data.table(
        project   = project,
        modality  = paste(combo, collapse = "_x_"),
        n_samples = length(Reduce(intersect, patient_sets))
      )))
    }
  }
  out <- rbind(counts, rbindlist(overlap_rows), fill = TRUE)
} else {
  out <- counts
}

fwrite(out, outfile, sep = "\t", quote = FALSE)
message("SUCCESS: wrote modality summary for ", project)
