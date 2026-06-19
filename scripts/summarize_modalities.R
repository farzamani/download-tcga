suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: project  raw_dir  outfile
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: summarize_modalities.R <project> <raw_dir> <outfile>")
}

project <- args[1]
raw_dir <- args[2]
outfile <- args[3]

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
# Collect unique barcodes per data modality (annotation excluded — it is a
# derived file, not an independent assay)
# ---------------------------------------------------------------------------
modalities <- list(
  mrna        = file.path(raw_dir, "mrna.tsv"),
  mirna       = file.path(raw_dir, "mirna.tsv"),
  methylation = file.path(raw_dir, "methylation.tsv"),
  cnv         = file.path(raw_dir, "cnv.tsv")
)

sample_sets <- setNames(lapply(names(modalities), function(mod) {
  read_col(modalities[[mod]], "barcode")
}), names(modalities))

# ---------------------------------------------------------------------------
# Per-modality sample counts
# ---------------------------------------------------------------------------
counts <- data.table(
  project   = project,
  modality  = names(sample_sets),
  n_samples = vapply(sample_sets, length, integer(1))
)

# ---------------------------------------------------------------------------
# N-way patient-level overlap (2-way through 4-way, for modalities with data)
# ---------------------------------------------------------------------------
active <- names(which(vapply(sample_sets, length, integer(1)) > 0))

if (length(active) >= 2) {
  overlap_rows <- list()
  for (k in 2:min(length(active), 4)) {
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
