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

id_col <- c(rna = "sample_id", mirna = "sample_id", methylation = "sample_id",
            cnv = "barcode", annotation = "patient_id")

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
# Pairwise patient-level overlap (only for modalities with data)
# ---------------------------------------------------------------------------
active <- names(which(vapply(sample_sets, length, integer(1)) > 0))

if (length(active) >= 2) {
  pairs <- combn(active, 2, simplify = FALSE)
  overlap_rows <- lapply(pairs, function(p) {
    data.table(
      project   = project,
      modality  = paste(p[1], p[2], sep = "_x_"),
      n_samples = length(intersect(to_patient(sample_sets[[p[1]]]),
                                   to_patient(sample_sets[[p[2]]]))
      )
    )
  })
  out <- rbind(counts, rbindlist(overlap_rows), fill = TRUE)
} else {
  out <- counts
}

fwrite(out, outfile, sep = "\t", quote = FALSE)
message("SUCCESS: wrote modality summary for ", project)
