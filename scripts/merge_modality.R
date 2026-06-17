suppressPackageStartupMessages({
  library(data.table)
})

# ---------------------------------------------------------------------------
# Args: format  outfile  file1  file2  ...
#
# format: "wide"  — samples × features matrix (RNA, miRNA, methylation)
#                   metadata columns are the first 5 (barcode … sample_type)
#                   feature columns are intersected across projects so every
#                   merged row has a value for every column
#         "long"  — one observation per row (CNV, annotation)
#                   columns are union-merged with NA fill (fill = TRUE)
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: merge_modality.R <wide|long> <outfile> <file1> [file2 ...]")
}

fmt     <- args[1]
outfile <- args[2]
infiles <- args[-(1:2)]

META_COLS <- c("barcode", "patient_id", "sample_id", "project", "sample_type")

# ---------------------------------------------------------------------------
# Helper: read one file, skip if empty
# ---------------------------------------------------------------------------
read_one <- function(path) {
  if (!file.exists(path)) {
    message("SKIP (not found): ", path)
    return(NULL)
  }
  dt <- tryCatch(
    fread(path, showProgress = FALSE),
    error = function(e) { message("SKIP (read error): ", path, " — ", conditionMessage(e)); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0) {
    message("SKIP (empty): ", path)
    return(NULL)
  }
  dt
}

# ---------------------------------------------------------------------------
# Wide merge: intersect feature columns, rbind
# ---------------------------------------------------------------------------
merge_wide <- function(files) {
  tables <- lapply(files, read_one)
  tables <- Filter(Negate(is.null), tables)

  if (length(tables) == 0) {
    message("WARNING: all input files empty; writing header-only output")
    return(data.table())
  }

  # Feature columns present in each file (excluding metadata)
  feature_sets <- lapply(tables, function(dt) {
    setdiff(colnames(dt), META_COLS)
  })

  # Intersection: only features measured in every project
  common_features <- Reduce(intersect, feature_sets)
  if (length(common_features) == 0) {
    message("WARNING: no common features across projects; writing metadata only")
    common_features <- character(0)
  } else {
    message("Common features retained: ", length(common_features),
            " (out of max ", max(vapply(feature_sets, length, integer(1))), ")")
  }

  keep_cols <- c(META_COLS, common_features)
  tables    <- lapply(tables, function(dt) {
    existing <- intersect(keep_cols, colnames(dt))
    dt[, ..existing]
  })

  rbindlist(tables, use.names = TRUE, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Long merge: simple row-bind with NA fill for any missing columns
# ---------------------------------------------------------------------------
merge_long <- function(files) {
  tables <- lapply(files, read_one)
  tables <- Filter(Negate(is.null), tables)

  if (length(tables) == 0) {
    message("WARNING: all input files empty; writing header-only output")
    return(data.table())
  }

  rbindlist(tables, use.names = TRUE, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
out <- if (fmt == "wide") {
  merge_wide(infiles)
} else if (fmt == "long") {
  merge_long(infiles)
} else {
  stop("Unknown format '", fmt, "'; expected 'wide' or 'long'")
}

message("Writing ", nrow(out), " rows × ", ncol(out), " columns to ", outfile)
fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
message("Done.")
