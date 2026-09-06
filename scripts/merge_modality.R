suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

# ---------------------------------------------------------------------------
# Args: format  outfile  [--top-var=N]  file1  file2  ...
#
# format: "wide"  — samples × features matrix (RNA, miRNA, methylation, CNV gene)
#                   metadata columns are the first 5 (barcode … sample_type)
#                   feature columns are intersected across projects so every
#                   merged row has a value for every column
#         "long"  — one observation per row (annotation)
#                   columns are union-merged with NA fill (fill = TRUE)
#
# --top-var=N: (wide only) after merging, keep only the N common feature
#              columns with the highest variance across pooled samples.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

top_var <- NULL
flag_idx <- grep("^--top-var=", args)
if (length(flag_idx) > 0) {
  top_var <- as.integer(sub("^--top-var=", "", args[flag_idx]))
  args <- args[-flag_idx]
}

if (length(args) < 3) {
  stop("Usage: merge_modality.R <wide|long> <outfile> [--top-var=N] <file1> [file2 ...]")
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
    # Use by-reference delete + reorder rather than dt[, ..cols] — with tens
    # of thousands of columns (e.g. full transcriptome), data.table's NSE
    # column-select path overflows R's PROTECT stack ("protection stack
    # overflow"). :=NULL and setcolorder() operate on plain vectors in C
    # and don't hit that limit.
    drop_cols <- setdiff(colnames(dt), keep_cols)
    if (length(drop_cols) > 0) dt[, (drop_cols) := NULL]
    setcolorder(dt, intersect(keep_cols, colnames(dt)))
    dt
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

if (fmt == "wide" && !is.null(top_var) && nrow(out) > 0) {
  meta_present <- intersect(META_COLS, colnames(out))
  feature_cols <- setdiff(colnames(out), META_COLS)
  if (length(feature_cols) > top_var) {
    # setDF() + base data.frame indexing, not dt[, ..cols] — with tens of
    # thousands of columns, data.table's NSE column-select path overflows
    # R's PROTECT stack (see rbindlist note above; same root cause).
    setDF(out)
    mat  <- as.matrix(out[, feature_cols, drop = FALSE])
    vars <- colVars(mat, na.rm = TRUE)
    keep <- feature_cols[order(vars, decreasing = TRUE)[seq_len(top_var)]]
    out  <- out[, c(meta_present, keep), drop = FALSE]
    message("Kept top ", top_var, " variable common features out of ", length(feature_cols))
  } else {
    message("--top-var=", top_var, " requested but only ", length(feature_cols),
            " common features available; keeping all of them")
  }
}

message("Writing ", nrow(out), " rows × ", ncol(out), " columns to ", outfile)
fwrite(out, outfile, sep = "\t", quote = FALSE, na = "NA")
message("Done.")
