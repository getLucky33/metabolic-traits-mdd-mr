args_map <- function(x = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- x[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

arg_required <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || identical(value, "")) stop("Missing required --", key)
  value
}

arg_value <- function(args, key, default = NULL, mode = c("character", "integer", "numeric", "logical")) {
  mode <- match.arg(mode)
  value <- args[[key]]
  if (is.null(value)) return(default)
  switch(mode,
    character = as.character(value),
    integer = as.integer(value),
    numeric = as.numeric(value),
    logical = tolower(as.character(value)) %in% c("1", "true", "yes", "y")
  )
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Restore the recorded environment before running this stage.")
  }
}

assert_file <- function(path, label = "input") {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

assert_columns <- function(x, required, label = "table") {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

assert_numeric <- function(x, columns, label = "table") {
  for (column in columns) {
    value <- x[[column]]
    if (!is.numeric(value)) stop(label, "$", column, " must be numeric")
    if (anyNA(value) || any(!is.finite(value))) stop(label, "$", column, " contains missing or non-finite values")
  }
  invisible(TRUE)
}

classify_full_screen <- function(x, p_column = "ivw_p", n_tests = 249L) {
  if (!"trait" %in% names(x)) stop("Full screen requires a trait column")
  if (!p_column %in% names(x)) stop("Full screen is missing ", p_column)
  if (nrow(x) != n_tests || data.table::uniqueN(x$trait) != n_tests) {
    stop("Full screen requires exactly ", n_tests, " unique traits")
  }
  p <- x[[p_column]]
  if (!is.numeric(p) || anyNA(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("Full-screen P values must be finite and within [0,1]")
  }
  x[, p_fdr_bh := stats::p.adjust(get(p_column), method = "BH", n = n_tests)]
  x[, p_bonf := stats::p.adjust(get(p_column), method = "bonferroni", n = n_tests)]
  x[, screen_level := data.table::fcase(
    p_bonf < 0.05, "bonferroni_hit",
    p_fdr_bh < 0.05, "fdr_only",
    get(p_column) < 0.05, "nominal",
    default = "null"
  )]
  x
}

window_overlaps_mhc <- function(chr, window_start, window_end,
                                mhc_start = 25000000, mhc_end = 34000000) {
  chr <- suppressWarnings(as.integer(chr))
  window_start <- suppressWarnings(as.numeric(window_start))
  window_end <- suppressWarnings(as.numeric(window_end))
  valid <- is.finite(chr) & is.finite(window_start) & is.finite(window_end) &
    window_start <= window_end
  valid & chr == 6L & window_end >= mhc_start & window_start <= mhc_end
}

i2gx_from_bse <- function(beta, se) {
  keep <- is.finite(beta) & is.finite(se) & se > 0
  beta <- as.numeric(beta[keep])
  se <- as.numeric(se[keep])
  k <- length(beta)
  if (k < 2L) return(NA_real_)
  w <- 1 / se^2
  mean_w <- sum(w * beta) / sum(w)
  qx <- sum(w * (beta - mean_w)^2)
  if (!is.finite(qx) || qx <= 0) return(0)
  max(0, (qx - (k - 1L)) / qx)
}

safe_trait <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]+", "_", x)
  ifelse(nzchar(y), y, "trait")
}

csv_values <- function(x) {
  if (is.null(x) || identical(x, "")) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) return(NA_character_)
  cmd <- if (.Platform$OS.type == "windows") "certutil" else "sha256sum"
  if (.Platform$OS.type == "windows") {
    z <- system2(cmd, c("-hashfile", shQuote(path), "SHA256"), stdout = TRUE, stderr = TRUE)
    hit <- grep("^[0-9A-Fa-f]{64}$", trimws(z), value = TRUE)
    if (!length(hit)) return(NA_character_)
    tolower(trimws(hit[[1]]))
  } else {
    z <- system2(cmd, shQuote(path), stdout = TRUE, stderr = TRUE)
    if (!length(z)) return(NA_character_)
    strsplit(z[[1]], "[[:space:]]+")[[1]][[1]]
  }
}

write_run_metadata <- function(out_dir, stage, inputs, parameters = list()) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  input_paths <- unlist(inputs, use.names = TRUE)
  input_meta <- data.frame(
    role = names(input_paths),
    file = basename(input_paths),
    sha256 = vapply(input_paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.table(input_meta, file.path(out_dir, paste0(stage, "_input_hashes.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  param_meta <- data.frame(parameter = names(parameters), value = unlist(parameters, use.names = FALSE),
                           stringsAsFactors = FALSE)
  utils::write.table(param_meta, file.path(out_dir, paste0(stage, "_parameters.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
}
