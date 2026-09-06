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
