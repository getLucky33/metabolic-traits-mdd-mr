#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "psych"))

harm_dir <- arg_required(args, "harm-dir")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
pilot <- isTRUE(args[["pilot"]])
prevalences <- as.numeric(csv_values(arg_value(args, "prevalences", "0.08,0.15,0.20")))
presso_nb <- arg_value(args, "presso-nb", 10000L, "integer")
presso_seed <- arg_value(args, "presso-seed", 20260815L, "integer")
run_presso <- !isTRUE(args[["skip-presso"]])
workers <- arg_value(args, "workers", 4L, "integer")
pleio_threshold <- arg_value(args, "pleio-trait-count", 5L, "integer")
exposure_n_file <- args[["exposure-n-file"]]
exposure_n_max <- arg_value(args, "exposure-n-max", 599249L, "integer")
iv_manifest_file <- args[["iv-manifest"]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(harm_dir, pattern = "^harmonised_main_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("^harmonised_main_a3_", basename(files))]
traits <- sub("^harmonised_main_(.*)\\.rds$", "\\1", basename(files))
if (length(traits_keep)) { keep <- traits %in% safe_trait(traits_keep); files <- files[keep]; traits <- traits[keep] }
if (!length(files)) stop("No harmonised_main_*.rds inputs")
if (anyDuplicated(traits)) stop("Harmonised inputs contain duplicate trait identifiers")
if (!pilot) {
  if (length(traits) != 15L || data.table::uniqueN(traits) != 15L) {
    stop("Formal robustness mode requires exactly 15 unique candidate traits; use --pilot for a subset")
  }
  if (!run_presso || presso_nb < 10000L) {
    stop("Formal robustness mode requires MR-PRESSO with at least 10,000 simulations")
  }
  if (is.null(iv_manifest_file)) {
    stop("Formal robustness mode requires the complete 249-trait --iv-manifest")
  }
  if (is.null(exposure_n_file)) {
    stop("Formal robustness mode requires --exposure-n-file with source-GWAS per-variant N")
  }
}
if (presso_nb < 1L) stop("--presso-nb must be positive")
if (is.na(workers) || workers < 1L || workers > 8L) stop("--workers must be between 1 and 8")
if (is.na(exposure_n_max) || exposure_n_max < 1L) stop("--exposure-n-max must be positive")
if (!length(prevalences) || any(!is.finite(prevalences)) ||
    any(prevalences <= 0 | prevalences >= 1) || anyDuplicated(prevalences)) {
  stop("--prevalences must contain unique finite values strictly between 0 and 1")
}

run_steiger <- !is.null(exposure_n_file)
exposure_n <- NULL
if (run_steiger) {
  exposure_n_file <- assert_file(exposure_n_file, "per-variant exposure N table")
  exposure_n <- data.table::fread(exposure_n_file)
  assert_columns(exposure_n, c("trait", "rsid", "variant_id", "n_exposure"),
                 "per-variant exposure N table")
  if (anyNA(exposure_n[, .(trait, rsid, variant_id)]) ||
      any(!nzchar(exposure_n$trait)) || any(!nzchar(exposure_n$rsid)) ||
      any(!nzchar(exposure_n$variant_id))) {
    stop("Per-variant exposure N keys must be complete and non-empty")
  }
  exposure_n[, trait_key := safe_trait(trait)]
  if (data.table::uniqueN(exposure_n$trait_key) != data.table::uniqueN(exposure_n$trait)) {
    stop("Per-variant exposure N trait identifiers collide after filename sanitization")
  }
  if (anyDuplicated(exposure_n[, .(trait_key, rsid)]) ||
      anyDuplicated(exposure_n[, .(trait_key, variant_id)])) {
    stop("Per-variant exposure N table contains duplicate trait+variant keys")
  }
  exposure_n[, n_exposure := suppressWarnings(as.numeric(n_exposure))]
  if (anyNA(exposure_n$n_exposure) || any(!is.finite(exposure_n$n_exposure)) ||
      any(exposure_n$n_exposure <= 0) || any(exposure_n$n_exposure > exposure_n_max) ||
      any(exposure_n$n_exposure != floor(exposure_n$n_exposure))) {
    stop("Per-variant exposure N values must be positive integers no greater than --exposure-n-max")
  }
  missing_n_traits <- setdiff(traits, unique(exposure_n$trait_key))
  if (length(missing_n_traits)) {
    stop("Requested traits lack per-variant exposure N: ", paste(missing_n_traits, collapse = ", "))
  }
}

pleio_rsids <- hotspot_rsids <- character()
iv_manifest_scope <- "not supplied"
if (!is.null(iv_manifest_file)) {
  ivm <- data.table::fread(assert_file(iv_manifest_file, "IV manifest"))
  assert_columns(ivm, c("trait", "iv_file"), "IV manifest")
  if (anyDuplicated(ivm$trait) || any(!nzchar(ivm$trait))) {
    stop("IV manifest must contain unique, non-empty trait identifiers")
  }
  if (!pilot && (nrow(ivm) != 249L || data.table::uniqueN(ivm$trait) != 249L)) {
    stop("Formal broad-pleiotropy analysis requires the complete 249-trait IV manifest")
  }
  if (!pilot && data.table::uniqueN(safe_trait(ivm$trait)) != 249L) {
    stop("The 249-trait IV manifest has colliding filename-safe trait identifiers")
  }
  if (!pilot && !all(traits %in% safe_trait(ivm$trait))) {
    stop("Every formal candidate must be represented in the 249-trait IV manifest")
  }
  iv_manifest_scope <- if (nrow(ivm) == 249L) "complete 249-trait manifest" else "pilot subset manifest"
  all_iv <- data.table::rbindlist(lapply(seq_len(nrow(ivm)), function(i) {
    z <- data.table::fread(assert_file(ivm$iv_file[[i]], "IV file"), select = c("rsid", "variant_id"))
    z[, trait := ivm$trait[[i]]]
  }))
  count <- unique(all_iv)[, .(n_traits = data.table::uniqueN(trait)), by = rsid]
  pleio_rsids <- count[n_traits > pleio_threshold, rsid]
  data.table::fwrite(count[n_traits > pleio_threshold],
                     file.path(out_dir, paste0(if (pilot) "pilot_" else "", "broad_pleiotropy_snps.tsv")),
                     sep = "\t")
  all_iv[, `:=`(
    chr = suppressWarnings(as.integer(sub("_.*", "", variant_id))),
    pos = suppressWarnings(as.integer(sub("^[^_]*_([0-9]+)_.*", "\\1", variant_id)))
  )]
  hotspots <- data.table::data.table(
    gene = c("FADS", "GCKR", "APOE", "LPA", "CETP"),
    chr = c(11L, 2L, 19L, 6L, 16L),
    start = c(61300000L, 27600000L, 44900000L, 160500000L, 56900000L),
    end = c(61800000L, 27800000L, 45500000L, 161000000L, 57100000L)
  )
  hotspot_hits <- all_iv[hotspots, on = .(chr, pos >= start, pos <= end), nomatch = NULL,
                         .(trait, rsid, variant_id, gene = i.gene)]
  hotspot_rsids <- unique(hotspot_hits$rsid)
  data.table::fwrite(unique(hotspot_hits),
                     file.path(out_dir, paste0(if (pilot) "pilot_" else "", "hotspot_region_snps.tsv")),
                     sep = "\t")
}

ivw <- function(h) {
  h <- h[h$mr_keep %in% TRUE, , drop = FALSE]
  if (nrow(h) < 2L) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  r <- tryCatch(fit_ivw_random(h, "metabolic trait", "major depressive disorder"), error = function(e) NULL)
  if (is.null(r)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  c(b = r@Estimate, se = r@StdError, p = r@Pvalue)
}

checkpoint_dir <- file.path(out_dir, paste0(if (pilot) "pilot_" else "", "presso_checkpoints"))
checkpoint_schema <- 2L
if (run_presso) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
input_hashes <- vapply(files, sha256_file, character(1))
if (run_presso && anyNA(input_hashes)) stop("Could not hash every harmonised input for checkpoint validation")

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite an existing checkpoint: ", path)
  temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Could not finalize checkpoint: ", path)
}

read_checkpoint <- function(path, trait, input_hash) {
  if (!file.exists(path)) return(NULL)
  checkpoint <- tryCatch(readRDS(path), error = function(e) {
    stop("Unreadable MR-PRESSO checkpoint for ", trait, ": ", conditionMessage(e))
  })
  metadata <- checkpoint$metadata
  if (is.null(metadata) || !identical(metadata$schema, checkpoint_schema) ||
      !identical(metadata$trait, trait) ||
      !identical(as.integer(metadata$presso_nb), as.integer(presso_nb)) ||
      !identical(as.integer(metadata$presso_seed), as.integer(presso_seed)) ||
      !identical(metadata$input_sha256, input_hash) || is.null(checkpoint$presso)) {
    stop("MR-PRESSO checkpoint metadata mismatch for ", trait,
         "; move the stale checkpoint before rerunning")
  }
  summary <- data.table::as.data.table(checkpoint$summary)
  if (nrow(summary) != 1L || summary$trait[[1]] != trait) {
    stop("Invalid MR-PRESSO checkpoint summary for ", trait)
  }
  summary
}

presso_columns <- c(
  "trait", "n_iv", "presso_global_p", "presso_global_p_display",
  "presso_global_rssobs", "presso_global_exceedance_n", "presso_outlier_n",
  "presso_distortion_p", "presso_distortion_p_display",
  "presso_distortion_coefficient", "presso_raw_b", "presso_raw_se",
  "presso_raw_df", "presso_raw_t_stat", "presso_raw_ci_lo",
  "presso_raw_ci_hi", "presso_raw_p", "presso_corrected_b",
  "presso_corrected_se", "presso_corrected_df", "presso_corrected_t_stat",
  "presso_corrected_ci_lo", "presso_corrected_ci_hi", "presso_corrected_p",
  "presso_regression_p_sidedness", "presso_empirical_test_df",
  "presso_se_status", "presso_nb", "presso_seed"
)

validate_presso_summary <- function(x, trait, n_iv) {
  assert_columns(x, presso_columns, paste0("MR-PRESSO summary for ", trait))
  if (nrow(x) != 1L || x$trait[[1]] != trait || x$n_iv[[1]] != n_iv ||
      x$presso_nb[[1]] != presso_nb || x$presso_seed[[1]] != presso_seed) {
    stop("MR-PRESSO summary identity or parameter mismatch for ", trait)
  }
  global_display_p <- num_p(x$presso_global_p_display[[1]])
  if (!is.finite(x$presso_global_p[[1]]) || !is.finite(global_display_p) ||
      !is.finite(x$presso_global_rssobs[[1]]) || x$presso_global_rssobs[[1]] < 0 ||
      !isTRUE(all.equal(x$presso_global_p[[1]], global_display_p, tolerance = 1e-12))) {
    stop("MR-PRESSO global reporting fields are incomplete or inconsistent for ", trait)
  }
  global_is_bound <- startsWith(x$presso_global_p_display[[1]], "<")
  if ((global_is_bound && !identical(as.integer(x$presso_global_exceedance_n[[1]]), 0L)) ||
      (!global_is_bound && !is.na(x$presso_global_exceedance_n[[1]]))) {
    stop("MR-PRESSO global exceedance count is inconsistent for ", trait)
  }
  if (!identical(x$presso_regression_p_sidedness[[1]], "two-sided") ||
      !identical(x$presso_empirical_test_df[[1]], "not_applicable")) {
    stop("MR-PRESSO test metadata are invalid for ", trait)
  }
  distortion_display_p <- num_p(x$presso_distortion_p_display[[1]])
  if (is.finite(x$presso_distortion_p[[1]])) {
    if (!is.finite(distortion_display_p) ||
        !isTRUE(all.equal(x$presso_distortion_p[[1]], distortion_display_p,
                          tolerance = 1e-12)) ||
        !is.finite(x$presso_distortion_coefficient[[1]])) {
      stop("MR-PRESSO distortion reporting fields are incomplete or inconsistent for ", trait)
    }
  } else if (!is.na(x$presso_distortion_p_display[[1]]) ||
             is.finite(x$presso_distortion_coefficient[[1]])) {
    stop("MR-PRESSO distortion reporting fields are inconsistent for ", trait)
  }
  if (!is.finite(x$presso_raw_b[[1]]) || !is.finite(x$presso_raw_se[[1]]) ||
      x$presso_raw_se[[1]] <= 0 || x$presso_raw_df[[1]] != n_iv - 1L ||
      !is.finite(x$presso_raw_t_stat[[1]]) ||
      !is.finite(x$presso_raw_ci_lo[[1]]) || !is.finite(x$presso_raw_ci_hi[[1]])) {
    stop("MR-PRESSO raw estimate did not retain valid Sd/t-interval fields for ", trait)
  }
  raw_critical <- stats::qt(0.975, df = x$presso_raw_df[[1]])
  if (!isTRUE(all.equal(x$presso_raw_ci_lo[[1]],
                        x$presso_raw_b[[1]] - raw_critical * x$presso_raw_se[[1]],
                        tolerance = 1e-10)) ||
      !isTRUE(all.equal(x$presso_raw_ci_hi[[1]],
                        x$presso_raw_b[[1]] + raw_critical * x$presso_raw_se[[1]],
                        tolerance = 1e-10))) {
    stop("MR-PRESSO raw t interval is inconsistent for ", trait)
  }
  raw_p_check <- 2 * stats::pt(
    -abs(x$presso_raw_b[[1]] / x$presso_raw_se[[1]]),
    df = x$presso_raw_df[[1]]
  )
  if (!isTRUE(all.equal(x$presso_raw_p[[1]], raw_p_check, tolerance = 1e-10))) {
    stop("MR-PRESSO raw P-value/Sd consistency check failed for ", trait)
  }
  if (!isTRUE(all.equal(x$presso_raw_t_stat[[1]],
                        x$presso_raw_b[[1]] / x$presso_raw_se[[1]],
                        tolerance = 1e-10))) {
    stop("MR-PRESSO raw T-stat/Sd consistency check failed for ", trait)
  }
  if (is.finite(x$presso_corrected_b[[1]])) {
    expected_df <- n_iv - x$presso_outlier_n[[1]] - 1L
    corrected_critical <- stats::qt(0.975, df = expected_df)
    if (!is.finite(x$presso_corrected_se[[1]]) || x$presso_corrected_se[[1]] <= 0 ||
        !is.finite(x$presso_corrected_t_stat[[1]]) ||
        x$presso_corrected_df[[1]] != expected_df || expected_df < 1L ||
        !isTRUE(all.equal(x$presso_corrected_ci_lo[[1]],
                          x$presso_corrected_b[[1]] - corrected_critical * x$presso_corrected_se[[1]],
                          tolerance = 1e-10)) ||
        !isTRUE(all.equal(x$presso_corrected_ci_hi[[1]],
                          x$presso_corrected_b[[1]] + corrected_critical * x$presso_corrected_se[[1]],
                          tolerance = 1e-10))) {
      stop("MR-PRESSO corrected estimate did not retain valid Sd/t-interval fields for ", trait)
    }
    corrected_p_check <- 2 * stats::pt(
      -abs(x$presso_corrected_b[[1]] / x$presso_corrected_se[[1]]),
      df = x$presso_corrected_df[[1]]
    )
    if (!isTRUE(all.equal(x$presso_corrected_p[[1]], corrected_p_check,
                          tolerance = 1e-10))) {
      stop("MR-PRESSO corrected P-value/Sd consistency check failed for ", trait)
    }
    if (!isTRUE(all.equal(x$presso_corrected_t_stat[[1]],
                          x$presso_corrected_b[[1]] / x$presso_corrected_se[[1]],
                          tolerance = 1e-10))) {
      stop("MR-PRESSO corrected T-stat/Sd consistency check failed for ", trait)
    }
  } else if (is.finite(x$presso_corrected_t_stat[[1]])) {
    stop("MR-PRESSO corrected T statistic exists without a corrected estimate for ", trait)
  }
  invisible(TRUE)
}

run_one <- function(i) {
  tr <- traits[[i]]
  h <- readRDS(files[[i]])
  h <- h[h$mr_keep %in% TRUE, , drop = FALSE]
  if (!nrow(h) || anyDuplicated(h$SNP)) stop("Invalid or duplicate harmonised SNPs for ", tr)
  full <- ivw(h)
  keep_broad <- if (length(pleio_rsids)) h[!h$SNP %in% pleio_rsids, , drop = FALSE] else h
  keep_hotspot <- if (length(hotspot_rsids)) h[!h$SNP %in% hotspot_rsids, , drop = FALSE] else h
  after_broad <- ivw(keep_broad)
  after_hotspot <- ivw(keep_hotspot)
  row <- data.table::data.table(
    trait = tr, n_iv_full = nrow(h),
    n_removed_broad_pleiotropy = nrow(h) - nrow(keep_broad),
    n_removed_hotspot = nrow(h) - nrow(keep_hotspot),
    full_b = full[["b"]], full_se = full[["se"]],
    full_ci_lo = full[["b"]] - 1.96 * full[["se"]],
    full_ci_hi = full[["b"]] + 1.96 * full[["se"]], full_p = full[["p"]],
    after_broad_b = after_broad[["b"]], after_broad_se = after_broad[["se"]],
    after_broad_ci_lo = after_broad[["b"]] - 1.96 * after_broad[["se"]],
    after_broad_ci_hi = after_broad[["b"]] + 1.96 * after_broad[["se"]],
    after_broad_p = after_broad[["p"]],
    after_broad_direction_same = sign(full[["b"]]) == sign(after_broad[["b"]]),
    after_hotspot_b = after_hotspot[["b"]], after_hotspot_se = after_hotspot[["se"]],
    after_hotspot_ci_lo = after_hotspot[["b"]] - 1.96 * after_hotspot[["se"]],
    after_hotspot_ci_hi = after_hotspot[["b"]] + 1.96 * after_hotspot[["se"]],
    after_hotspot_p = after_hotspot[["p"]],
    after_hotspot_direction_same = sign(full[["b"]]) == sign(after_hotspot[["b"]])
  )

  if (run_steiger) {
    required <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome",
                  "eaf.outcome", "ncase", "ncontrol")
    if (!all(required %in% names(h))) {
      stop("Steiger input is missing columns for ", tr, ": ",
           paste(setdiff(required, names(h)), collapse = ", "))
    }
    n_map <- exposure_n[trait_key == tr]
    matched <- match(h$SNP, n_map$rsid)
    if (!nrow(n_map) || anyNA(matched)) {
      stop("Per-variant exposure N match is incomplete for ", tr,
           "; missing harmonised SNPs: ", sum(is.na(matched)))
    }
    n_exp <- n_map$n_exposure[matched]
    if (anyNA(h$ncase) || anyNA(h$ncontrol) ||
        any(!is.finite(h$ncase)) || any(!is.finite(h$ncontrol)) ||
        any(h$ncase <= 0) || any(h$ncontrol <= 0)) {
      stop("Invalid per-variant outcome case/control counts for Steiger analysis: ", tr)
    }
    h$samplesize.exposure <- n_exp
    h$r.exposure <- TwoSampleMR::get_r_from_bsen(
      h$beta.exposure, h$se.exposure, h$samplesize.exposure
    )
    h$samplesize.outcome <- h$ncase + h$ncontrol
    if (any(!is.finite(h$r.exposure)) || any(!is.finite(h$samplesize.outcome))) {
      stop("Non-finite exposure R or outcome sample size for Steiger analysis: ", tr)
    }
    row[, `:=`(
      exposure_n_mode = "source_variant_specific",
      exposure_n_min = min(n_exp), exposure_n_max = max(n_exp),
      exposure_n_unique = data.table::uniqueN(n_exp),
      steiger_test_distribution = "standard_normal",
      steiger_p_sidedness = "two-sided",
      steiger_df = "not_applicable",
      steiger_p_display = NA_character_
    )]
    steiger_p_values <- numeric(length(prevalences))
    for (k in seq_along(prevalences)) {
      prev <- prevalences[[k]]
      hh <- h
      hh$r.outcome <- TwoSampleMR::get_r_from_lor(
        hh$beta.outcome, hh$eaf.outcome,
        ncase = hh$ncase, ncontrol = hh$ncontrol, prevalence = prev
      )
      if (any(!is.finite(hh$r.outcome))) {
        stop("Non-finite outcome R for Steiger analysis: ", tr, " at prevalence ", prev)
      }
      st <- TwoSampleMR::directionality_test(hh)
      needed <- c("correct_causal_direction", "steiger_pval", "snp_r2.exposure", "snp_r2.outcome")
      if (nrow(st) != 1L || !all(needed %in% names(st)) ||
          anyNA(st[, needed, drop = FALSE])) {
        stop("Incomplete Steiger result for ", tr, " at prevalence ", prev)
      }
      r_exp <- sqrt(sum(abs(hh$r.exposure)^2))
      r_out <- sqrt(sum(abs(hh$r.outcome)^2))
      r_test <- psych::r.test(
        n = mean(hh$samplesize.exposure), n2 = mean(hh$samplesize.outcome),
        r12 = r_exp, r34 = r_out
      )
      z_stat <- as.numeric(r_test[["z"]])
      if (length(z_stat) != 1L || !is.finite(z_stat)) {
        stop("Non-finite Steiger z statistic for ", tr, " at prevalence ", prev)
      }
      p_check <- 2 * stats::pnorm(-abs(z_stat))
      log10_p <- (log(2) + stats::pnorm(-abs(z_stat), log.p = TRUE)) / log(10)
      if (!is.finite(log10_p) ||
          !isTRUE(all.equal(p_check, st$steiger_pval[[1]], tolerance = 1e-12))) {
        stop("Steiger z/P consistency check failed for ", tr, " at prevalence ", prev)
      }
      suffix <- gsub("\\.", "", sprintf("%0.2f", prev))
      row[[paste0("steiger_correct_", suffix)]] <- st$correct_causal_direction[[1]]
      row[[paste0("steiger_p_", suffix)]] <- st$steiger_pval[[1]]
      row[[paste0("steiger_z_", suffix)]] <- z_stat
      row[[paste0("steiger_log10_p_", suffix)]] <- log10_p
      steiger_p_values[[k]] <- st$steiger_pval[[1]]
      if (abs(prev - 0.15) < 1e-12) {
        row[["r2_exposure"]] <- st$snp_r2.exposure[[1]]
        row[["r2_outcome"]] <- st$snp_r2.outcome[[1]]
      }
    }
    if (length(steiger_p_values) == 3L && all(steiger_p_values == 0)) {
      row[, steiger_p_display := "P<1e-300 (double-precision underflow; see log10 P)"]
    }
  }

  resumed <- FALSE
  if (run_presso) {
    checkpoint_file <- file.path(checkpoint_dir, paste0(tr, ".rds"))
    presso_row <- read_checkpoint(checkpoint_file, tr, input_hashes[[i]])
    if (is.null(presso_row)) {
      suite <- run_mr_suite(h, tr, "major depressive disorder", presso_nb, presso_seed, TRUE)
      if (is.null(suite$presso)) {
        detail <- if (is.na(suite$presso_error)) "no result returned" else suite$presso_error
        stop("MR-PRESSO failed for ", tr, ": ", detail)
      }
      presso_row <- summarise_suite(suite, tr, presso_nb)[, ..presso_columns]
      validate_presso_summary(presso_row, tr, nrow(h))
      atomic_save_rds(list(
        metadata = list(schema = checkpoint_schema, trait = tr, presso_nb = presso_nb,
                        presso_seed = presso_seed, input_sha256 = input_hashes[[i]],
                        completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
        summary = as.data.frame(presso_row), presso = suite$presso
      ), checkpoint_file)
    } else {
      resumed <- TRUE
      validate_presso_summary(presso_row, tr, nrow(h))
    }
  } else {
    suite <- run_mr_suite(h, tr, "major depressive disorder", presso_nb, presso_seed, FALSE)
    presso_row <- summarise_suite(suite, tr, presso_nb)[, ..presso_columns]
  }
  presso_row[, checkpoint_resumed := resumed]
  list(robust = row, presso = presso_row)
}

effective_workers <- min(workers, length(files))
if (effective_workers == 1L) {
  results <- lapply(seq_along(files), run_one)
} else {
  require_packages(c("future", "future.apply"))
  old_plan <- future::plan()
  results <- tryCatch({
    future::plan(future::multisession, workers = effective_workers)
    future.apply::future_lapply(
      seq_along(files), run_one, future.seed = TRUE,
      future.packages = c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "psych")
    )
  }, finally = future::plan(old_plan))
}
robust_out <- data.table::rbindlist(lapply(results, `[[`, "robust"), fill = TRUE)
presso_out <- data.table::rbindlist(lapply(results, `[[`, "presso"), fill = TRUE)
data.table::setorder(robust_out, trait)
data.table::setorder(presso_out, trait)
if (!pilot && (nrow(robust_out) != 15L || nrow(presso_out) != 15L ||
               data.table::uniqueN(robust_out$trait) != 15L ||
               data.table::uniqueN(presso_out$trait) != 15L ||
               !setequal(robust_out$trait, presso_out$trait))) {
  stop("Formal robustness outputs must each contain all 15 candidate traits")
}
prefix <- if (pilot) "pilot_" else ""
data.table::fwrite(robust_out, file.path(out_dir, paste0(prefix, "robustness_summary.tsv")), sep = "\t")
data.table::fwrite(
  presso_out,
  file.path(out_dir, paste0(prefix, "presso_rerun_", presso_nb, ".tsv")),
  sep = "\t"
)
metadata_inputs <- setNames(as.list(files), paste0("harmonised_", traits))
if (run_steiger) metadata_inputs$exposure_variant_n <- exposure_n_file
write_run_metadata(out_dir, "robustness", metadata_inputs, list(
  prevalences = paste(prevalences, collapse = ","), presso_nb = presso_nb,
  presso_seed = presso_seed, presso_run = run_presso, workers = workers,
  checkpoint_resume = if (run_presso) "automatic, validated by parameters and input SHA-256" else "not applicable",
  checkpoint_resumed_n = sum(presso_out$checkpoint_resumed),
  analysis_mode = if (pilot) "pilot" else "formal_15_trait",
  iv_manifest_scope = iv_manifest_scope,
  broad_pleiotropy_definition = paste0(">", pleio_threshold, " traits"),
  steiger_run = run_steiger,
  exposure_n_min = if (run_steiger) min(robust_out$exposure_n_min) else NA_real_,
  exposure_n_max = if (run_steiger) max(robust_out$exposure_n_max) else NA_real_,
  steiger_exposure_n_note = if (run_steiger) {
    "source-GWAS per-variant N; incomplete or invalid matching fails closed"
  } else {
    "not requested; permitted only in pilot mode"
  }
))
