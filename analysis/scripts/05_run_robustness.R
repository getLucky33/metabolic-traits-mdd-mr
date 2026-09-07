#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO"))

harm_dir <- arg_required(args, "harm-dir")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
pilot <- isTRUE(args[["pilot"]])
prevalences <- as.numeric(csv_values(arg_value(args, "prevalences", "0.08,0.15,0.20")))
presso_nb <- arg_value(args, "presso-nb", 10000L, "integer")
presso_seed <- arg_value(args, "presso-seed", 20260815L, "integer")
run_presso <- !isTRUE(args[["skip-presso"]])
pleio_threshold <- arg_value(args, "pleio-trait-count", 5L, "integer")
exposure_n_approx <- arg_value(args, "exposure-n", 599249L, "integer")
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
}
if (presso_nb < 1L) stop("--presso-nb must be positive")

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

robust <- presso <- list()
for (i in seq_along(files)) {
  tr <- traits[[i]]
  h <- readRDS(files[[i]])
  h <- h[h$mr_keep %in% TRUE, , drop = FALSE]
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
  if (all(c("ncase", "ncontrol") %in% names(h))) {
    if (anyNA(h$ncase) || anyNA(h$ncontrol) || any(h$ncase <= 0) || any(h$ncontrol <= 0)) {
      stop("Invalid ncase/ncontrol for Steiger analysis: ", tr)
    }
    h$samplesize.exposure <- exposure_n_approx
    h$r.exposure <- TwoSampleMR::get_r_from_bsen(
      h$beta.exposure, h$se.exposure, exposure_n_approx
    )
    h$samplesize.outcome <- h$ncase + h$ncontrol
    for (prev in prevalences) {
      hh <- h
      hh$r.outcome <- tryCatch(
        TwoSampleMR::get_r_from_lor(hh$beta.outcome, hh$eaf.outcome,
                                    ncase = hh$ncase, ncontrol = hh$ncontrol, prevalence = prev),
        error = function(e) rep(NA_real_, nrow(hh))
      )
      st <- tryCatch(TwoSampleMR::directionality_test(hh), error = function(e) NULL)
      suffix <- gsub("\\.", "", sprintf("%0.2f", prev))
      row[[paste0("steiger_correct_", suffix)]] <- if (is.null(st)) NA else st$correct_causal_direction[[1]]
      row[[paste0("steiger_p_", suffix)]] <- if (is.null(st)) NA_real_ else st$steiger_pval[[1]]
      if (abs(prev - 0.15) < 1e-12) {
        row[["r2_exposure"]] <- if (is.null(st)) NA_real_ else st$snp_r2.exposure[[1]]
        row[["r2_outcome"]] <- if (is.null(st)) NA_real_ else st$snp_r2.outcome[[1]]
      }
    }
  }
  robust[[i]] <- row

  s <- run_mr_suite(h, tr, "major depressive disorder", presso_nb, presso_seed, run_presso)
  presso[[i]] <- summarise_suite(s, tr, presso_nb)[, .(
    trait, n_iv, presso_global_p, presso_outlier_n, presso_distortion_p,
    presso_raw_b, presso_raw_se, presso_raw_ci_lo, presso_raw_ci_hi, presso_raw_p,
    presso_corrected_b, presso_corrected_se, presso_corrected_ci_lo,
    presso_corrected_ci_hi, presso_corrected_p, presso_se_status, presso_nb
  )]
}
robust_out <- data.table::rbindlist(robust, fill = TRUE)
presso_out <- data.table::rbindlist(presso, fill = TRUE)
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
write_run_metadata(out_dir, "robustness", setNames(as.list(files), paste0("harmonised_", traits)), list(
  prevalences = paste(prevalences, collapse = ","), presso_nb = presso_nb,
  presso_seed = presso_seed, presso_run = run_presso,
  analysis_mode = if (pilot) "pilot" else "formal_15_trait",
  iv_manifest_scope = iv_manifest_scope,
  broad_pleiotropy_definition = paste0(">", pleio_threshold, " traits"),
  exposure_sample_size_approximation = exposure_n_approx,
  steiger_exposure_n_note = "fixed maximum meta-analysis N; per-SNP N unavailable"
))
