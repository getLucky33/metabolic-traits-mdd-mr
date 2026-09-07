#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages("data.table")

manifest_file <- assert_file(arg_required(args, "iv-manifest"), "IV manifest")
out_file <- arg_required(args, "out")
screen_file <- args[["screen"]]
selection_qc_file <- args[["selection-qc"]]
trait_labels_file <- assert_file(
  arg_value(args, "trait-labels", "data/derived/trait_display_dictionary.tsv"),
  "trait display dictionary"
)
expected_traits <- arg_value(args, "expected-traits", 249L, "integer")
f_min <- arg_value(args, "f-min", 10, "numeric")

manifest <- data.table::fread(manifest_file)
assert_columns(manifest, c("trait", "iv_file"), "IV manifest")
if (nrow(manifest) != expected_traits || data.table::uniqueN(manifest$trait) != expected_traits) {
  stop("IV manifest requires exactly ", expected_traits, " unique traits")
}

rows <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  trait <- manifest$trait[[i]]
  iv <- data.table::fread(assert_file(manifest$iv_file[[i]], paste0("IV file for ", trait)))
  assert_columns(iv, c("beta", "standard_error", "effect_allele_frequency"), trait)
  assert_numeric(iv, c("beta", "standard_error", "effect_allele_frequency"), trait)
  if (any(iv$standard_error <= 0) ||
      any(iv$effect_allele_frequency <= 0 | iv$effect_allele_frequency >= 1)) {
    stop(trait, ": invalid SE or effect-allele frequency")
  }
  f_stat <- (iv$beta / iv$standard_error)^2
  if ("F" %in% names(iv) && any(abs(iv$F - f_stat) > pmax(1e-8, abs(f_stat) * 1e-8))) {
    stop(trait, ": stored and recalculated F statistics differ")
  }
  maf <- pmin(iv$effect_allele_frequency, 1 - iv$effect_allele_frequency)
  r2_approx <- 2 * iv$effect_allele_frequency *
    (1 - iv$effect_allele_frequency) * iv$beta^2
  q <- stats::quantile(maf, c(0, 0.25, 0.5, 0.75, 1), names = FALSE, type = 7)
  i2gx <- i2gx_from_bse(iv$beta, iv$standard_error)
  rows[[i]] <- data.table::data.table(
    trait = trait,
    n_iv = nrow(iv),
    min_f = min(f_stat),
    median_f = stats::median(f_stat),
    cumulative_r2_approx = sum(r2_approx),
    maf_min = q[[1]], maf_q25 = q[[2]], maf_median = q[[3]],
    maf_q75 = q[[4]], maf_max = q[[5]],
    i2gx = i2gx,
    nome_i2gx_ge_0_90 = if (is.finite(i2gx)) i2gx >= 0.90 else NA
  )
}
out <- data.table::rbindlist(rows)
if (any(out$min_f < f_min)) stop("At least one final instrument has F below ", f_min)

labels <- data.table::fread(trait_labels_file)
assert_columns(labels, c("trait_id", "trait_display"), "trait display dictionary")
if (anyDuplicated(labels$trait_id) || anyDuplicated(labels$trait_display) ||
    anyNA(labels$trait_id) || anyNA(labels$trait_display)) {
  stop("Trait display dictionary has missing or duplicate identifiers")
}
out <- merge(labels[, .(trait = trait_id, trait_display)], out, by = "trait", all.y = TRUE)
if (anyNA(out$trait_display)) stop("Trait display name is missing for at least one IV trait")
data.table::setnames(out, "trait", "trait_id")

if (!is.null(screen_file)) {
  screen <- data.table::fread(assert_file(screen_file, "screen table"))
  assert_columns(screen, c("trait", "screen_level"), "screen table")
  if (anyDuplicated(screen$trait)) stop("Screen table has duplicate traits")
  data.table::setnames(screen, "trait", "trait_id")
  out <- merge(out, screen[, .(trait_id, screen_level)], by = "trait_id", all.x = TRUE)
  if (anyNA(out$screen_level)) stop("Screen level is missing for at least one IV trait")
}

out[, `:=`(
  ld_panel_missing_n = NA_integer_,
  ld_panel_eligible_n = NA_integer_,
  ld_panel_missing_fraction = NA_real_,
  ld_panel_missingness_note = "Selection-stage LD-panel denominator was not supplied; final-IV files alone cannot reconstruct missingness"
)]
if (!is.null(selection_qc_file)) {
  selection <- data.table::fread(assert_file(selection_qc_file, "instrument-selection QC"))
  assert_columns(selection, c("trait", "autosomal_biallelic", "in_ld_panel", "ld_panel_missing"),
                 "instrument-selection QC")
  selection <- selection[, .(
    trait_id = trait, ld_panel_missing_n = as.integer(ld_panel_missing),
    ld_panel_eligible_n = as.integer(autosomal_biallelic),
    ld_panel_missing_fraction = ld_panel_missing / autosomal_biallelic
  )]
  out[, c("ld_panel_missing_n", "ld_panel_eligible_n", "ld_panel_missing_fraction") := NULL]
  out <- merge(out, selection, by = "trait_id", all.x = TRUE)
  out[!is.na(ld_panel_missing_n),
      ld_panel_missingness_note := "Computed from the instrument-selection QC table"]
}

data.table::setcolorder(out, c(
  "trait_id", "trait_display", intersect("screen_level", names(out)),
  "n_iv", "min_f", "median_f",
  "cumulative_r2_approx", "maf_min", "maf_q25", "maf_median", "maf_q75", "maf_max",
  "i2gx", "nome_i2gx_ge_0_90", "ld_panel_missing_n", "ld_panel_eligible_n",
  "ld_panel_missing_fraction", "ld_panel_missingness_note"
))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, out_file, sep = "\t", na = "NA", quote = FALSE)
