#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR", "MendelianRandomization"))

finngen_file <- assert_file(arg_required(args, "finngen"), "FinnGen summary statistics")
iv_manifest_file <- assert_file(arg_required(args, "iv-manifest"), "IV manifest")
pgc_summary_file <- assert_file(arg_required(args, "pgc-summary"), "PGC MR summary")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
n_tests <- arg_value(args, "n-tests", 15L, "integer")
pilot <- length(traits_keep) > 0L || isTRUE(args[["pilot"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ivm <- data.table::fread(iv_manifest_file)
assert_columns(ivm, c("trait", "iv_file"), "IV manifest")
if (anyDuplicated(ivm$trait) || any(!nzchar(ivm$trait))) {
  stop("IV manifest must contain unique, non-empty trait identifiers")
}
if (length(traits_keep)) {
  missing_traits <- setdiff(traits_keep, ivm$trait)
  if (length(missing_traits)) stop("Requested trait absent from IV manifest: ", missing_traits[[1]])
  ivm <- ivm[trait %in% traits_keep]
}
if (!nrow(ivm)) stop("No FinnGen traits selected")
if (!pilot && (n_tests != 15L || nrow(ivm) != 15L || data.table::uniqueN(ivm$trait) != 15L)) {
  stop("Formal FinnGen mode requires exactly 15 unique candidate traits and --n-tests 15; use --traits for a pilot subset")
}
pgc <- data.table::fread(pgc_summary_file)
assert_columns(pgc, c("trait", "ivw_b", "ivw_se", "ivw_p"), "PGC MR summary")
pgc_selected <- pgc[trait %in% ivm$trait]
if (anyDuplicated(pgc_selected$trait) || nrow(pgc_selected) != nrow(ivm) ||
    !setequal(pgc_selected$trait, ivm$trait)) {
  stop("PGC summary must contain exactly one row for every selected FinnGen trait")
}
assert_numeric(pgc_selected, c("ivw_b", "ivw_se", "ivw_p"), "selected PGC summary")
if (any(pgc_selected$ivw_se <= 0) || any(pgc_selected$ivw_p < 0 | pgc_selected$ivw_p > 1)) {
  stop("Selected PGC summary failed SE/P range checks")
}

header <- data.table::fread(finngen_file, nrows = 0L)
expected_header <- c("#chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
                     "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls")
if (!identical(names(header), expected_header)) stop("FinnGen header does not match the locked 13-column R13 schema")
required_fg <- c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt")
assert_columns(header, required_fg, "FinnGen file")
message("Reading the FinnGen file once")
fg <- data.table::fread(finngen_file, select = required_fg, showProgress = TRUE)
assert_numeric(fg, c("pos", "pval", "beta", "sebeta", "af_alt"), "FinnGen file")
if (any(fg$sebeta <= 0) || any(fg$pval < 0 | fg$pval > 1) || any(fg$af_alt < 0 | fg$af_alt > 1)) {
  stop("FinnGen numeric-range admission check failed")
}
fg[, poskey := paste(`#chrom`, pos, sep = ":")]

cross_het <- function(b1, se1, b2, se2) {
  w1 <- 1 / se1^2; w2 <- 1 / se2^2
  bbar <- (w1 * b1 + w2 * b2) / (w1 + w2)
  q <- w1 * (b1 - bbar)^2 + w2 * (b2 - bbar)^2
  c(Q = q, p = pchisq(q, 1, lower.tail = FALSE), I2 = if (q > 1) (q - 1) / q * 100 else 0)
}

iv_rows <- vector("list", nrow(ivm))
for (i in seq_len(nrow(ivm))) {
  tr <- ivm$trait[[i]]
  iv <- data.table::fread(assert_file(ivm$iv_file[[i]], paste0("IV set for ", tr)))
  assert_columns(iv, c("rsid", "variant_id", "effect_allele", "other_allele", "beta",
                       "standard_error", "effect_allele_frequency", "p"), tr)
  iv[, c("chr", "pos") := data.table::tstrsplit(variant_id, "_", fixed = TRUE)[1:2]]
  iv[, `:=`(trait = tr, chr = as.integer(chr), pos = as.integer(pos))]
  if (anyNA(iv$chr) || anyNA(iv$pos) || anyDuplicated(iv$rsid)) {
    stop(tr, ": invalid coordinates or duplicate rsIDs in the IV set")
  }
  iv[, poskey := paste(chr, pos, sep = ":")]
  iv_rows[[i]] <- iv[, .(
    trait, rsid, poskey, effect_allele, other_allele,
    effect_allele_frequency, beta_exp = beta, se_exp = standard_error, p_exp = p
  )]
}
iv_all <- data.table::rbindlist(iv_rows)

cand <- fg[poskey %in% unique(iv_all$poskey)]
cand[, fg_rsid := strsplit(rsids, ",", fixed = TRUE)]
cand <- cand[, .(
  rsid = unlist(fg_rsid)
), by = .(poskey, ref, alt, pval, beta_fg = beta, se_fg = sebeta, af_alt)]
matched <- merge(cand, iv_all, by = c("poskey", "rsid"))
matched[, allele_ok := paste(pmin(ref, alt), pmax(ref, alt)) ==
                           paste(pmin(effect_allele, other_allele),
                                 pmax(effect_allele, other_allele))]
per_iv <- matched[, .(
  n_rows = .N, n_ok = sum(allele_ok), n_bad = sum(!allele_ok)
), by = .(trait, rsid)]
per_iv[, match_status := data.table::fcase(
  n_rows > 1L, "duplicate_conflict",
  n_rows == 1L & n_bad == 1L, "allele_conflict",
  n_rows == 1L & n_ok == 1L, "matched",
  default = "not_found"
)]
iv_status <- merge(unique(iv_all[, .(trait, rsid)]), per_iv,
                   by = c("trait", "rsid"), all.x = TRUE)
iv_status[is.na(match_status), match_status := "not_found"]
if (any(iv_status$match_status %in% c("allele_conflict", "duplicate_conflict"))) {
  stop("FinnGen matching found allele or duplicate conflicts; fail closed")
}
matched_ok <- merge(
  matched[allele_ok == TRUE],
  iv_status[match_status == "matched", .(trait, rsid)],
  by = c("trait", "rsid")
)
if (anyDuplicated(matched_ok[, .(trait, rsid)])) stop("Duplicate FinnGen match; fail closed")

tool_loss <- iv_status[, .(
  n_iv = .N,
  n_matched = sum(match_status == "matched"),
  n_allele_conflict = sum(match_status == "allele_conflict"),
  n_duplicate_conflict = sum(match_status == "duplicate_conflict"),
  n_not_found = sum(match_status == "not_found")
), by = trait]

rows <- methods <- harmonisation_qc <- trait_status <- list()
for (i in seq_len(nrow(ivm))) {
  tr <- ivm$trait[[i]]
  m <- matched_ok[trait == tr]
  if (!nrow(m)) {
    harmonisation_qc[[i]] <- data.table::data.table(
      trait = tr, n_harmonise_used = 0L, palindromic = NA_integer_,
      flipped = NA_integer_, removed = NA_integer_
    )
    trait_status[[i]] <- data.table::data.table(
      trait = tr, analysis_status = "incomplete", failure_reason = "no_matched_instruments"
    )
    next
  }
  exp <- TwoSampleMR::format_data(as.data.frame(m), type = "exposure", snp_col = "rsid",
    beta_col = "beta_exp", se_col = "se_exp", effect_allele_col = "effect_allele",
    other_allele_col = "other_allele", eaf_col = "effect_allele_frequency", pval_col = "p_exp")
  out <- TwoSampleMR::format_data(as.data.frame(m), type = "outcome", snp_col = "rsid",
    beta_col = "beta_fg", se_col = "se_fg", effect_allele_col = "alt",
    other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
  h_all <- TwoSampleMR::harmonise_data(exp, out, action = 2)
  h <- h_all[h_all$mr_keep %in% TRUE, , drop = FALSE]
  log <- attr(h_all, "log")
  flipped <- if (is.null(log)) NA_integer_ else
    sum(log$switched_alleles, log$flipped_alleles_basic,
        log$flipped_alleles_palindrome, na.rm = TRUE)
  harmonisation_qc[[i]] <- data.table::data.table(
    trait = tr, n_harmonise_used = nrow(h),
    palindromic = sum(h_all$palindromic %in% TRUE, na.rm = TRUE),
    flipped = as.integer(flipped), removed = nrow(m) - nrow(h)
  )
  if (nrow(h) < 2L) {
    trait_status[[i]] <- data.table::data.table(
      trait = tr, analysis_status = "incomplete",
      failure_reason = "fewer_than_two_harmonised_instruments"
    )
    next
  }
  s <- run_mr_suite(h, tr, "FinnGen R13 register-based depression", 0L, 20260815L, FALSE)
  rr <- summarise_suite(s, tr, 0L)
  prow <- pgc[trait == tr]
  if (nrow(prow) != 1L) stop("PGC summary must contain exactly one row for ", tr)
  if (any(!is.finite(c(rr$ivw_b, rr$ivw_se, rr$ivw_p))) || rr$ivw_se <= 0 ||
      rr$ivw_p < 0 || rr$ivw_p > 1) {
    trait_status[[i]] <- data.table::data.table(
      trait = tr, analysis_status = "incomplete", failure_reason = "ivw_estimate_unavailable"
    )
    next
  }
  ht <- cross_het(rr$ivw_b, rr$ivw_se, prow$ivw_b[[1]], prow$ivw_se[[1]])
  ml <- data.table::as.data.table(s$mr_long)
  if (nrow(ml)) { ml[, trait := tr]; methods[[i]] <- ml }
  result_row <- rr[, .(
    trait, n_iv, fg_b = ivw_b, fg_se = ivw_se, fg_p = ivw_p,
    fg_or = ivw_or, fg_or_lci = ivw_or_lci, fg_or_uci = ivw_or_uci
  )]
  result_row[, `:=`(
    pgc_b = prow$ivw_b[[1]], pgc_se = prow$ivw_se[[1]], pgc_p = prow$ivw_p[[1]]
  )]
  result_row[, `:=`(
    dir_same = sign(fg_b) == sign(pgc_b),
    overlap_uncorrected_Q = ht[["Q"]], overlap_uncorrected_Q_p = ht[["p"]],
    overlap_uncorrected_I2_pct = ht[["I2"]]
  )]
  rows[[i]] <- result_row
  trait_status[[i]] <- data.table::data.table(
    trait = tr, analysis_status = "success", failure_reason = "none"
  )
}

tool_loss <- merge(tool_loss, data.table::rbindlist(harmonisation_qc, fill = TRUE),
                   by = "trait", all.x = TRUE)
status <- data.table::rbindlist(trait_status, fill = TRUE)
tool_loss <- merge(tool_loss, status, by = "trait", all.x = TRUE)
res <- data.table::rbindlist(rows, fill = TRUE)
prefix <- if (pilot) "pilot_" else ""
data.table::fwrite(data.table::rbindlist(methods, fill = TRUE),
                   file.path(out_dir, paste0(prefix, "finngen_r13_mr_methods.tsv")), sep = "\t")
data.table::fwrite(tool_loss,
                   file.path(out_dir, paste0(prefix, "finngen_r13_admission_and_tool_loss.tsv")), sep = "\t")
data.table::fwrite(status,
                   file.path(out_dir, paste0(prefix, "finngen_r13_trait_status.tsv")), sep = "\t")
admission <- data.table::data.table(
  check = c("header_13cols", "se_positive", "p_in_range", "af_in_range", "finite_numeric"),
  value = c("PASS", "PASS", "PASS", "PASS", "PASS")
)
data.table::fwrite(admission, file.path(out_dir, paste0(prefix, "finngen_r13_admission_qc.tsv")), sep = "\t")
if (!pilot && (nrow(res) != 15L || nrow(status) != 15L || any(status$analysis_status != "success"))) {
  stop("Formal FinnGen classification aborted: all 15 traits require successful MR; inspect trait-status and tool-loss outputs")
}
if (pilot) {
  if (nrow(res)) {
    res[, pilot_nominal_direction_label := data.table::fcase(
      dir_same & fg_p < 0.05, "directionally-concordant-nominal",
      default = "no-directionally-concordant-nominal-signal"
    )]
  }
} else {
  res[, `:=`(
    p_bh = p.adjust(fg_p, "BH", n = 15L),
    p_bonf = p.adjust(fg_p, "bonferroni", n = 15L)
  )]
  res[, cross_outcome_support_label := data.table::fcase(
    dir_same & p_bh < 0.05, "FDR-supported",
    dir_same & fg_p < 0.05, "nominally-supported",
    default = "no-nominal-support"
  )]
}
data.table::fwrite(res, file.path(out_dir, paste0(prefix, "finngen_r13_mr_results.tsv")), sep = "\t")
write_run_metadata(out_dir, "finngen", list(finngen = finngen_file, iv_manifest = iv_manifest_file,
                                             pgc_summary = pgc_summary_file),
                   list(harmonise_action = 2L,
                        analysis_mode = if (pilot) "pilot subset" else "formal 15-trait family",
                        n_tests = if (pilot) "not applied" else 15L,
                        cross_outcome_support_rule = if (pilot) {
                          "pilot nominal-direction label only; no formal BH or cross-outcome support label"
                        } else {
                          "direction concordance plus FinnGen BH or nominal P; Q excluded"
                        }))
