#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
stopifnot(
  window_overlaps_mhc(6L, 24000000, 25000000),
  window_overlaps_mhc(6L, 34000000, 35000000),
  !window_overlaps_mhc(6L, 24000000, 24999999),
  abs(i2gx_from_bse(c(-2, 0, 2), c(1, 1, 1)) - 0.75) < 1e-12
)
needed <- c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO", "coloc")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Smoke test requires: ", paste(missing, collapse = ", "))

root <- normalizePath(getwd(), winslash = "/")
rscript <- file.path(R.home("bin"), "Rscript")
work <- tempfile("mdd_mr_smoke_")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

# MR smoke: build one harmonised synthetic trait with the same TwoSampleMR path.
set.seed(20260815)
n <- 20L
snps <- paste0("rs", 100001L + seq_len(n))
alleles <- data.frame(SNP = snps, ea = rep(c("A", "C"), length.out = n),
                      oa = rep(c("G", "T"), length.out = n))
exp <- data.frame(alleles, beta = runif(n, 0.04, 0.12), se = runif(n, 0.008, 0.018),
                  eaf = runif(n, 0.10, 0.45), p = runif(n, 1e-12, 1e-9))
out <- data.frame(alleles, beta = 0.20 * exp$beta + rnorm(n, 0, 0.01),
                  se = runif(n, 0.012, 0.025), eaf = exp$eaf, p = runif(n, 0.01, 0.9))
ef <- TwoSampleMR::format_data(exp, "exposure", snp_col = "SNP", beta_col = "beta", se_col = "se",
  effect_allele_col = "ea", other_allele_col = "oa", eaf_col = "eaf", pval_col = "p")
of <- TwoSampleMR::format_data(out, "outcome", snp_col = "SNP", beta_col = "beta", se_col = "se",
  effect_allele_col = "ea", other_allele_col = "oa", eaf_col = "eaf", pval_col = "p")
h2 <- TwoSampleMR::harmonise_data(ef, of, action = 2)
h3 <- TwoSampleMR::harmonise_data(ef, of, action = 3)
h2$ncase <- 10000L; h2$ncontrol <- 20000L
h3$ncase <- 10000L; h3$ncontrol <- 20000L
presso_probe <- run_mr_suite(h2, "synthetic exposure", "synthetic outcome",
                             presso_nb = 500L, presso_seed = 20260815L,
                             run_presso = TRUE)
presso_summary <- summarise_suite(presso_probe, "Synthetic_trait", 500L)
stopifnot(is.finite(presso_summary$presso_raw_se),
          is.finite(presso_summary$presso_raw_ci_lo),
          is.finite(presso_summary$presso_raw_ci_hi),
          presso_summary$presso_se_status == "available")
harm_dir <- file.path(work, "harm")
dir.create(harm_dir)
saveRDS(h2, file.path(harm_dir, "harmonised_main_Synthetic_trait.rds"))
saveRDS(h2, file.path(harm_dir, "harmonised_noukbb_Synthetic_trait.rds"))
saveRDS(h3, file.path(harm_dir, "harmonised_main_a3_Synthetic_trait.rds"))
mr_out <- file.path(work, "mr")
code <- system2(rscript, c("analysis/scripts/03_run_forward_mr.R", "--harm-dir", shQuote(harm_dir),
  "--out-dir", shQuote(mr_out), "--traits", "Synthetic_trait", "--skip-presso"))
stopifnot(code == 0L)
mr <- data.table::fread(file.path(mr_out, "mr_main_summary.tsv"))
stopifnot(nrow(mr) == 1L, mr$n_iv == n, is.finite(mr$ivw_b), is.finite(mr$ivw_p))

# Robustness pilot: a subset is explicit and its MR-PRESSO filename reflects NbDistribution.
robust_out <- file.path(work, "robustness")
code <- system2(rscript, c(
  "analysis/scripts/05_run_robustness.R", "--harm-dir", shQuote(harm_dir),
  "--out-dir", shQuote(robust_out), "--traits", "Synthetic_trait", "--pilot",
  "--presso-nb", "500", "--presso-seed", "20260815"
))
stopifnot(code == 0L,
          file.exists(file.path(robust_out, "pilot_robustness_summary.tsv")),
          file.exists(file.path(robust_out, "pilot_presso_rerun_500.tsv")),
          !file.exists(file.path(robust_out, "presso_rerun_10000.tsv")))
formal_reject <- suppressWarnings(system2(rscript, c(
  "analysis/scripts/05_run_robustness.R", "--harm-dir", shQuote(harm_dir),
  "--out-dir", shQuote(file.path(work, "robustness-invalid")),
  "--traits", "Synthetic_trait", "--presso-nb", "10000"
), stdout = FALSE, stderr = FALSE))
stopifnot(formal_reject != 0L)
formal_harm <- file.path(work, "harm-formal-guards")
dir.create(formal_harm)
for (j in seq_len(15L)) {
  saveRDS(h2, file.path(formal_harm, paste0("harmonised_main_Formal_trait_", j, ".rds")))
}
bad_nb <- suppressWarnings(system2(rscript, c(
  "analysis/scripts/05_run_robustness.R", "--harm-dir", shQuote(formal_harm),
  "--out-dir", shQuote(file.path(work, "robustness-bad-nb")), "--presso-nb", "9999"
), stdout = FALSE, stderr = FALSE))
no_manifest <- suppressWarnings(system2(rscript, c(
  "analysis/scripts/05_run_robustness.R", "--harm-dir", shQuote(formal_harm),
  "--out-dir", shQuote(file.path(work, "robustness-no-manifest")), "--presso-nb", "10000"
), stdout = FALSE, stderr = FALSE))
incomplete_manifest <- file.path(work, "iv_manifest_incomplete.tsv")
data.table::fwrite(data.table::data.table(
  trait = "Formal_trait_1", iv_file = file.path(work, "not-read.tsv")
), incomplete_manifest, sep = "\t")
bad_manifest <- suppressWarnings(system2(rscript, c(
  "analysis/scripts/05_run_robustness.R", "--harm-dir", shQuote(formal_harm),
  "--iv-manifest", shQuote(incomplete_manifest),
  "--out-dir", shQuote(file.path(work, "robustness-bad-manifest")), "--presso-nb", "10000"
), stdout = FALSE, stderr = FALSE))
stopifnot(bad_nb != 0L, no_manifest != 0L, bad_manifest != 0L)

# FinnGen smoke: exact 13-column schema, complete matching and neutral support label.
iv_file <- file.path(work, "synthetic_ivs.tsv")
data.table::fwrite(data.table::data.table(
  rsid = snps, variant_id = paste(1L, seq_len(n) + 100000L, alleles$oa, alleles$ea, sep = "_"),
  effect_allele = alleles$ea, other_allele = alleles$oa,
  beta = exp$beta, standard_error = exp$se, effect_allele_frequency = exp$eaf,
  p = exp$p, F = (exp$beta / exp$se)^2
), iv_file, sep = "\t")
iv_manifest <- file.path(work, "iv_manifest.tsv")
data.table::fwrite(data.table::data.table(trait = "Synthetic_trait", iv_file = iv_file),
                   iv_manifest, sep = "\t")
fg_file <- file.path(work, "finngen.tsv")
data.table::fwrite(data.table::data.table(
  `#chrom` = 1L, pos = seq_len(n) + 100000L, ref = alleles$oa, alt = alleles$ea,
  rsids = snps, nearest_genes = "GENE", pval = out$p,
  mlogp = -log10(out$p), beta = out$beta, sebeta = out$se,
  af_alt = out$eaf, af_alt_cases = out$eaf, af_alt_controls = out$eaf
), fg_file, sep = "\t")
pgc_file <- file.path(work, "pgc.tsv")
data.table::fwrite(mr[, .(trait, ivw_b, ivw_se, ivw_p)], pgc_file, sep = "\t")
fg_out <- file.path(work, "finngen")
code <- system2(rscript, c(
  "analysis/scripts/06_run_finngen.R", "--finngen", shQuote(fg_file),
  "--iv-manifest", shQuote(iv_manifest), "--pgc-summary", shQuote(pgc_file),
  "--out-dir", shQuote(fg_out), "--traits", "Synthetic_trait"
))
stopifnot(code == 0L)
fg_qc <- data.table::fread(file.path(fg_out, "pilot_finngen_r13_admission_and_tool_loss.tsv"))
fg_result <- data.table::fread(file.path(fg_out, "pilot_finngen_r13_mr_results.tsv"))
stopifnot(fg_qc$n_iv == n, fg_qc$n_matched == n, fg_qc$n_not_found == 0L,
          fg_qc$n_harmonise_used == n,
          "pilot_nominal_direction_label" %in% names(fg_result),
          !"p_bh" %in% names(fg_result),
          !"cross_outcome_support_label" %in% names(fg_result),
          !file.exists(file.path(fg_out, "finngen_r13_mr_results.tsv")))

formal_traits <- paste0("Synthetic_trait_", seq_len(15L))
formal_manifest <- file.path(work, "iv_manifest_15.tsv")
data.table::fwrite(data.table::data.table(trait = formal_traits, iv_file = iv_file),
                   formal_manifest, sep = "\t")
formal_pgc <- file.path(work, "pgc_15.tsv")
data.table::fwrite(data.table::data.table(
  trait = formal_traits, ivw_b = mr$ivw_b[[1]], ivw_se = mr$ivw_se[[1]], ivw_p = mr$ivw_p[[1]]
), formal_pgc, sep = "\t")
formal_fg_out <- file.path(work, "finngen-formal")
code <- system2(rscript, c(
  "analysis/scripts/06_run_finngen.R", "--finngen", shQuote(fg_file),
  "--iv-manifest", shQuote(formal_manifest), "--pgc-summary", shQuote(formal_pgc),
  "--out-dir", shQuote(formal_fg_out), "--n-tests", "15"
))
stopifnot(code == 0L)
formal_fg <- data.table::fread(file.path(formal_fg_out, "finngen_r13_mr_results.tsv"))
stopifnot(nrow(formal_fg) == 15L, data.table::uniqueN(formal_fg$trait) == 15L,
          isTRUE(all.equal(formal_fg$p_bh, p.adjust(formal_fg$fg_p, "BH", n = 15L))),
          all(formal_fg$cross_outcome_support_label %in%
                c("FDR-supported", "nominally-supported", "no-nominal-support")))
formal_reject <- suppressWarnings(system2(rscript, c(
  "analysis/scripts/06_run_finngen.R", "--finngen", shQuote(fg_file),
  "--iv-manifest", shQuote(iv_manifest), "--pgc-summary", shQuote(pgc_file),
  "--out-dir", shQuote(file.path(work, "finngen-invalid")), "--n-tests", "15"
), stdout = FALSE, stderr = FALSE))
stopifnot(formal_reject != 0L)

# Colocalisation smoke: one 80-SNP locus and two case-control outcomes.
nloc <- 80L
rs <- paste0("rs", 200001L + seq_len(nloc))
ea <- rep(c("A", "C"), length.out = nloc)
oa <- rep(c("G", "T"), length.out = nloc)
bx <- rnorm(nloc, 0, 0.08); bx[1] <- 0.35
window <- data.table::data.table(rsid = rs, beta = bx, se = rep(0.03, nloc),
  eaf = seq(0.1, 0.45, length.out = nloc), effect_allele = ea, other_allele = oa)
window_file <- file.path(work, "window.tsv")
data.table::fwrite(window, window_file, sep = "\t")

make_outcome <- function(scale, path) {
  by <- scale * bx + rnorm(nloc, 0, 0.015)
  z <- data.table::data.table(SNP = rs, A1 = ea, A2 = oa, OR = exp(by), SE = 0.025,
    FRQ_A_10000 = window$eaf, FRQ_U_20000 = window$eaf, Nca = 10000L, Nco = 20000L)
  data.table::fwrite(z, path, sep = "\t")
}
main_file <- file.path(work, "main.tsv")
noukbb_file <- file.path(work, "noukbb.tsv")
make_outcome(0.25, main_file); make_outcome(0.23, noukbb_file)
locus_manifest <- file.path(work, "loci.tsv")
data.table::fwrite(data.table::data.table(
  locus_key = "Synthetic_trait__locus001", trait = "Synthetic_trait", locus_id = "locus001",
  chr = 1L, lead_rsid = rs[[1]], lead_pos = 1000000L,
  window_start = 500000L, window_end = 1500000L, window_file = window_file
), locus_manifest, sep = "\t")
outcome_manifest <- file.path(work, "outcomes.tsv")
data.table::fwrite(data.table::data.table(
  outcome = c("main", "noUKBB"), file = c(main_file, noukbb_file),
  ncase = 10000L, ncontrol = 20000L, frq_case_col = "FRQ_A_10000",
  frq_control_col = "FRQ_U_20000"
), outcome_manifest, sep = "\t")
coloc_out <- file.path(work, "coloc")
code <- system2(rscript, c("analysis/scripts/07_run_coloc_abf.R",
  "--locus-manifest", shQuote(locus_manifest), "--outcome-manifest", shQuote(outcome_manifest),
  "--out-dir", shQuote(coloc_out), "--min-snps", "50"))
stopifnot(code == 0L)
cr <- data.table::fread(file.path(coloc_out, "loci_coloc_results.tsv"))
cs <- data.table::fread(file.path(coloc_out, "loci_status.tsv"))
if (!nrow(cr)) print(cs)
stopifnot(nrow(cr) == 8L, nrow(cs) == 2L, all(cs$analysis_status == "success"),
          all(abs(cr$PP.H0 + cr$PP.H1 + cr$PP.H2 + cr$PP.H3 + cr$PP.H4 - 1) < 0.01))

code <- system2(rscript, c("analysis/scripts/08_classify_coloc.R",
  "--locus-manifest", shQuote(locus_manifest),
  "--results", shQuote(file.path(coloc_out, "loci_coloc_results.tsv")),
  "--status", shQuote(file.path(coloc_out, "loci_status.tsv")),
  "--final-snps", shQuote(file.path(coloc_out, "loci_final_snps.tsv")),
  "--exposure-lead-qc", shQuote(file.path(coloc_out, "exposure_lead_qc.tsv")),
  "--out-dir", shQuote(coloc_out)))
stopifnot(code == 0L)
cc <- data.table::fread(file.path(coloc_out, "loci_classification.tsv"))
stopifnot(nrow(cc) == 1L, cc$locus_key == "Synthetic_trait__locus001",
          cc$abf_class %in% c("robust_coloc", "prior_sensitive_coloc", "distinct_signal",
                              "trait_specific_or_low_power", "inconclusive"))

cat("PASS: synthetic MR and colocalisation smoke tests\n")
