#!/usr/bin/env Rscript
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
harm_dir <- file.path(work, "harm")
dir.create(harm_dir)
saveRDS(h2, file.path(harm_dir, "harmonised_main_Synthetic_trait.rds"))
saveRDS(h2, file.path(harm_dir, "harmonised_noukbb_Synthetic_trait.rds"))
saveRDS(h3, file.path(harm_dir, "harmonised_main_a3_Synthetic_trait.rds"))
mr_out <- file.path(work, "mr")
code <- system2(rscript, c("analysis/scripts/03_run_forward_mr.R", "--harm-dir", shQuote(harm_dir),
  "--out-dir", shQuote(mr_out), "--skip-presso"))
stopifnot(code == 0L)
mr <- data.table::fread(file.path(mr_out, "mr_main_summary.tsv"))
stopifnot(nrow(mr) == 1L, mr$n_iv == n, is.finite(mr$ivw_b), is.finite(mr$ivw_p))

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
  chr = 1L, lead_rsid = rs[[1]], lead_pos = 1000000L, window_file = window_file
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
