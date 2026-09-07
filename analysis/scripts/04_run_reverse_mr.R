#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO"))

assoc_file <- assert_file(arg_required(args, "assoc"), "reverse association table")
out_dir <- arg_required(args, "out-dir")
n_tests <- arg_value(args, "n-tests", 249L, "integer")
presso_nb <- arg_value(args, "presso-nb", 1000L, "integer")
presso_seed <- arg_value(args, "presso-seed", 20260815L, "integer")
run_presso <- !isTRUE(args[["skip-presso"]])
traits_keep <- csv_values(args[["traits"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dall <- data.table::fread(assoc_file)
req <- c("rsid", "mdd_A1", "mdd_A2", "mdd_beta", "mdd_se", "mdd_p", "mdd_eaf",
         "trait", "met_effect", "met_other", "met_beta", "met_se", "met_eaf", "met_p")
assert_columns(dall, req, "reverse association table")
traits <- sort(unique(dall$trait))
if (length(traits_keep)) traits <- intersect(traits, traits_keep)

summary_rows <- diagnostic_rows <- long_rows <- loo_rows <- a3_rows <- qc_rows <- list()
for (i in seq_along(traits)) {
  tr <- traits[[i]]
  d <- dall[trait == tr]
  exp <- TwoSampleMR::format_data(as.data.frame(d), type = "exposure", snp_col = "rsid",
    beta_col = "mdd_beta", se_col = "mdd_se", effect_allele_col = "mdd_A1",
    other_allele_col = "mdd_A2", eaf_col = "mdd_eaf", pval_col = "mdd_p")
  out <- TwoSampleMR::format_data(as.data.frame(d), type = "outcome", snp_col = "rsid",
    beta_col = "met_beta", se_col = "met_se", effect_allele_col = "met_effect",
    other_allele_col = "met_other", eaf_col = "met_eaf", pval_col = "met_p")
  h2 <- TwoSampleMR::harmonise_data(exp, out, action = 2)
  h3 <- TwoSampleMR::harmonise_data(exp, out, action = 3)
  s <- run_mr_suite(h2, "major depressive disorder", tr, presso_nb, presso_seed, run_presso)
  row <- summarise_suite(s, tr, presso_nb)
  row[, `:=`(ivw_or = NULL, ivw_or_lci = NULL, ivw_or_uci = NULL)]
  summary_rows[[i]] <- row
  diagnostic_rows[[i]] <- suite_diagnostics(s, tr)
  if (!is.null(s$mr_long)) { z <- data.table::as.data.table(s$mr_long); z[, trait := tr]; long_rows[[i]] <- z }
  if (!is.null(s$loo)) { z <- data.table::as.data.table(s$loo); z[, trait := tr]; loo_rows[[i]] <- z }
  s3 <- run_mr_suite(h3, "major depressive disorder", tr, 0L, presso_seed, FALSE)
  a3_rows[[i]] <- summarise_suite(s3, tr, 0L)[, .(trait, n_iv, ivw_b, ivw_se, ivw_p)]
  qc_rows[[i]] <- data.table::data.table(
    trait = tr, n_assoc = nrow(d), action2_keep = sum(h2$mr_keep %in% TRUE),
    action3_keep = sum(h3$mr_keep %in% TRUE),
    action2_palindromic = sum(h2$palindromic %in% TRUE)
  )
}

main <- data.table::rbindlist(summary_rows, fill = TRUE)
if (!length(traits_keep)) {
  main <- classify_full_screen(main, "ivw_p", n_tests)
  data.table::fwrite(main, file.path(out_dir, "reverse_mr_summary_screen.tsv"), sep = "\t")
} else {
  main[, `:=`(
    p_bonf = p.adjust(ivw_p, method = "bonferroni", n = n_tests),
    pilot_fdr = p.adjust(ivw_p, method = "BH")
  )]
  main[, sig_level := ifelse(p_bonf < 0.05, "bonferroni_hit", "not_bonferroni_significant")]
}
data.table::fwrite(main, file.path(out_dir, "reverse_mr_summary.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(long_rows, fill = TRUE), file.path(out_dir, "reverse_mr_long.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(diagnostic_rows, fill = TRUE), file.path(out_dir, "reverse_mr_sensitivity.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(loo_rows, fill = TRUE), file.path(out_dir, "reverse_leaveoneout.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(a3_rows, fill = TRUE), file.path(out_dir, "reverse_mr_a3_summary.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(qc_rows), file.path(out_dir, "reverse_harmonisation_qc.tsv"), sep = "\t")
write_run_metadata(out_dir, "reverse_mr", list(association_table = assoc_file), list(
  n_tests = n_tests, presso_nb = presso_nb, presso_seed = presso_seed,
  harmonise_primary_action = 2L, harmonise_sensitivity_action = 3L
))
