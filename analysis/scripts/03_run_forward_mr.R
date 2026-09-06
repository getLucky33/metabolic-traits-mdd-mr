#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
source(file.path("analysis", "R", "mr_engine.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR", "MendelianRandomization", "MRPRESSO"))

harm_dir <- arg_required(args, "harm-dir")
out_dir <- arg_required(args, "out-dir")
n_tests <- arg_value(args, "n-tests", 249L, "integer")
presso_nb <- arg_value(args, "presso-nb", 1000L, "integer")
presso_seed <- arg_value(args, "presso-seed", 20260815L, "integer")
run_presso <- !isTRUE(args[["skip-presso"]])
traits_keep <- csv_values(args[["traits"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(harm_dir, pattern = "^harmonised_main_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("^harmonised_main_a3_", basename(files))]
traits <- sub("^harmonised_main_(.*)\\.rds$", "\\1", basename(files))
if (length(traits_keep)) {
  keep <- traits %in% safe_trait(traits_keep)
  files <- files[keep]; traits <- traits[keep]
}
if (!length(files)) stop("No harmonised_main_*.rds inputs")

summary_rows <- diagnostic_rows <- long_rows <- loo_rows <- list()
noukbb_rows <- noukbb_long_rows <- a3_rows <- list()
for (i in seq_along(files)) {
  tr <- traits[[i]]
  message("Forward MR: ", tr)
  s <- run_mr_suite(readRDS(files[[i]]), tr, "major depressive disorder",
                    presso_nb, presso_seed, run_presso)
  summary_rows[[i]] <- summarise_suite(s, tr, presso_nb)
  diagnostic_rows[[i]] <- suite_diagnostics(s, tr)
  if (!is.null(s$mr_long)) { z <- data.table::as.data.table(s$mr_long); z[, trait := tr]; long_rows[[i]] <- z }
  if (!is.null(s$loo)) { z <- data.table::as.data.table(s$loo); z[, trait := tr]; loo_rows[[i]] <- z }

  fn <- file.path(harm_dir, paste0("harmonised_noukbb_", tr, ".rds"))
  if (file.exists(fn)) {
    sn <- run_mr_suite(readRDS(fn), tr, "major depressive disorder (UK Biobank excluded)",
                       presso_nb, presso_seed, run_presso)
    noukbb_rows[[i]] <- summarise_suite(sn, tr, presso_nb)
    if (!is.null(sn$mr_long)) { z <- data.table::as.data.table(sn$mr_long); z[, trait := tr]; noukbb_long_rows[[i]] <- z }
  }
  fa3 <- file.path(harm_dir, paste0("harmonised_main_a3_", tr, ".rds"))
  if (file.exists(fa3)) {
    sa3 <- run_mr_suite(readRDS(fa3), tr, "major depressive disorder", 0L, presso_seed, FALSE)
    a3_rows[[i]] <- summarise_suite(sa3, tr, 0L)[, .(trait, n_iv, ivw_b, ivw_se, ivw_p)]
  }
}

main <- data.table::rbindlist(summary_rows, fill = TRUE)
main[, `:=`(p_bonf = p.adjust(ivw_p, method = "bonferroni", n = n_tests),
            pilot_fdr = p.adjust(ivw_p, method = "BH"))]
main[, sig_level := ifelse(p_bonf < 0.05, "bonferroni_hit", "not_bonferroni_significant")]
data.table::fwrite(main, file.path(out_dir, "mr_main_summary.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(long_rows, fill = TRUE), file.path(out_dir, "mr_main_long.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(diagnostic_rows, fill = TRUE), file.path(out_dir, "mr_sensitivity.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(loo_rows, fill = TRUE), file.path(out_dir, "mr_leaveoneout.tsv"), sep = "\t")
if (length(noukbb_rows)) data.table::fwrite(data.table::rbindlist(noukbb_rows, fill = TRUE), file.path(out_dir, "mr_noukbb_summary.tsv"), sep = "\t")
if (length(noukbb_long_rows)) data.table::fwrite(data.table::rbindlist(noukbb_long_rows, fill = TRUE), file.path(out_dir, "mr_noukbb_long.tsv"), sep = "\t")
if (length(a3_rows)) data.table::fwrite(data.table::rbindlist(a3_rows, fill = TRUE), file.path(out_dir, "mr_a3_summary.tsv"), sep = "\t")
write_run_metadata(out_dir, "forward_mr", setNames(as.list(files), paste0("harmonised_", traits)), list(
  n_tests = n_tests, presso_nb = presso_nb, presso_seed = presso_seed,
  harmonise_primary_action = 2L, harmonise_sensitivity_action = 3L
))
