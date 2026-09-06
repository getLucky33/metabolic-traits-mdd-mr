#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages("data.table")

manifest_file <- assert_file(arg_required(args, "locus-manifest"), "locus manifest")
results_file <- assert_file(arg_required(args, "results"), "coloc results")
status_file <- assert_file(arg_required(args, "status"), "coloc status")
final_snps_file <- assert_file(arg_required(args, "final-snps"), "final SNP table")
exposure_qc_file <- assert_file(arg_required(args, "exposure-lead-qc"), "exposure lead QC")
out_dir <- arg_required(args, "out-dir")
plink_bin <- args[["plink"]]
ld_bfile <- args[["ld-bfile"]]
if (!is.null(plink_bin)) plink_bin <- assert_file(plink_bin, "PLINK executable")
if (!is.null(ld_bfile) && (!file.exists(paste0(ld_bfile, ".bed")) ||
    !file.exists(paste0(ld_bfile, ".bim")) || !file.exists(paste0(ld_bfile, ".fam")))) {
  stop("LD reference requires .bed/.bim/.fam: ", ld_bfile)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

classify_full <- function(st_main, st_noukbb, pp_valid_main, pp_valid_noukbb,
                          qc_exp, qc_main, qc_noukbb, is_mhc,
                          h1m, h2m, h3m, h4m, h3n, h4n,
                          h4_1e6_m, h4_1e6_n, h3_1e5_m, h4_1e5_m) {
  if (st_main != "success") return(c("inconclusive", st_main))
  if (st_noukbb != "success") return(c("inconclusive", st_noukbb))
  if (!pp_valid_main || !pp_valid_noukbb) return(c("inconclusive", "posterior_missing_invalid"))
  if (!qc_exp || !qc_main || !qc_noukbb) return(c("inconclusive", "lead_missing"))
  if (is_mhc) return(c("inconclusive", "mhc_excluded"))
  ratio_pass <- function(h4, h3) h4 > 0.90 && ((h3 == 0 && h4 > 0) || (h3 > 0 && h4 / h3 > 3))
  if (ratio_pass(h4m, h3m) && ratio_pass(h4n, h3n) && h4_1e6_m > 0.80 && h4_1e6_n > 0.80)
    return(c("robust_coloc", "none"))
  if (ratio_pass(h4m, h3m) || ratio_pass(h4_1e5_m, h3_1e5_m))
    return(c("prior_sensitive_coloc", "none"))
  if (h3m > 0.80 && h4m < 0.10 && h3n > 0.80 && h4n < 0.10)
    return(c("distinct_signal", "none"))
  if (h1m + h2m > 0.80 && h4m < 0.10)
    return(c("trait_specific_or_low_power", "none"))
  c("inconclusive", "posterior_indeterminate")
}

man <- data.table::fread(manifest_file)
res <- data.table::fread(results_file)
sts <- data.table::fread(status_file)
fs <- data.table::fread(final_snps_file)
eqc <- data.table::fread(exposure_qc_file)
assert_columns(man, c("locus_key", "trait", "locus_id", "chr", "lead_rsid", "lead_pos"), "locus manifest")
assert_columns(res, c("locus_key", "outcome", "p12", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"), "results")

res[, pp_sum := PP.H0 + PP.H1 + PP.H2 + PP.H3 + PP.H4]
valid <- res[, .(pp_valid = .N == 4L && data.table::uniqueN(p12) == 4L &&
                   all(is.finite(c(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4))) &&
                   all(c(PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) >= 0) &&
                   all(pp_sum >= 0.99 & pp_sum <= 1.01)), by = .(locus_key, outcome)]

wide_st <- data.table::dcast(sts, locus_key ~ outcome, value.var = "analysis_status")
wide_valid <- data.table::dcast(valid, locus_key ~ outcome, value.var = "pp_valid", fill = FALSE)

resolve_lead <- function(key, outcome_name, lead_rsid) {
  final <- fs[locus_key == key & outcome == outcome_name, unique(rsid)]
  if (lead_rsid %in% final) return(data.table::data.table(
    locus_key = key, outcome = outcome_name, lead_present = TRUE,
    lead_resolution = "direct", proxy_rsid = NA_character_, proxy_r2 = NA_real_))
  if (is.null(plink_bin) || is.null(ld_bfile)) {
    stop("Lead ", lead_rsid, " is absent from ", key, "/", outcome_name,
         "; provide --plink and --ld-bfile for the frozen r2>0.8 proxy search")
  }
  prefix <- tempfile(pattern = "plink_proxy_", tmpdir = tempdir())
  on.exit(unlink(paste0(prefix, c(".ld", ".log", ".nosex"))), add = TRUE)
  code <- system2(plink_bin, c(
    "--bfile", shQuote(ld_bfile), "--r2", "--ld-snp", lead_rsid,
    "--ld-window-kb", "2000", "--ld-window", "999999", "--ld-window-r2", "0.8",
    "--out", shQuote(prefix)
  ), stdout = FALSE, stderr = FALSE)
  ld_file <- paste0(prefix, ".ld")
  if (code != 0L || !file.exists(ld_file)) return(data.table::data.table(
    locus_key = key, outcome = outcome_name, lead_present = FALSE,
    lead_resolution = "plink_failed_or_lead_absent", proxy_rsid = NA_character_, proxy_r2 = NA_real_))
  ld <- data.table::fread(ld_file)
  if (!nrow(ld)) return(data.table::data.table(
    locus_key = key, outcome = outcome_name, lead_present = FALSE,
    lead_resolution = "no_proxy_r2_gt_0.8", proxy_rsid = NA_character_, proxy_r2 = NA_real_))
  candidates <- unique(data.table::rbindlist(list(
    ld[SNP_A == lead_rsid, .(proxy_rsid = SNP_B, proxy_r2 = R2)],
    ld[SNP_B == lead_rsid, .(proxy_rsid = SNP_A, proxy_r2 = R2)]
  )))[proxy_rsid %in% final][order(-proxy_r2)]
  if (!nrow(candidates)) return(data.table::data.table(
    locus_key = key, outcome = outcome_name, lead_present = FALSE,
    lead_resolution = "proxy_not_in_final_set", proxy_rsid = NA_character_, proxy_r2 = NA_real_))
  data.table::data.table(locus_key = key, outcome = outcome_name, lead_present = TRUE,
                         lead_resolution = "proxy", proxy_rsid = candidates$proxy_rsid[[1]],
                         proxy_r2 = candidates$proxy_r2[[1]])
}

lead_rows <- list()
for (i in seq_len(nrow(man))) {
  for (outcome_name in c("main", "noUKBB")) {
    lead_rows[[length(lead_rows) + 1L]] <- resolve_lead(
      man$locus_key[[i]], outcome_name, man$lead_rsid[[i]])
  }
}
lead <- data.table::rbindlist(lead_rows)
data.table::fwrite(lead, file.path(out_dir, "lead_proxy_qc.tsv"), sep = "\t")
wide_lead <- data.table::dcast(lead, locus_key ~ outcome, value.var = "lead_present", fill = FALSE)

x <- merge(man, eqc[, .(locus_key, qc_exp_lead_present)], by = "locus_key", all.x = TRUE)
x <- merge(x, wide_st, by = "locus_key", all.x = TRUE, suffixes = c("", "_status"))
data.table::setnames(x, intersect(c("main", "noUKBB"), names(x)), paste0("status_", intersect(c("main", "noUKBB"), names(x))))
x <- merge(x, wide_valid, by = "locus_key", all.x = TRUE, suffixes = c("", "_valid"))
if ("main" %in% names(x)) data.table::setnames(x, "main", "pp_valid_main")
if ("noUKBB" %in% names(x)) data.table::setnames(x, "noUKBB", "pp_valid_noUKBB")
x <- merge(x, wide_lead, by = "locus_key", all.x = TRUE, suffixes = c("", "_lead"))
if ("main" %in% names(x)) data.table::setnames(x, "main", "lead_main")
if ("noUKBB" %in% names(x)) data.table::setnames(x, "noUKBB", "lead_noUKBB")
x[, qc_exp_lead_present := qc_exp_lead_present | lead_main | lead_noUKBB]
get_pp <- function(key, outcome_name, prior, column) {
  z <- res[locus_key == key & outcome == outcome_name & abs(p12 - prior) < prior * 1e-8]
  if (nrow(z) != 1L) return(NA_real_)
  as.numeric(z[[column]][[1]])
}

rows <- vector("list", nrow(x))
for (i in seq_len(nrow(x))) {
  r <- x[i]
  st_m <- if ("status_main" %in% names(r) && !is.na(r$status_main)) r$status_main else "missing"
  st_n <- if ("status_noUKBB" %in% names(r) && !is.na(r$status_noUKBB)) r$status_noUKBB else "missing"
  z <- classify_full(
    st_m, st_n,
    isTRUE(r$pp_valid_main), isTRUE(r$pp_valid_noUKBB),
    isTRUE(r$qc_exp_lead_present), isTRUE(r$lead_main), isTRUE(r$lead_noUKBB),
    r$chr == 6L && r$lead_pos >= 25e6 && r$lead_pos <= 34e6,
    get_pp(r$locus_key, "main", 5e-6, "PP.H1"), get_pp(r$locus_key, "main", 5e-6, "PP.H2"),
    get_pp(r$locus_key, "main", 5e-6, "PP.H3"), get_pp(r$locus_key, "main", 5e-6, "PP.H4"),
    get_pp(r$locus_key, "noUKBB", 5e-6, "PP.H3"), get_pp(r$locus_key, "noUKBB", 5e-6, "PP.H4"),
    get_pp(r$locus_key, "main", 1e-6, "PP.H4"), get_pp(r$locus_key, "noUKBB", 1e-6, "PP.H4"),
    get_pp(r$locus_key, "main", 1e-5, "PP.H3"), get_pp(r$locus_key, "main", 1e-5, "PP.H4")
  )
  rows[[i]] <- data.table::data.table(
    locus_key = r$locus_key, trait = r$trait, locus = r$locus_id,
    abf_class = z[[1]], failure_reason = z[[2]],
    qc_exp_lead_present = isTRUE(r$qc_exp_lead_present),
    qc_mdd_lead_present_main = isTRUE(r$lead_main),
    qc_mdd_lead_present_noukbb = isTRUE(r$lead_noUKBB)
  )
}
out <- data.table::rbindlist(rows)
if (nrow(out) != nrow(man) || anyDuplicated(out$locus_key)) stop("Classification did not preserve one row per manifest locus")
data.table::fwrite(out, file.path(out_dir, "loci_classification.tsv"), sep = "\t")
data.table::fwrite(out[, .N, by = abf_class][order(abf_class)], file.path(out_dir, "class_summary.tsv"), sep = "\t")
