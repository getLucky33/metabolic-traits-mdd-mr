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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ivm <- data.table::fread(iv_manifest_file)
assert_columns(ivm, c("trait", "iv_file"), "IV manifest")
if (length(traits_keep)) ivm <- ivm[trait %in% traits_keep]
pgc <- data.table::fread(pgc_summary_file)
assert_columns(pgc, c("trait", "ivw_b", "ivw_se", "ivw_p"), "PGC MR summary")

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

rows <- methods <- qc <- list()
for (i in seq_len(nrow(ivm))) {
  tr <- ivm$trait[[i]]
  iv <- data.table::fread(assert_file(ivm$iv_file[[i]], paste0("IV set for ", tr)))
  assert_columns(iv, c("rsid", "variant_id", "effect_allele", "other_allele", "beta",
                       "standard_error", "effect_allele_frequency", "p"), tr)
  iv[, c("chr", "pos") := data.table::tstrsplit(variant_id, "_", fixed = TRUE)[1:2]]
  iv[, poskey := paste(chr, pos, sep = ":")]
  cand <- fg[poskey %in% iv$poskey]
  cand[, rsid := strsplit(rsids, ",", fixed = TRUE)]
  cand <- cand[, .(rsid = unlist(rsid)), by = .(`#chrom`, pos, ref, alt, pval, beta, sebeta, af_alt, poskey)]
  m <- merge(iv, cand, by = c("poskey", "rsid"))
  m[, allele_ok := paste(pmin(ref, alt), pmax(ref, alt)) ==
                         paste(pmin(effect_allele, other_allele), pmax(effect_allele, other_allele))]
  m <- m[allele_ok]
  if (anyDuplicated(m$rsid)) stop(tr, ": duplicate FinnGen match; fail closed")
  exp <- TwoSampleMR::format_data(as.data.frame(m), type = "exposure", snp_col = "rsid",
    beta_col = "beta.x", se_col = "standard_error", effect_allele_col = "effect_allele",
    other_allele_col = "other_allele", eaf_col = "effect_allele_frequency", pval_col = "p")
  out <- TwoSampleMR::format_data(as.data.frame(m), type = "outcome", snp_col = "rsid",
    beta_col = "beta.y", se_col = "sebeta", effect_allele_col = "alt",
    other_allele_col = "ref", eaf_col = "af_alt", pval_col = "pval")
  h <- TwoSampleMR::harmonise_data(exp, out, action = 2)
  h <- h[h$mr_keep %in% TRUE, , drop = FALSE]
  qc[[i]] <- data.table::data.table(trait = tr, n_iv = nrow(iv), n_position_match = nrow(cand),
                                     n_allele_match = nrow(m), n_harmonised = nrow(h))
  if (nrow(h) < 2L) next
  s <- run_mr_suite(h, tr, "FinnGen R13 register-based depression", 0L, 20260815L, FALSE)
  rr <- summarise_suite(s, tr, 0L)
  prow <- pgc[trait == tr]
  if (nrow(prow) != 1L) stop("PGC summary must contain exactly one row for ", tr)
  ht <- cross_het(rr$ivw_b, rr$ivw_se, prow$ivw_b[[1]], prow$ivw_se[[1]])
  ml <- data.table::as.data.table(s$mr_long)
  if (nrow(ml)) { ml[, trait := tr]; methods[[i]] <- ml }
  rows[[i]] <- rr[, .(trait, n_iv, fg_b = ivw_b, fg_se = ivw_se, fg_p = ivw_p,
                       fg_or = ivw_or, fg_or_lci = ivw_or_lci, fg_or_uci = ivw_or_uci)][,
    `:=`(pgc_b = prow$ivw_b[[1]], pgc_se = prow$ivw_se[[1]], pgc_p = prow$ivw_p[[1]],
         dir_same = sign(fg_b) == sign(pgc_b), het_Q = ht[["Q"]], het_p = ht[["p"]], I2_pct = ht[["I2"]])]
}

res <- data.table::rbindlist(rows, fill = TRUE)
res[, `:=`(p_bh = p.adjust(fg_p, "BH"), p_bonf = p.adjust(fg_p, "bonferroni", n = n_tests))]
res[, replication_label := data.table::fcase(
  dir_same & p_bh < 0.05 & (fg_or_lci > 1 | fg_or_uci < 1) & het_p >= 0.05, "strong",
  dir_same & fg_p < 0.05 & het_p >= 0.05, "supportive",
  default = "not_replicated"
)]
data.table::fwrite(res, file.path(out_dir, "finngen_r13_mr_results.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(methods, fill = TRUE), file.path(out_dir, "finngen_r13_mr_methods.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(qc, fill = TRUE), file.path(out_dir, "finngen_r13_admission_and_tool_loss.tsv"), sep = "\t")
write_run_metadata(out_dir, "finngen", list(finngen = finngen_file, iv_manifest = iv_manifest_file,
                                             pgc_summary = pgc_summary_file),
                   list(harmonise_action = 2L, n_tests = n_tests))
