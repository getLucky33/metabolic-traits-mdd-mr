#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages(c("data.table", "TwoSampleMR"))

iv_manifest_file <- assert_file(arg_required(args, "iv-manifest"), "IV manifest")
mdd_main_file <- assert_file(arg_required(args, "mdd-main"), "main MDD GWAS")
mdd_noukbb_file <- assert_file(arg_required(args, "mdd-noukbb"), "UK Biobank-excluded MDD GWAS")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.table::fread(iv_manifest_file)
assert_columns(manifest, c("trait", "iv_file"), "IV manifest")
if (length(traits_keep)) manifest <- manifest[trait %in% traits_keep]
if (!nrow(manifest)) stop("No IV sets selected")

read_mdd <- function(file) {
  h <- data.table::fread(file, nrows = 0L)
  frq_a <- grep("^FRQ_A_", names(h), value = TRUE)[1]
  frq_u <- grep("^FRQ_U_", names(h), value = TRUE)[1]
  cols <- c("SNP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco", frq_a, frq_u)
  if (anyNA(c(frq_a, frq_u))) stop("MDD GWAS needs FRQ_A_* and FRQ_U_* columns")
  d <- data.table::fread(file, select = cols, showProgress = TRUE)
  assert_columns(d, cols, basename(file))
  d[, EAF := (get(frq_a) * Nca + get(frq_u) * Nco) / (Nca + Nco)]
  d[, .(ID = SNP, EA = A1, NEA = A2, BETA = log(OR), SE, PVAL = P, EAF, Nca, Nco)]
}

message("Reading each MDD GWAS once")
main <- read_mdd(mdd_main_file)
noukbb <- read_mdd(mdd_noukbb_file)
data.table::setkey(main, ID)
data.table::setkey(noukbb, ID)

qc_stat <- function(h) {
  lg <- attr(h, "log")
  data.table::data.table(
    n_total = nrow(h), n_keep = sum(h$mr_keep %in% TRUE),
    n_palindromic = sum(h$palindromic %in% TRUE),
    n_ambiguous = if (is.null(lg)) NA_integer_ else lg$ambiguous_alleles,
    n_incompatible = if (is.null(lg)) NA_integer_ else lg$incompatible_alleles
  )
}

harmonise_one <- function(iv, out, action) {
  matched <- out[iv$rsid, nomatch = NULL]
  exp <- TwoSampleMR::format_data(
    as.data.frame(iv), type = "exposure", snp_col = "rsid", beta_col = "beta",
    se_col = "standard_error", effect_allele_col = "effect_allele",
    other_allele_col = "other_allele", eaf_col = "effect_allele_frequency", pval_col = "p"
  )
  y <- TwoSampleMR::format_data(
    as.data.frame(matched), type = "outcome", snp_col = "ID", beta_col = "BETA",
    se_col = "SE", effect_allele_col = "EA", other_allele_col = "NEA",
    eaf_col = "EAF", pval_col = "PVAL"
  )
  h <- TwoSampleMR::harmonise_data(exp, y, action = action)
  h$ncase <- matched$Nca[match(h$SNP, matched$ID)]
  h$ncontrol <- matched$Nco[match(h$SNP, matched$ID)]
  h
}

qc <- list()
for (i in seq_len(nrow(manifest))) {
  tr <- manifest$trait[[i]]
  iv <- data.table::fread(assert_file(manifest$iv_file[[i]], paste0("IV set for ", tr)))
  assert_columns(iv, c("rsid", "effect_allele", "other_allele", "beta", "standard_error",
                       "effect_allele_frequency", "p"), tr)
  for (dataset in c("main", "noukbb")) {
    out <- if (dataset == "main") main else noukbb
    for (action in c(2L, 3L)) {
      h <- harmonise_one(iv, out, action)
      tag <- paste(dataset, if (action == 3L) "a3" else NULL, sep = "_")
      tag <- sub("_$", "", tag)
      saveRDS(h, file.path(out_dir, paste0("harmonised_", tag, "_", safe_trait(tr), ".rds")))
      q <- qc_stat(h)
      q[, `:=`(trait = tr, dataset = dataset, action = action, input_iv = nrow(iv))]
      qc[[length(qc) + 1L]] <- q
    }
  }
}
data.table::fwrite(data.table::rbindlist(qc, fill = TRUE),
                   file.path(out_dir, "harmonisation_qc.tsv"), sep = "\t")
write_run_metadata(out_dir, "forward_harmonisation",
                   list(iv_manifest = iv_manifest_file, mdd_main = mdd_main_file,
                        mdd_noukbb = mdd_noukbb_file),
                   list(action_primary = 2L, action_palindrome_sensitivity = 3L))
