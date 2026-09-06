#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages(c("data.table", "ieugwasr"))

manifest_file <- assert_file(arg_required(args, "manifest"), "source manifest")
ld_bfile <- arg_required(args, "ld-bfile")
plink_bin <- assert_file(arg_required(args, "plink"), "PLINK executable")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p_common <- arg_value(args, "p-common", 5e-8, "numeric")
p_rare <- arg_value(args, "p-rare", 6.25e-10, "numeric")
rare_maf <- arg_value(args, "rare-maf", 0.001, "numeric")
clump_r2 <- arg_value(args, "clump-r2", 0.001, "numeric")
clump_kb <- arg_value(args, "clump-kb", 10000L, "integer")
f_min <- arg_value(args, "f-min", 10, "numeric")

if (!file.exists(paste0(ld_bfile, ".bim"))) stop("LD reference .bim not found: ", ld_bfile)
panel_rsids <- unique(data.table::fread(paste0(ld_bfile, ".bim"), select = 2L, header = FALSE)[[1]])
manifest <- data.table::fread(manifest_file)
assert_columns(manifest, c("trait", "source_file", "source_type"), "source manifest")
if (!"rsid_map_file" %in% names(manifest)) manifest[, rsid_map_file := ""]
if (length(traits_keep)) manifest <- manifest[trait %in% traits_keep]
if (!nrow(manifest)) stop("No traits selected")

read_standard_source <- function(row) {
  f <- assert_file(row$source_file, paste0("GWAS for ", row$trait))
  type <- tolower(row$source_type)
  if (type == "metabolite") {
    cols <- c("variant_id", "effect_allele", "other_allele", "beta", "standard_error",
              "effect_allele_frequency", "neg_log_10_p_value")
    d <- data.table::fread(f, select = cols, showProgress = TRUE)
    assert_columns(d, cols, row$trait)
    d[, p := 10^(-neg_log_10_p_value)]
    if (!nzchar(row$rsid_map_file)) stop(row$trait, ": rsid_map_file is required")
    mp <- data.table::fread(assert_file(row$rsid_map_file, "rsID map"))
    assert_columns(mp, c("variant_id", "rsid", "rsid_status"), "rsID map")
    mp <- mp[rsid_status == "matched", .(variant_id, rsid)]
    if (anyDuplicated(mp$variant_id)) stop(row$trait, ": duplicate variant_id in rsID map")
    d <- merge(d, mp, by = "variant_id", all.x = FALSE)
  } else if (type == "mdd") {
    d0 <- data.table::fread(f, showProgress = TRUE)
    assert_columns(d0, c("CHR", "BP", "SNP", "A1", "A2", "OR", "SE", "P"), row$trait)
    eaf_col <- grep("^FRQ_A_", names(d0), value = TRUE)[1]
    if (is.na(eaf_col)) stop(row$trait, ": no FRQ_A_* column")
    d <- d0[, .(
      variant_id = paste(CHR, BP, A2, A1, sep = "_"), rsid = SNP,
      effect_allele = A1, other_allele = A2, beta = log(OR),
      standard_error = SE, effect_allele_frequency = get(eaf_col), p = P
    )]
    rm(d0); gc()
  } else stop("Unsupported source_type for ", row$trait, ": ", row$source_type)
  d
}

qc <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  tr <- row$trait
  message("Selecting instruments for ", tr)
  d <- read_standard_source(row)
  n_input <- nrow(d)
  d <- d[is.finite(beta) & is.finite(standard_error) & standard_error > 0 &
           is.finite(p) & p >= 0 & p <= 1 &
           effect_allele_frequency > 0 & effect_allele_frequency < 1]
  d[, maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]
  if (tolower(row$source_type) == "metabolite") {
    d <- d[(maf >= rare_maf & p < p_common) | (maf < rare_maf & p < p_rare)]
  } else {
    d <- d[p < p_common]
  }
  n_after_significance_maf <- nrow(d)
  d[, chr := suppressWarnings(as.integer(sub("_.*", "", variant_id)))]
  d <- d[chr %in% 1:22 & nchar(effect_allele) == 1L & nchar(other_allele) == 1L]
  n_autosomal_biallelic <- nrow(d)
  missing_panel <- d[!rsid %in% panel_rsids]
  if (nrow(missing_panel)) {
    data.table::fwrite(missing_panel[, .(trait = tr, rsid, variant_id, p, maf)],
                       file.path(out_dir, paste0(safe_trait(tr), "_ld_panel_missing.tsv")), sep = "\t")
  }
  d <- d[rsid %in% panel_rsids]
  n_in_panel <- nrow(d)
  if (anyDuplicated(d$rsid)) d <- d[order(p)][!duplicated(rsid)]
  if (nrow(d) < 2L) stop(tr, ": fewer than two eligible variants before clumping")
  clumped <- ieugwasr::ld_clump(
    data.frame(rsid = d$rsid, pval = d$p, id = tr),
    bfile = ld_bfile, plink_bin = plink_bin,
    clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = 1
  )
  d <- d[rsid %in% clumped$rsid]
  n_clumped <- nrow(d)
  d[, F := (beta / standard_error)^2]
  d <- d[F >= f_min]
  d[, r2_approx := 2 * maf * (1 - maf) * beta^2]
  out <- d[, .(rsid, variant_id, effect_allele, other_allele, beta,
                standard_error, effect_allele_frequency, p, maf, F, r2_approx)]
  out_file <- file.path(out_dir, paste0(safe_trait(tr), "_final_ivs.tsv"))
  data.table::fwrite(out, out_file, sep = "\t")
  qc[[length(qc) + 1L]] <- data.table::data.table(
    trait = tr, input_rows = n_input, after_significance_maf = n_after_significance_maf,
    autosomal_biallelic = n_autosomal_biallelic, in_ld_panel = n_in_panel,
    ld_panel_missing = nrow(missing_panel), after_clumping = n_clumped, final_iv_f_ge_10 = nrow(out),
    source_sha256 = sha256_file(row$source_file)
  )
}
data.table::fwrite(data.table::rbindlist(qc), file.path(out_dir, "instrument_selection_qc.tsv"), sep = "\t")
write_run_metadata(out_dir, "instrument_selection", list(manifest = manifest_file), list(
  p_common = p_common, p_rare = p_rare, rare_maf = rare_maf,
  clump_r2 = clump_r2, clump_kb = clump_kb, f_min = f_min
))
