#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages("data.table")

locus_manifest_file <- assert_file(arg_required(args, "locus-manifest"), "locus manifest")
source_manifest_file <- assert_file(arg_required(args, "source-manifest"), "source manifest")
out_dir <- arg_required(args, "out-dir")
traits_keep <- csv_values(args[["traits"]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

loci <- data.table::fread(locus_manifest_file)
assert_columns(loci, c("locus_key", "trait", "chr", "window_start", "window_end"),
               "locus manifest")
sources <- data.table::fread(source_manifest_file)
assert_columns(sources, c("trait", "source_file", "rsid_map_file"), "source manifest")
if (anyDuplicated(sources$trait)) stop("Source manifest has duplicate traits")
if (length(traits_keep)) loci <- loci[trait %in% traits_keep]
if (!nrow(loci)) stop("No loci selected")

source_columns <- c(
  "chromosome", "base_pair_location", "variant_id", "effect_allele", "other_allele",
  "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value"
)
qc <- list()
for (trait_name in unique(loci$trait)) {
  spec <- sources[trait == trait_name]
  if (nrow(spec) != 1L) stop("Source manifest must contain one row for ", trait_name)
  source_file <- assert_file(spec$source_file[[1]], paste0("GWAS source for ", trait_name))
  rsid_map_file <- assert_file(spec$rsid_map_file[[1]], paste0("rsID map for ", trait_name))
  message("Reading one GWAS source for ", trait_name)
  source <- data.table::fread(source_file, select = source_columns, showProgress = TRUE)
  assert_columns(source, source_columns, trait_name)
  data.table::setkey(source, chromosome, base_pair_location)
  rsid_map <- data.table::fread(rsid_map_file)
  assert_columns(rsid_map, c("variant_id", "rsid", "rsid_status"), "rsID map")
  rsid_map <- rsid_map[rsid_status == "matched", .(variant_id, rsid)]
  if (anyDuplicated(rsid_map$variant_id)) stop(trait_name, ": duplicate mapped variant IDs")
  trait_loci <- loci[trait == trait_name]
  for (j in seq_len(nrow(trait_loci))) {
    row <- trait_loci[j]
    window <- source[
      chromosome == row$chr & base_pair_location >= row$window_start &
        base_pair_location < row$window_end
    ]
    source_n <- nrow(window)
    window <- merge(window, rsid_map, by = "variant_id")
    window <- window[, .(
      rsid, variant_id, chr = chromosome, pos = base_pair_location,
      effect_allele, other_allele, beta, se = standard_error,
      eaf = effect_allele_frequency, neg_log10p = neg_log_10_p_value
    )]
    out_file <- file.path(out_dir, paste0(row$locus_key, "_window.tsv"))
    data.table::fwrite(window, out_file, sep = "\t")
    qc[[length(qc) + 1L]] <- data.table::data.table(
      locus_key = row$locus_key, trait = trait_name,
      source_window_rows = source_n, matched_rsid_rows = nrow(window)
    )
  }
  rm(source, rsid_map)
  gc()
}
data.table::fwrite(data.table::rbindlist(qc),
                   file.path(out_dir, "locus_window_extraction_qc.tsv"), sep = "\t")
