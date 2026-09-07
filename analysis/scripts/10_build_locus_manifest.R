#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages("data.table")

screen_file <- assert_file(arg_required(args, "forward-screen"), "forward screen")
iv_manifest_file <- assert_file(arg_required(args, "iv-manifest"), "IV manifest")
main_coords_file <- assert_file(arg_required(args, "mdd-main-coords"), "main MDD coordinates")
noukbb_coords_file <- assert_file(arg_required(args, "mdd-noukbb-coords"), "no-UKBB MDD coordinates")
r_tier_file <- assert_file(arg_required(args, "r-tier"), "reverse-tier table")
out_file <- arg_required(args, "out")
window_dir <- arg_value(args, "window-dir", "data-local/coloc/windows")
expected_candidates <- arg_value(args, "expected-candidates", 15L, "integer")
expected_loci <- arg_value(args, "expected-loci", 0L, "integer")
merge_dist <- arg_value(args, "merge-dist", 1000000, "numeric")
flank <- arg_value(args, "flank", 500000, "numeric")

screen <- data.table::fread(screen_file)
assert_columns(screen, c("trait", "screen_level"), "forward screen")
candidates <- sort(screen[screen_level == "bonferroni_hit", unique(trait)])
if (length(candidates) != expected_candidates) {
  stop("Expected ", expected_candidates, " Bonferroni candidates; found ", length(candidates))
}
ivm <- data.table::fread(iv_manifest_file)
assert_columns(ivm, c("trait", "iv_file"), "IV manifest")
if (anyDuplicated(ivm$trait)) stop("IV manifest has duplicate traits")
tiers <- data.table::fread(r_tier_file)
assert_columns(tiers, c("trait", "r_tier"), "reverse-tier table")
if (anyDuplicated(tiers$trait)) stop("Reverse-tier table has duplicate traits")
tier_map <- setNames(tiers$r_tier, tiers$trait)

read_coords <- function(path, source_name) {
  x <- data.table::fread(path)
  assert_columns(x, c("rsid", "chr", "pos"), paste0(source_name, " MDD coordinates"))
  x[, source := source_name]
  x[, .(rsid, chr = as.integer(chr), pos = as.integer(pos), source)]
}
r_all <- data.table::rbindlist(list(
  read_coords(main_coords_file, "main"),
  read_coords(noukbb_coords_file, "noUKBB")
))
if (anyNA(r_all$chr) || anyNA(r_all$pos)) stop("MDD coordinates contain missing values")
if (nrow(r_all[, .(n_rsid = data.table::uniqueN(rsid)), by = .(chr, pos)][n_rsid > 1L]) > 0L ||
    nrow(r_all[, .(n_pos = data.table::uniqueN(paste(chr, pos))), by = rsid][n_pos > 1L]) > 0L) {
  stop("MDD lead coordinates conflict; fail closed")
}
r_leads <- unique(r_all[, .(rsid, chr, pos)])

rows <- list()
for (trait_name in candidates) {
  iv_file <- ivm[trait == trait_name, iv_file]
  if (length(iv_file) != 1L) stop("IV manifest must contain one file for ", trait_name)
  iv <- data.table::fread(assert_file(iv_file, paste0("IV set for ", trait_name)))
  assert_columns(iv, c("rsid", "variant_id"), trait_name)
  iv[, `:=`(
    chr = suppressWarnings(as.integer(sub("_.*", "", variant_id))),
    pos = suppressWarnings(as.integer(sub("^[^_]*_([0-9]+)_.*", "\\1", variant_id)))
  )]
  if (anyNA(iv$chr) || anyNA(iv$pos)) stop(trait_name, ": invalid GRCh38 variant_id")
  leads <- iv[, .(rsid, chr, pos, lead_source = "forward")]
  tier <- unname(tier_map[trait_name])
  if (!length(tier) || is.na(tier)) stop("Missing reverse tier for ", trait_name)
  if (tier %in% c("confirmatory", "sensitivity")) {
    leads <- data.table::rbindlist(list(
      leads, r_leads[, .(rsid, chr, pos, lead_source = "reverse")]
    ))
  }
  leads <- leads[, .(
    forward = any(lead_source == "forward"),
    reverse = any(lead_source == "reverse")
  ), by = .(chr, pos, rsid)]
  if (nrow(leads[, .(n_rsid = data.table::uniqueN(rsid)), by = .(chr, pos)][n_rsid > 1L]) > 0L ||
      nrow(leads[, .(n_pos = data.table::uniqueN(paste(chr, pos))), by = rsid][n_pos > 1L]) > 0L) {
    stop(trait_name, ": lead-coordinate conflict; fail closed")
  }
  data.table::setorder(leads, chr, pos, rsid)
  leads[, cluster := 0L]
  current_chr <- current_anchor <- NA_integer_
  cluster <- 0L
  for (j in seq_len(nrow(leads))) {
    if (is.na(current_chr) || leads$chr[[j]] != current_chr ||
        leads$pos[[j]] - current_anchor > merge_dist) {
      cluster <- cluster + 1L
      current_chr <- leads$chr[[j]]
      current_anchor <- leads$pos[[j]]
    }
    data.table::set(leads, j, "cluster", cluster)
  }
  loci <- leads[, .(
    chr = chr[[1]], lead_rsid = rsid[[1]], lead_pos = pos[[1]], n_iv = .N,
    min_pos = min(pos), max_pos = max(pos),
    locus_source = if (any(forward) && any(reverse)) "both" else
      if (any(forward)) "forward" else "reverse"
  ), by = cluster]
  loci[, `:=`(
    window_start = pmax(1, min_pos - flank),
    window_end = max_pos + flank
  )]
  for (j in seq_len(nrow(loci))) {
    locus_id <- sprintf("locus%03d", loci$cluster[[j]])
    locus_key <- paste0(trait_name, "__", locus_id)
    rows[[length(rows) + 1L]] <- data.table::data.table(
      locus_key = locus_key, trait = trait_name, r_tier = tier, locus_id = locus_id,
      chr = loci$chr[[j]], lead_rsid = loci$lead_rsid[[j]], lead_pos = loci$lead_pos[[j]],
      n_iv = loci$n_iv[[j]], window_start = loci$window_start[[j]],
      window_end = loci$window_end[[j]],
      span = loci$window_end[[j]] - loci$window_start[[j]],
      locus_source = loci$locus_source[[j]], complex_locus = loci$n_iv[[j]] >= 2L,
      window_file = gsub("\\\\", "/", file.path(window_dir, paste0(locus_key, "_window.tsv")))
    )
  }
}

out <- data.table::rbindlist(rows)
if (anyDuplicated(out$locus_key) || any(out$span > 2 * merge_dist)) {
  stop("Locus manifest violated key uniqueness or the locked span limit")
}
if (expected_loci > 0L && nrow(out) != expected_loci) {
  stop("Expected ", expected_loci, " loci; found ", nrow(out))
}
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(out, out_file, sep = "\t")
