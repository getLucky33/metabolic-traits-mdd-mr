#!/usr/bin/env Rscript
source(file.path("analysis", "R", "common.R"))
args <- args_map()
require_packages(c("data.table", "coloc"))

locus_manifest_file <- assert_file(arg_required(args, "locus-manifest"), "locus manifest")
outcome_manifest_file <- assert_file(arg_required(args, "outcome-manifest"), "outcome manifest")
out_dir <- arg_required(args, "out-dir")
limit <- arg_value(args, "limit", 0L, "integer")
min_snps <- arg_value(args, "min-snps", 50L, "integer")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

p12_grid <- c(1e-6, 5e-6, 1e-5, 1e-4)
p1 <- 1e-4; p2 <- 1e-4
comp <- c(A = "T", T = "A", C = "G", G = "C")
is_palin <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") |
  (a == "C" & b == "G") | (a == "G" & b == "C")

loci <- data.table::fread(locus_manifest_file)
assert_columns(loci, c("locus_key", "trait", "locus_id", "lead_rsid", "window_file"), "locus manifest")
if (limit > 0L) loci <- head(loci, limit)
if (!nrow(loci)) stop("No loci selected")
outcomes <- data.table::fread(outcome_manifest_file)
assert_columns(outcomes, c("outcome", "file", "ncase", "ncontrol", "frq_case_col", "frq_control_col"),
               "outcome manifest")

needed <- character()
lead_qc <- list()
for (i in seq_len(nrow(loci))) {
  w <- data.table::fread(assert_file(loci$window_file[[i]], "colocalization window"), select = "rsid")
  needed <- c(needed, w$rsid[w$rsid != "."])
  lead_qc[[i]] <- data.table::data.table(
    locus_key = loci$locus_key[[i]], lead_rsid = loci$lead_rsid[[i]],
    qc_exp_lead_present = loci$lead_rsid[[i]] %in% w$rsid
  )
}
needed <- unique(needed)

outcome_data <- list()
for (i in seq_len(nrow(outcomes))) {
  spec <- outcomes[i]
  f <- assert_file(spec$file, paste0("outcome ", spec$outcome))
  cols <- c("SNP", "A1", "A2", "OR", "SE", spec$frq_case_col, spec$frq_control_col, "Nca", "Nco")
  message("Reading outcome once and retaining needed rsIDs: ", spec$outcome)
  d <- data.table::fread(f, select = cols, showProgress = TRUE)
  assert_columns(d, cols, spec$outcome)
  data.table::setnames(d, c("SNP", spec$frq_case_col, spec$frq_control_col), c("rsid", "frq_case", "frq_control"))
  d <- d[rsid %in% needed]
  if (anyDuplicated(d$rsid)) {
    dup <- d[duplicated(rsid) | duplicated(rsid, fromLast = TRUE)]
    conflict <- dup[, data.table::uniqueN(paste(A1, A2, OR, SE, frq_case, frq_control)), by = rsid][V1 > 1L]
    if (nrow(conflict)) stop(spec$outcome, ": conflicting duplicate rsIDs")
    d <- unique(d, by = "rsid")
  }
  d[, `:=`(
    beta_mdd = log(OR), se_mdd = SE,
    eaf_mdd = (frq_case * Nca + frq_control * Nco) / (Nca + Nco)
  )]
  outcome_data[[spec$outcome]] <- d
}

harmonise_one <- function(exp, out) {
  n_window <- nrow(exp)
  m <- merge(exp, out, by = "rsid")
  n_merged <- nrow(m)
  m <- m[effect_allele %in% names(comp) & other_allele %in% names(comp) &
           A1 %in% names(comp) & A2 %in% names(comp)]
  m[, palindromic := is_palin(effect_allele, other_allele)]
  m[, direct := (effect_allele == A1 & other_allele == A2) |
                  (effect_allele == A2 & other_allele == A1)]
  m[, complement := (effect_allele == comp[A1] & other_allele == comp[A2]) |
                      (effect_allele == comp[A2] & other_allele == comp[A1])]
  strand_drop <- sum(!m$direct & !m$complement)
  m <- m[direct | complement]
  m[complement & !direct, `:=`(A1 = comp[A1], A2 = comp[A2])]
  m[, eaf_ambig := palindromic & (abs(eaf - 0.5) < 0.08 | abs(eaf_mdd - 0.5) < 0.08)]
  pal_drop <- sum(m$eaf_ambig, na.rm = TRUE)
  m <- m[eaf_ambig == FALSE]
  m[, flip_mdd := palindromic & abs(eaf - eaf_mdd) > abs(eaf - (1 - eaf_mdd))]
  m[flip_mdd == TRUE, `:=`(A1 = A2, A2 = A1, beta_mdd = -beta_mdd, eaf_mdd = 1 - eaf_mdd)]
  m[effect_allele == A2, `:=`(beta = -beta, eaf = 1 - eaf,
                               effect_allele = A1, other_allele = A2)]
  m[, `:=`(maf_exp = pmin(eaf, 1 - eaf), maf_mdd = pmin(eaf_mdd, 1 - eaf_mdd))]
  m <- m[is.finite(beta) & is.finite(se) & se > 0 & is.finite(beta_mdd) &
           is.finite(se_mdd) & se_mdd > 0 & maf_exp > 0 & maf_exp <= 0.5 &
           maf_mdd > 0 & maf_mdd <= 0.5]
  if (anyDuplicated(m$rsid)) {
    dup <- m[duplicated(rsid) | duplicated(rsid, fromLast = TRUE)]
    conflict <- dup[, data.table::uniqueN(paste(effect_allele, other_allele, beta, se)), by = rsid][V1 > 1L]
    if (nrow(conflict)) stop("conflicting duplicate exposure rsIDs")
    m <- unique(m, by = "rsid")
  }
  list(h = m, flow = data.table::data.table(
    n_window = n_window, n_merged = n_merged, n_strand_drop = strand_drop,
    n_palindromic_ambiguous_drop = pal_drop, n_after_dedup = nrow(m)
  ))
}

results <- status <- flow <- final_snps <- list()
ri <- si <- fi <- qi <- 0L
for (i in seq_len(nrow(loci))) {
  spec <- loci[i]
  exp <- data.table::fread(spec$window_file)
  assert_columns(exp, c("rsid", "beta", "se", "eaf", "effect_allele", "other_allele"), spec$locus_key)
  exp <- exp[rsid != "."]
  for (outcome_name in names(outcome_data)) {
    hh <- tryCatch(harmonise_one(exp, outcome_data[[outcome_name]]), error = function(e) list(error = conditionMessage(e)))
    if (!is.null(hh$error)) {
      si <- si + 1L; status[[si]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, analysis_status = "model_error", failure_reason = hh$error)
      next
    }
    h <- hh$h
    qi <- qi + 1L; flow[[qi]] <- hh$flow[, `:=`(
      locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id, outcome = outcome_name)]
    if (!nrow(h)) {
      si <- si + 1L; status[[si]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, analysis_status = "no_overlap", failure_reason = "0 comparable SNPs")
      next
    }
    if (nrow(h) < min_snps) {
      si <- si + 1L; status[[si]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, analysis_status = "nsnp_insufficient",
        failure_reason = paste0("n_after_dedup=", nrow(h), " < ", min_snps))
      next
    }
    d1 <- list(snp = h$rsid, beta = h$beta, varbeta = h$se^2, type = "quant", sdY = 1)
    ospec <- outcomes[outcome == outcome_name]
    d2 <- list(snp = h$rsid, beta = h$beta_mdd, varbeta = h$se_mdd^2, type = "cc",
               s = ospec$ncase[[1]] / (ospec$ncase[[1]] + ospec$ncontrol[[1]]))
    ck <- tryCatch(c(coloc::check_dataset(d1), coloc::check_dataset(d2)), error = function(e) conditionMessage(e))
    if (is.character(ck) && length(ck)) {
      si <- si + 1L; status[[si]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, analysis_status = "dataset_invalid", failure_reason = paste(ck, collapse = " | "))
      next
    }
    cache <- list(); failed <- NULL
    for (p12 in p12_grid) {
      z <- tryCatch(coloc::coloc.abf(d1, d2, p1 = p1, p2 = p2, p12 = p12),
                    error = function(e) { failed <<- conditionMessage(e); NULL })
      if (is.null(z)) break
      pp <- z$summary
      cache[[length(cache) + 1L]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, p12 = p12, nsnp = nrow(h),
        PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]],
        PP.H2 = pp[["PP.H2.abf"]], PP.H3 = pp[["PP.H3.abf"]], PP.H4 = pp[["PP.H4.abf"]]
      )
    }
    if (length(cache) != length(p12_grid)) {
      si <- si + 1L; status[[si]] <- data.table::data.table(
        locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
        outcome = outcome_name, analysis_status = "model_error", failure_reason = failed)
      next
    }
    for (z in cache) { ri <- ri + 1L; results[[ri]] <- z }
    fi <- fi + 1L; final_snps[[fi]] <- data.table::data.table(
      locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
      outcome = outcome_name, rsid = h$rsid)
    si <- si + 1L; status[[si]] <- data.table::data.table(
      locus_key = spec$locus_key, trait = spec$trait, locus = spec$locus_id,
      outcome = outcome_name, analysis_status = "success", failure_reason = "none")
  }
}
data.table::fwrite(data.table::rbindlist(results, fill = TRUE), file.path(out_dir, "loci_coloc_results.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(status, fill = TRUE), file.path(out_dir, "loci_status.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(flow, fill = TRUE), file.path(out_dir, "loci_flow.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(final_snps, fill = TRUE), file.path(out_dir, "loci_final_snps.tsv"), sep = "\t")
data.table::fwrite(data.table::rbindlist(lead_qc), file.path(out_dir, "exposure_lead_qc.tsv"), sep = "\t")
write_run_metadata(out_dir, "coloc_abf", list(locus_manifest = locus_manifest_file,
                                               outcome_manifest = outcome_manifest_file),
                   list(p1 = p1, p2 = p2, p12_grid = paste(p12_grid, collapse = ","),
                        primary_p12 = 5e-6, min_snps = min_snps,
                        palindrome_eaf_exclusion = 0.08))
