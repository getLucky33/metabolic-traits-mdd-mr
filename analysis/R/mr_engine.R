num_p <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  suppressWarnings(as.numeric(gsub("<", "", as.character(x), fixed = TRUE)))
}

fit_ivw_random <- function(harm, exposure_name, outcome_name) {
  mri <- MendelianRandomization::mr_input(
    bx = harm$beta.exposure,
    bxse = harm$se.exposure,
    by = harm$beta.outcome,
    byse = harm$se.outcome,
    exposure = exposure_name,
    outcome = outcome_name
  )
  MendelianRandomization::mr_ivw(mri, model = "random")
}

run_mr_suite <- function(harm, exposure_name, outcome_name,
                         presso_nb = 1000L, presso_seed = 20260815L,
                         run_presso = TRUE) {
  harm <- harm[harm$mr_keep %in% TRUE, , drop = FALSE]
  if (nrow(harm) < 2L) return(list(error = "fewer_than_two_instruments", harm = harm))
  out <- list(harm = harm)
  out$ivw_random <- tryCatch(
    fit_ivw_random(harm, exposure_name, outcome_name),
    error = function(e) NULL
  )
  out$mr_long <- tryCatch(
    TwoSampleMR::mr(harm, method_list = c(
      "mr_ivw", "mr_egger_regression", "mr_weighted_median",
      "mr_weighted_mode", "mr_simple_mode"
    )),
    error = function(e) NULL
  )
  out$q <- tryCatch(TwoSampleMR::mr_heterogeneity(harm), error = function(e) NULL)
  out$egger <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harm), error = function(e) NULL)
  out$loo <- tryCatch(TwoSampleMR::mr_leaveoneout(harm), error = function(e) NULL)
  out$presso <- NULL
  if (isTRUE(run_presso) && nrow(harm) >= 4L) {
    set.seed(presso_seed)
    out$presso <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = as.data.frame(harm), NbDistribution = presso_nb,
        SignifThreshold = 0.05
      ),
      error = function(e) NULL
    )
  }
  out
}

summarise_suite <- function(suite, trait, presso_nb) {
  harm <- suite$harm
  ivw <- suite$ivw_random
  b <- if (is.null(ivw)) NA_real_ else ivw@Estimate
  se <- if (is.null(ivw)) NA_real_ else ivw@StdError
  p <- if (is.null(ivw)) NA_real_ else ivw@Pvalue
  global_p <- outlier_n <- distortion_p <- raw_b <- raw_se <- raw_p <- NA_real_
  corrected_b <- corrected_se <- corrected_p <- NA_real_
  presso_se_status <- "mr_presso_not_run_or_failed"
  if (!is.null(suite$presso)) {
    pr <- suite$presso$`MR-PRESSO results`
    global_p <- num_p(pr$`Global Test`$Pvalue)
    ot <- pr$`Outlier Test`
    if (!is.null(ot) && nrow(ot) && "Pvalue" %in% names(ot)) {
      outlier_n <- sum(grepl("^<", ot$Pvalue) |
                         suppressWarnings(as.numeric(ot$Pvalue)) <= 0.05, na.rm = TRUE)
    }
    distortion_p <- num_p(pr$`Distortion Test`$Pvalue)
    mrow <- suite$presso$`Main MR results`
    if (!is.null(mrow) && nrow(mrow)) {
      raw_b <- mrow$`Causal Estimate`[[1]]
      if ("Sd" %in% names(mrow)) raw_se <- as.numeric(mrow$Sd[[1]])
      raw_p <- num_p(mrow$`P-value`[[1]])
      presso_se_status <- if (is.finite(raw_se)) "available" else "main_mr_results_sd_missing"
      last <- nrow(mrow)
      if (last >= 2L && grepl("Outlier", mrow$`MR Analysis`[[last]], ignore.case = TRUE)) {
        corrected_b <- mrow$`Causal Estimate`[[last]]
        if ("Sd" %in% names(mrow)) corrected_se <- as.numeric(mrow$Sd[[last]])
        corrected_p <- num_p(mrow$`P-value`[[last]])
      }
    }
  }
  data.table::data.table(
    trait = trait, n_iv = nrow(harm), ivw_b = b, ivw_se = se, ivw_p = p,
    ivw_or = exp(b), ivw_or_lci = exp(b - 1.96 * se), ivw_or_uci = exp(b + 1.96 * se),
    presso_global_p = global_p, presso_outlier_n = as.integer(outlier_n),
    presso_distortion_p = distortion_p,
    presso_raw_b = raw_b, presso_raw_se = raw_se,
    presso_raw_ci_lo = raw_b - 1.96 * raw_se, presso_raw_ci_hi = raw_b + 1.96 * raw_se,
    presso_raw_p = raw_p,
    presso_corrected_b = corrected_b, presso_corrected_se = corrected_se,
    presso_corrected_ci_lo = corrected_b - 1.96 * corrected_se,
    presso_corrected_ci_hi = corrected_b + 1.96 * corrected_se,
    presso_corrected_p = corrected_p, presso_se_status = presso_se_status,
    presso_nb = as.integer(presso_nb)
  )
}

suite_diagnostics <- function(suite, trait) {
  q <- suite$q
  qrow <- if (!is.null(q) && nrow(q)) q[q$method == "Inverse variance weighted", , drop = FALSE] else NULL
  egger <- suite$egger
  data.table::data.table(
    trait = trait,
    q_ivw = if (is.null(qrow) || !nrow(qrow)) NA_real_ else qrow$Q[[1]],
    q_p = if (is.null(qrow) || !nrow(qrow)) NA_real_ else qrow$Q_pval[[1]],
    egger_intercept = if (is.null(egger) || !nrow(egger)) NA_real_ else egger$egger_intercept[[1]],
    egger_intercept_se = if (is.null(egger) || !nrow(egger)) NA_real_ else egger$se[[1]],
    egger_intercept_p = if (is.null(egger) || !nrow(egger)) NA_real_ else egger$pval[[1]]
  )
}
