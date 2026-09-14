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
  out$presso_seed <- as.integer(presso_seed)
  presso_error <- NA_character_
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
      error = function(e) {
        presso_error <<- conditionMessage(e)
        NULL
      }
    )
  }
  out$presso_error <- presso_error
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
  global_rssobs <- distortion_coefficient <- raw_t_stat <- corrected_t_stat <- NA_real_
  raw_df <- corrected_df <- NA_integer_
  global_exceedance_n <- NA_integer_
  raw_ci_lo <- raw_ci_hi <- corrected_ci_lo <- corrected_ci_hi <- NA_real_
  global_p_display <- distortion_p_display <- NA_character_
  regression_p_sidedness <- empirical_test_df <- NA_character_
  presso_se_status <- "mr_presso_not_run_or_failed"
  if (!is.null(suite$presso)) {
    pr <- suite$presso$`MR-PRESSO results`
    global_test <- pr$`Global Test`
    global_p_display <- if (is.null(global_test$Pvalue) || !length(global_test$Pvalue)) {
      NA_character_
    } else {
      as.character(global_test$Pvalue[[1]])
    }
    global_p <- num_p(global_p_display)
    global_rssobs <- if (is.null(global_test$RSSobs) || !length(global_test$RSSobs)) {
      NA_real_
    } else {
      as.numeric(global_test$RSSobs[[1]])
    }
    global_exceedance_n <- if (!is.na(global_p_display) && startsWith(global_p_display, "<")) 0L else NA_integer_
    ot <- pr$`Outlier Test`
    if (!is.null(ot) && nrow(ot) && "Pvalue" %in% names(ot)) {
      outlier_n <- sum(grepl("^<", ot$Pvalue) |
                         suppressWarnings(as.numeric(ot$Pvalue)) <= 0.05, na.rm = TRUE)
    }
    distortion <- pr$`Distortion Test`
    distortion_p_display <- if (is.null(distortion$Pvalue) || !length(distortion$Pvalue)) {
      NA_character_
    } else {
      as.character(distortion$Pvalue[[1]])
    }
    distortion_p <- num_p(distortion_p_display)
    distortion_coefficient <- if (is.null(distortion$`Distortion Coefficient`) ||
                                      !length(distortion$`Distortion Coefficient`)) {
      NA_real_
    } else {
      as.numeric(distortion$`Distortion Coefficient`[[1]])
    }
    mrow <- suite$presso$`Main MR results`
    if (!is.null(mrow) && nrow(mrow)) {
      raw_b <- mrow$`Causal Estimate`[[1]]
      if ("Sd" %in% names(mrow)) raw_se <- as.numeric(mrow$Sd[[1]])
      if ("T-stat" %in% names(mrow)) raw_t_stat <- as.numeric(mrow$`T-stat`[[1]])
      raw_p <- num_p(mrow$`P-value`[[1]])
      if (is.finite(raw_se) && raw_se > 0) {
        raw_df <- nrow(harm) - 1L
        raw_ci_lo <- raw_b - stats::qt(0.975, df = raw_df) * raw_se
        raw_ci_hi <- raw_b + stats::qt(0.975, df = raw_df) * raw_se
        presso_se_status <- "available"
      } else {
        presso_se_status <- "main_mr_results_sd_missing"
      }
      last <- nrow(mrow)
      if (last >= 2L && grepl("Outlier", mrow$`MR Analysis`[[last]], ignore.case = TRUE)) {
        corrected_b <- as.numeric(mrow$`Causal Estimate`[[last]])
        if ("Sd" %in% names(mrow)) corrected_se <- as.numeric(mrow$Sd[[last]])
        if ("T-stat" %in% names(mrow)) corrected_t_stat <- as.numeric(mrow$`T-stat`[[last]])
        corrected_p <- num_p(mrow$`P-value`[[last]])
        if (!is.finite(corrected_b)) {
          corrected_b <- corrected_se <- corrected_p <- corrected_t_stat <- NA_real_
        } else if (is.finite(corrected_se) && corrected_se > 0 && is.finite(outlier_n)) {
          corrected_df <- nrow(harm) - as.integer(outlier_n) - 1L
          if (corrected_df < 1L) stop("Invalid MR-PRESSO corrected residual df for ", trait)
          corrected_ci_lo <- corrected_b - stats::qt(0.975, df = corrected_df) * corrected_se
          corrected_ci_hi <- corrected_b + stats::qt(0.975, df = corrected_df) * corrected_se
        } else {
          presso_se_status <- "corrected_main_mr_results_sd_or_outlier_count_missing"
        }
      }
    }
    regression_p_sidedness <- "two-sided"
    empirical_test_df <- "not_applicable"
  }
  data.table::data.table(
    trait = trait, n_iv = nrow(harm), ivw_b = b, ivw_se = se, ivw_p = p,
    ivw_or = exp(b), ivw_or_lci = exp(b - 1.96 * se), ivw_or_uci = exp(b + 1.96 * se),
    presso_global_p = global_p, presso_global_p_display = global_p_display,
    presso_global_rssobs = global_rssobs,
    presso_global_exceedance_n = global_exceedance_n,
    presso_outlier_n = as.integer(outlier_n),
    presso_distortion_p = distortion_p,
    presso_distortion_p_display = distortion_p_display,
    presso_distortion_coefficient = distortion_coefficient,
    presso_raw_b = raw_b, presso_raw_se = raw_se, presso_raw_df = raw_df,
    presso_raw_t_stat = raw_t_stat,
    presso_raw_ci_lo = raw_ci_lo, presso_raw_ci_hi = raw_ci_hi,
    presso_raw_p = raw_p,
    presso_corrected_b = corrected_b, presso_corrected_se = corrected_se,
    presso_corrected_df = corrected_df,
    presso_corrected_t_stat = corrected_t_stat,
    presso_corrected_ci_lo = corrected_ci_lo,
    presso_corrected_ci_hi = corrected_ci_hi,
    presso_corrected_p = corrected_p,
    presso_regression_p_sidedness = regression_p_sidedness,
    presso_empirical_test_df = empirical_test_df,
    presso_se_status = presso_se_status,
    presso_nb = as.integer(presso_nb),
    presso_seed = if (is.null(suite$presso_seed)) NA_integer_ else as.integer(suite$presso_seed)
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
