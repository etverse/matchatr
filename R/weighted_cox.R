#' Fit the Samuelsen IPW-weighted Cox model for a nested case-control sample
#'
#' Implements the Samuelsen (1997) inverse-probability-weighted partial
#' likelihood for a nested case-control (NCC) sample that was drawn with
#' [sample_ncc()]. The fitting strategy breaks the matching — controls are
#' reused across all event times rather than contributing only at their matched
#' case's failure time — and upweights each sampled control by its Samuelsen
#' inclusion weight w_j = 1/π_j. Cases retain weight 1.
#'
#' @details
#' The weighted partial likelihood (Samuelsen 1997, Biometrika 84(2)):
#'
#'   L(β) = ∏_{k: cases} exp(β x_k) /
#'            (exp(β x_k) + Σ_{j: sampled, at risk at t_k} w_j exp(β x_j))
#'
#' is identical to a standard Cox partial likelihood with observation weights
#' w_j. This is fit via [survival::coxph()] with `robust = TRUE`, which
#' uses the full cohort survival time (not `risk_time`) as the time axis, so
#' controls contribute at every event time where they are still at risk rather
#' than only at their matched set's failure time.
#'
#' The `ipw_weight` and `.cohort_row` columns are added to the NCC data by
#' [sample_ncc()] with `incl_prob = TRUE`. Rows corresponding to the same
#' cohort subject (identified by `.cohort_row`) are deduplicated before fitting
#' — a control sampled into multiple risk sets appears once in the analysis.
#'
#' The robust (Lin-Wei) sandwich variance is returned by `coxph` when
#' `robust = TRUE`; `vcov()` on the fitted object returns it directly.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"ipw_cox"`,
#'   carrying the NCC analysis `data` with `.cohort_row` and `ipw_weight`
#'   columns, plus `outcome` / `exposure` / `confounders` / design `time` slots.
#' @returns The fitted [survival::coxph] object with robust variance.
#' @family estimators
#' @seealso [matcha()], [contrast()], [sample_ncc()], [survival::coxph()]
#' @noRd
fit_ipw_cox <- function(fit) {
  ipw_col <- "ipw_weight"
  row_col <- ".cohort_row"

  if (!ipw_col %in% names(fit$data)) {
    rlang::abort(
      c(
        "The `ipw_cox` estimator requires an `ipw_weight` column in `data`.",
        i = paste0(
          "Use `sample_ncc(..., incl_prob = TRUE)` to generate the Samuelsen ",
          "KM inclusion weights before calling `matcha()`."
        )
      ),
      class = c("matchatr_missing_ipw_weights", "matchatr_error")
    )
  }
  if (!row_col %in% names(fit$data)) {
    rlang::abort(
      c(
        "The `ipw_cox` estimator requires a `.cohort_row` column in `data`.",
        i = paste0(
          "Use `sample_ncc(..., incl_prob = TRUE)` to attach the cohort row ",
          "index before calling `matcha()`."
        )
      ),
      class = c("matchatr_missing_ipw_weights", "matchatr_error")
    )
  }

  time_col <- fit$design$time
  if (is.null(time_col)) {
    rlang::abort(
      c(
        "The `ipw_cox` estimator requires a `time` column specified in `nested_cc()`.",
        i = "Use `nested_cc(strata = \"set\", time = \"t\")` with the cohort time column name."
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  }

  # Break the matching: deduplicate by cohort row index so each unique subject
  # appears once, with cohort cases forced to weight 1. The same analysis sample
  # backs the IPW Breslow absolute risk (`ipw_breslow_ncc()`).
  dt_unique <- ncc_ipw_analysis_data(fit)

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  model_formula <- stats::reformulate(
    termlabels = c(fit$exposure, conf_terms),
    response = paste0("survival::Surv(", time_col, ", ", fit$outcome, ")")
  )

  n_rows <- nrow(dt_unique)
  # `ties = "breslow"` (not the coxph default Efron) so the partial-likelihood
  # coefficients and the Breslow cumulative baseline hazard used for absolute
  # risk (`ipw_breslow_ncc()`) are mutually consistent. With Efron coefficients
  # the plain Breslow baseline disagrees at tied event times; under the
  # incidence-density sampling of an NCC the failure times are typically
  # distinct, so this matches Efron in the no-tie case and only differs (by a
  # negligible, well-understood amount) when events are tied.
  model <- survival::coxph(
    model_formula,
    data = dt_unique,
    weights = dt_unique[[ipw_col]],
    robust = TRUE,
    ties = "breslow"
  )

  # coxph's na.action silently drops rows with missing values. `model$n` is the
  # rows actually used.
  n_dropped <- n_rows - model$n
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the IPW Cox fit."
        ),
        i = "The hazard ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the IPW NCC hazard-ratio contrast
#'
#' Turns a fitted Samuelsen IPW weighted Cox model into a `matchatr_result`
#' reporting the exposure's hazard ratio with a Lin-Wei robust Wald confidence
#' interval. Only `type = "hr"` is identified; the odds ratio, risk difference,
#' and risk ratio are rejected.
#'
#' @details
#' The IPW NCC weighted Cox estimator (Samuelsen 1997) is consistent for the
#' Cox hazard ratio from the source cohort. The `coxph(robust = TRUE)` call
#' returns the Lin-Wei (1989) robust sandwich variance via `vcov()`, so
#' `ci_method = "model"` returns the robust interval. Using
#' `ci_method = "sandwich"` would re-compute the sandwich via
#' `sandwich::sandwich()` on the weighted coxph object, which is an
#' approximation; the `coxph(robust = TRUE)` variance is preferred and is the
#' default. Bootstrap is rejected as it does not account for the weight
#' estimation uncertainty from the NCC sampling.
#'
#' @param fit A `matchatr_fit` whose `model` is a `coxph` fitted by
#'   `fit_ipw_cox()`.
#' @param type Character contrast scale; only `"hr"` is computed.
#' @param ci_method Character variance source; `"model"` uses the robust
#'   sandwich from `coxph(robust = TRUE)`; `"sandwich"` re-computes via
#'   `sandwich::sandwich()` and is accepted for consistency; `"bootstrap"` is
#'   rejected.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the hazard ratio(s) for the exposure
#'   term with the Lin-Wei robust variance.
#' @family estimators
#' @seealso [contrast()], `fit_ipw_cox()`
#' @noRd
contrast_ipw_cox <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  if (identical(type, "or")) {
    rlang::abort(
      c(
        "The IPW NCC weighted Cox estimator reports hazard ratios, not odds ratios.",
        i = 'Use `type = "hr"` (the default for `ipw_cox`).'
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (identical(ci_method, "bootstrap")) {
    rlang::abort(
      c(
        '`ci_method = "bootstrap"` is not available for the IPW Cox estimator.',
        i = paste0(
          'Use `ci_method = "model"` (the robust Lin-Wei interval from ',
          '`coxph(robust = TRUE)`) or `ci_method = "sandwich"`.'
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # coxph(robust = TRUE) stores the robust sandwich in vcov(); we call
  # conditional_or_result() with robust = FALSE so it reads vcov() directly
  # rather than calling sandwich::sandwich() a second time.
  robust_flag <- identical(ci_method, "sandwich")

  conditional_or_result(
    fit,
    model = fit$model,
    robust = robust_flag,
    ci_method = ci_method,
    conf_level = conf_level,
    type = "hr",
    estimand = "hazard ratio",
    n = fit$model$n,
    call = call
  )
}

#' Deduplicated, case-weighted analysis sample for an IPW nested case-control fit
#'
#' Breaks the matching of an NCC dataset drawn by `sample_ncc(incl_prob = TRUE)`:
#' a control sampled into several risk sets carries identical cohort-level rows,
#' so the data is deduplicated by `.cohort_row` to keep each unique subject once,
#' and every subject ascertained with probability 1 is forced to weight 1
#' regardless of which duplicate row deduplication retained.
#'
#' @details
#' A subject is ascertained with probability 1 (hence weight 1) when it is
#' either a case of the analysed endpoint (`outcome == 1`) **or** the failing
#' subject of some sampled risk set (`case == 1` in any set). The second clause
#' matters for the multiple-endpoint reuse: when one control set is reused to fit
#' a *different* endpoint, the primary endpoint's cases are competing events for
#' the new analysis — they are not the analysed outcome, but they were ascertained
#' by the sampling (each anchors a risk set), so they must keep weight 1 rather
#' than the control weight 1/π_j they would carry on a row where they happened to
#' be sampled as a control. The ascertained set is read from the *pre*-dedup data
#' so the weight is correct even when deduplication retains a subject's control
#' row over its case row. For the single-endpoint analysis the two clauses
#' coincide (the sampling cases are exactly the outcome cases), so this is a
#' no-op generalisation of the original case-weight-1 rule.
#'
#' This is the sample both `fit_ipw_cox()` (the weighted Cox) and
#' `ipw_breslow_ncc()` (the IPW Breslow cumulative baseline hazard) operate on,
#' so it is computed once here.
#'
#' @param fit A `matchatr_fit` whose `data` carries the `.cohort_row` and
#'   `ipw_weight` columns, a binary `outcome` column, and (for the reuse
#'   generalisation) the per-set `case` indicator.
#' @returns A data frame with one row per unique cohort subject, the `ipw_weight`
#'   column overridden to 1 for every ascertained subject.
#' @family estimators
#' @seealso `fit_ipw_cox()`, `ipw_breslow_ncc()`, [reuse_ncc_endpoint()]
#' @noRd
ncc_ipw_analysis_data <- function(fit) {
  dt <- fit$data
  # Cohort rows ascertained with probability 1: cases of the analysed endpoint,
  # plus the failing subject of any sampled risk set (the `case == 1` rows).
  # Computed before deduplication so a subject that also appears as a control is
  # still recognised as ascertained.
  outcome_pos <- as.logical(dt[[fit$outcome]])
  if ("case" %in% names(dt)) {
    sampling_case <- !is.na(dt[["case"]]) & dt[["case"]] == 1L
    ascertained <- outcome_pos | sampling_case
  } else {
    ascertained <- outcome_pos
  }
  ascertained_rows <- unique(dt[[".cohort_row"]][ascertained])

  dup_rows <- duplicated(dt[[".cohort_row"]])
  dt_unique <- dt[!dup_rows, , drop = FALSE]
  is_ascertained <- dt_unique[[".cohort_row"]] %in% ascertained_rows
  dt_unique[["ipw_weight"]][is_ascertained] <- 1.0
  dt_unique
}

#' Fit the counter-matched weighted partial likelihood via coxph + offset
#'
#' Fits the Langholz-Borgan (1995) counter-matched partial likelihood for a
#' counter-matched NCC sample via [survival::coxph()]: the model is
#' `outcome ~ exposure + confounders + strata(set) + offset(log_w)`, where
#' `log_w` carries the log-sampling-weight for each observation. The weights
#' are the inverse inclusion probabilities within each surrogate stratum:
#' the case represents its entire surrogate stratum (weight = stratum size),
#' and each sampled control from the opposite stratum represents that
#' stratum's at-risk count divided by the number of controls drawn (weight =
#' n_other / m). The weighted partial likelihood concentrates power near the
#' exposure-surrogate correlation and identifies the hazard ratio.
#'
#' @details
#' The log-weights are a Cox offset: the conditional probability that subject
#' k is the case given the sampled risk set S̃_j is
#'
#'   exp(beta x_k + log_w_k) / sum_{l in S̃_j} exp(beta x_l + log_w_l)
#'
#' which for the case equals exp(beta x_case + log(n_z_case)) and for a
#' control from the opposite stratum equals exp(beta x_ctrl + log(n_other/m)).
#' Under a Cox PH model this is a proper partial likelihood (Borgan, Goldstein
#' & Langholz 1995, Ann. Stat.) consistent for the hazard ratio.
#'
#' `survival::coxph()` is used directly (not the `clogit()` wrapper) because
#' clogit does not pass an `offset` argument through to its internal coxph
#' call. The response is coded `Surv(rep(1, n), case_ind)` — the same
#' constant-time construction clogit uses internally — so within each stratum
#' the partial likelihood conditions on who "fails" at time 1, exactly
#' reproducing the conditional likelihood structure without offsets for uniform
#' sampling, and the counter-matching weights for non-uniform sampling.
#'
#' Rows with missing outcome, exposure, confounder, or log-weight are dropped
#' by coxph's default na.action; a `matchatr_dropped_rows` warning reports
#' how many.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"weighted_cox"`,
#'   carrying the analysis `data`, the `outcome` / `exposure` column names,
#'   the `confounders` formula (or `NULL`), and the design's `strata`,
#'   `time`, and `weights` slots.
#' @returns The fitted `survival::coxph` object.
#' @family estimators
#' @seealso [matcha()], [contrast()], [sample_ncc_counter_matched()],
#'   [survival::coxph()]
#' @noRd
fit_weighted_cox <- function(fit) {
  weights_col <- fit$design$weights
  if (is.null(weights_col)) {
    rlang::abort(
      c(
        "A counter-matched design requires a `weights` column.",
        i = paste0(
          "Supply the name of the log-weight column via ",
          "`counter_matched(weights = \"log_w\")`. ",
          "`sample_ncc_counter_matched()` appends it automatically."
        )
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  }

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # `strata()` must be unqualified: coxph treats it as a formula special by
  # name; a namespaced `survival::strata` would not be recognised.
  strata_term <- paste0(
    "strata(",
    paste(fit$design$strata, collapse = ", "),
    ")"
  )
  offset_term <- paste0("offset(", weights_col, ")")

  # Response is Surv(rep(1, n), case_ind) — the same constant-time coding that
  # clogit() uses internally for each stratum's conditional likelihood. This
  # makes everyone at risk at t=1 within each stratum; the case is the event,
  # controls are censored. The offset then enters the counter-matching weights
  # into the denominator of the within-stratum partial likelihood.
  n_rows <- nrow(fit$data)
  model_formula <- stats::reformulate(
    termlabels = c(fit$exposure, conf_terms, strata_term, offset_term),
    response = paste0(
      "survival::Surv(rep(1, ",
      n_rows,
      "), ",
      fit$outcome,
      ")"
    )
  )

  model <- survival::coxph(model_formula, data = fit$data)

  # coxph's na.action silently drops rows with missing values. `model$n` is
  # the rows actually used (coxph's nobs() counts events, not rows).
  n_dropped <- n_rows - model$n
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The hazard ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the counter-matched hazard-ratio contrast
#'
#' Turns a fitted counter-matched weighted partial likelihood into a
#' `matchatr_result` reporting the exposure's hazard ratio with a Wald
#' confidence interval. The variance is the weighted partial-likelihood
#' information matrix returned by [survival::coxph()]. Only `type = "hr"` is
#' identified; the odds ratio, risk difference, and risk ratio are rejected.
#'
#' @details
#' The counter-matched partial likelihood (Langholz & Borgan 1995) is
#' consistent for the Cox hazard ratio; the variance returned by `coxph` is the
#' inverse of the weighted information matrix and is the correct asymptotic
#' variance for this estimator. The Wald interval is on the log scale and
#' exponentiated: `estimate +/- z * se` does not reconstruct the interval on
#' the HR scale. `ci_method = "sandwich"` and `"bootstrap"` are rejected
#' because the cluster-robust and resampling variance for reused-control designs
#' belong to the inclusion-weight / IPW analyses (Phase 7), not this classical
#' conditional analysis.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted `survival::coxph`.
#' @param type Character contrast scale; `"hr"` is computed, while `"or"`,
#'   `"difference"`, and `"ratio"` abort with `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source; only `"model"` (the
#'   weighted partial-likelihood information matrix) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the hazard ratio(s) for the exposure
#'   term, with the weighted partial-likelihood information-matrix variance.
#' @family estimators
#' @seealso [contrast()], `fit_weighted_cox()`
#' @noRd
contrast_weighted_cox <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  # Counter-matching identifies the hazard ratio. The risk difference and risk
  # ratio need q0; an odds ratio confuses the estimand from a risk-set design.
  reject_unidentified_rd_rr(type, call = call)
  if (identical(type, "or")) {
    rlang::abort(
      c(
        "A counter-matched design is reported on the hazard-ratio scale.",
        i = paste0(
          "The counter-matched partial likelihood identifies the hazard ratio ",
          "(Langholz & Borgan 1995). Use `type = \"hr\"` (the default)."
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the counter-matched estimator."
        ),
        i = paste0(
          "It reports the weighted partial-likelihood information-matrix ",
          "interval. Use `ci_method = \"model\"` (the default)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  conditional_or_result(
    fit,
    model = fit$model,
    robust = FALSE,
    ci_method = ci_method,
    conf_level = conf_level,
    type = "hr",
    estimand = "hazard ratio",
    n = fit$model$n,
    call = call
  )
}
