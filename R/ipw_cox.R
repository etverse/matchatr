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
  time_col <- require_ipw_ncc_columns(fit, "ipw_cox")

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
  # The weighted Cox identifies only the hazard ratio; an odds ratio, the AFT
  # time ratio, or the additive excess hazard are sibling NCC-IPW scales fit by
  # the other engines, not relabelings of this fit.
  if (!identical(type, "hr")) {
    rlang::abort(
      c(
        paste0(
          "The IPW NCC weighted Cox estimator reports hazard ratios, not `type = \"",
          type,
          "\"`."
        ),
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

#' Validate that a fit carries the columns an IPW nested case-control engine needs
#'
#' The Samuelsen IPW engines (`ipw_cox`, `ipw_aft`, `ipw_aalen`) all break the
#' matching and fit a weighted model on the deduplicated NCC analysis sample, so
#' they share the same data contract: an `ipw_weight` inclusion-weight column, a
#' `.cohort_row` index for deduplication (both attached by
#' `sample_ncc(incl_prob = TRUE)`), and a cohort `time` column declared on the
#' `nested_cc()` design.
#'
#' @param fit A `matchatr_fit` whose `data` should carry `ipw_weight` /
#'   `.cohort_row` and whose `design$time` should be set.
#' @param estimator Character scalar naming the estimator in the error messages
#'   (e.g. `"ipw_cox"`, `"ipw_aft"`, `"ipw_aalen"`).
#' @returns The cohort time column name (character scalar); aborts with
#'   `matchatr_missing_ipw_weights` (missing weight / row columns) or
#'   `matchatr_bad_design` (missing time) otherwise.
#' @family estimators
#' @seealso `fit_ipw_cox()`, `fit_ipw_aft()`, `fit_ipw_aalen()`
#' @noRd
require_ipw_ncc_columns <- function(fit, estimator) {
  if (!"ipw_weight" %in% names(fit$data)) {
    rlang::abort(
      c(
        paste0(
          "The `",
          estimator,
          "` estimator requires an `ipw_weight` column in `data`."
        ),
        i = paste0(
          "Use `sample_ncc(..., incl_prob = TRUE)` to generate the Samuelsen ",
          "KM inclusion weights before calling `matcha()`."
        )
      ),
      class = c("matchatr_missing_ipw_weights", "matchatr_error")
    )
  }
  if (!".cohort_row" %in% names(fit$data)) {
    rlang::abort(
      c(
        paste0(
          "The `",
          estimator,
          "` estimator requires a `.cohort_row` column in `data`."
        ),
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
        paste0(
          "The `",
          estimator,
          "` estimator requires a `time` column specified in `nested_cc()`."
        ),
        i = "Use `nested_cc(strata = \"set\", time = \"t\")` with the cohort time column name."
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  }
  time_col
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
