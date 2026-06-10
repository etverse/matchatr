#' Fit the Samuelsen IPW-weighted accelerated failure time model for NCC data
#'
#' Fits a weighted accelerated failure time (AFT) model on a nested case-control
#' sample drawn with [sample_ncc()] (`incl_prob = TRUE`). The matching is broken:
#' each unique cohort subject (deduplicated by `.cohort_row`) enters once,
#' upweighted by its Samuelsen inclusion weight 1/π_j (cases at weight 1), and a
#' parametric AFT is fit by weighted likelihood with the Lin-Wei robust sandwich
#' variance.
#'
#' @details
#' The AFT model parameterises the log failure time directly,
#'
#'   log T = β₀ + βᵀ x + σ ε,
#'
#' so exp(β_j) is the **time ratio** (acceleration factor): a unit increase in
#' covariate j multiplies the survival time by exp(β_j) (> 1 prolongs, < 1
#' shortens survival). This is the time-scale analogue of the Cox hazard ratio
#' and a different estimand. A Weibull baseline is used (`dist = "weibull"`): it
#' is the canonical AFT distribution and the only one that is simultaneously an
#' AFT and a proportional-hazards model. Estimation is by
#' [survival::survreg()] with the IPW observation weights and `robust = TRUE`,
#' which returns the Lin-Wei robust sandwich variance appropriate when the
#' weights reuse controls across the cohort (Kang, Lu & Liu 2017, Biometrics
#' 73(1)).
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"ipw_aft"`, carrying
#'   the NCC analysis `data` with `.cohort_row` / `ipw_weight` columns plus the
#'   `outcome` / `exposure` / `confounders` slots and the design `time`.
#' @returns The fitted [survival::survreg] object (Weibull, robust variance).
#' @family estimators
#' @seealso [matcha()], [contrast()], [sample_ncc()], [survival::survreg()]
#' @noRd
fit_ipw_aft <- function(fit) {
  time_col <- require_ipw_ncc_columns(fit, "ipw_aft")

  # Break the matching: the same deduplicated, case-weight-forced analysis
  # sample the IPW Cox uses.
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
  # Weibull baseline: the canonical AFT distribution (and the only AFT that is
  # also proportional-hazards), with the Lin-Wei robust sandwich (`robust =
  # TRUE`) for the reused-control IPW weights.
  model <- survival::survreg(
    model_formula,
    data = dt_unique,
    weights = dt_unique[["ipw_weight"]],
    dist = "weibull",
    robust = TRUE
  )

  # survreg's default na.action drops rows with missing values; the linear
  # predictor has one entry per row actually used.
  n_used <- length(model$linear.predictors)
  n_dropped <- n_rows - n_used
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the IPW AFT fit."
        ),
        i = "The time ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the IPW NCC accelerated-failure-time contrast
#'
#' Turns a fitted Samuelsen IPW weighted AFT model into a `matchatr_result`
#' reporting the exposure's **time ratio** exp(β) (acceleration factor) with a
#' Lin-Wei robust Wald confidence interval. Only `type = "af"` is identified; the
#' odds ratio, hazard ratio, risk difference, and risk ratio are rejected as
#' off-scale for the time-ratio estimand.
#'
#' @details
#' The time ratio is exp(β) with a Wald interval formed on the log (β) scale and
#' exponentiated, so it is asymmetric on the time-ratio scale. The variance is
#' the robust sandwich [survival::survreg()] stores under `robust = TRUE`, so
#' `ci_method = "model"` reads it directly; `"bootstrap"` is rejected (it does
#' not account for the NCC weight estimation), and `"sandwich"` is rejected
#' because `survreg`'s robust variance is already the sandwich and
#' `sandwich::sandwich()` has no survreg method.
#'
#' @param fit A `matchatr_fit` whose `model` is a `survreg` fit by
#'   `fit_ipw_aft()`.
#' @param type Character contrast scale; only `"af"` (acceleration factor / time
#'   ratio) is computed.
#' @param ci_method Character variance source; only `"model"` (the robust
#'   sandwich from `survreg(robust = TRUE)`) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the time ratio(s) for the exposure term
#'   with the Lin-Wei robust variance.
#' @family estimators
#' @seealso [contrast()], `fit_ipw_aft()`
#' @noRd
contrast_ipw_aft <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  if (!identical(type, "af")) {
    rlang::abort(
      c(
        paste0(
          "The IPW NCC accelerated failure time model reports a time ratio, not `type = \"",
          type,
          "\"`."
        ),
        i = 'Use `type = "af"` (the default for `ipw_aft`).'
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
          "\"` is not available for the IPW AFT estimator."
        ),
        i = paste0(
          "It reports the robust Lin-Wei sandwich `survreg(robust = TRUE)` ",
          "computes. Use `ci_method = \"model\"` (the default)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # survreg(robust = TRUE) stores the robust sandwich in vcov(); pass
  # robust = FALSE so conditional_or_result() reads it directly rather than
  # calling sandwich::sandwich() (which has no survreg method). vcov() carries
  # the extra log(scale) row/column, but estimable_vcov() indexes the leading
  # regression-coefficient block by position, so the time ratio reads cleanly.
  conditional_or_result(
    fit,
    model = fit$model,
    robust = FALSE,
    ci_method = ci_method,
    conf_level = conf_level,
    type = "af",
    estimand = "time ratio",
    n = length(fit$model$linear.predictors),
    call = call
  )
}
