#' Fit a matched case-control conditional logistic regression
#'
#' Fits the conditional maximum-likelihood (CMLE) odds ratio for an individually
#' or frequency matched case-control sample via [survival::clogit()]: the model
#' is `outcome ~ exposure + confounders + strata(set)`, with each matched set as
#' a stratum. Conditioning on the matched-set totals removes the
#' matching-variable nuisance parameters, so only the exposure / adjustment odds
#' ratios are reported; the matching variables are controlled implicitly and
#' have no estimable coefficient.
#'
#' @details
#' The conditional likelihood for a stratum with one case and M controls is
#'
#'   prod exp(x_case . beta) / sum_j exp(x_j . beta),
#'
#' which for 1:1 matching reduces to expit{(x_case - x_control) . beta}.
#' `survival::clogit` is the Cox partial likelihood with each matched set as a
#' stratum, so this is exactly the CMLE. Unconditional logistic regression on
#' matched-set indicators is never used: for 1:1 matching its MLE converges to
#' the squared odds ratio in large samples (Pike et al. 1980; Breslow & Day
#' 1980).
#'
#' The confounder formula's terms — transforms (`poly(age, 2)`), interactions
#' (`age:smoke`) — are carried through verbatim for adjustment of non-matching
#' covariates. Several matching columns are crossed into one `strata()` term
#' (frequency matching on, e.g., age group and sex). Matched sets with no case
#' or no control carry no information and are dropped by `clogit`; the
#' `matcha()` entry point already warns about them
#' (`matchatr_uninformative_stratum`). Rows with a missing outcome, exposure, or
#' confounder are dropped by the default `na.action`; a `matchatr_dropped_rows`
#' warning reports how many.
#'
#' The nested case-control risk-set analysis shares this conditional likelihood
#' (a sampled risk set is the stratum) but adds time-aware sampling and
#' inclusion weighting handled by its own design layer, so this engine fits the
#' matched case-control design only.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"clogit"`, carrying
#'   the analysis `data`, the `outcome` / `exposure` column names, the
#'   `confounders` formula (or `NULL`), and the design's matched-set `strata`.
#' @returns The fitted `survival::clogit` object (a `clogit` / `coxph`), or
#'   `NULL` when the design is not matched case-control (the conditional engine
#'   has nothing wired for it yet).
#' @family estimators
#' @seealso [matcha()], [contrast()], [survival::clogit()]
#' @noRd
fit_clogit <- function(fit) {
  # The conditional partial likelihood is the matched case-control analysis;
  # the nested case-control risk-set form (time-aware sampling, inclusion
  # weighting, hazard-ratio reporting) is wired by its own design layer, so the
  # clogit engine fits the matched design only and otherwise stays unestimated.
  if (!identical(fit$design$type, "matched_cc")) {
    return(NULL)
  }

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # `strata()` must be a bare, un-namespaced term: `survival::clogit` detects it
  # as a formula special by the name "strata", and a namespaced `survival::strata`
  # would not be recognised as the special. Several matching columns cross into
  # one stratifying factor (frequency matching).
  strata_term <- paste0(
    "strata(",
    paste(fit$design$strata, collapse = ", "),
    ")"
  )
  # With an effect modifier the exposure enters crossed with the modifier
  # (`exposure * modifier` = exposure + modifier + exposure:modifier); the
  # interaction coefficients carry the per-level shift in the exposure log OR
  # (the stratum-specific contrast). Without one the exposure is a plain main
  # effect whose coefficient(s) identify the conditional log OR.
  exposure_term <- if (is.null(fit$effect_modifier)) {
    fit$exposure
  } else {
    paste0(fit$exposure, " * ", fit$effect_modifier)
  }
  # outcome ~ exposure (* modifier) + confounders + strata(set). reformulate()
  # preserves the confounder transforms / interactions.
  model_formula <- stats::reformulate(
    termlabels = c(exposure_term, conf_terms, strata_term),
    response = fit$outcome
  )
  # Fit on a copy whose modifier is coerced to a factor with its unused levels
  # dropped: per-level odds ratios need discrete levels (and the model's
  # `xlevels`), and an empty/unused factor level would otherwise contribute an
  # all-zero interaction column aliased to NA, which would wrongly mark the
  # whole stratum-specific OR unestimable. droplevels() keeps the order of the
  # remaining declared levels, so a user-set reference level is preserved.
  fit_data <- fit$data
  em <- fit$effect_modifier
  if (!is.null(em)) {
    fit_data[[em]] <- droplevels(as.factor(fit_data[[em]]))
  }
  # `clogit` rewrites its own call to an unqualified `coxph(Surv(...) ~ ... +
  # strata(...))` and evaluates it in this frame, so those three survival names
  # must resolve here (the imports below) even though the entry point is called
  # qualified.
  model <- survival::clogit(model_formula, data = fit_data)

  # clogit's default na.action silently drops rows with a missing outcome,
  # exposure, or confounder. `model$n` is the rows actually used (coxph's
  # nobs() counts events, not rows), so the dropped count is taken against it.
  n_dropped <- nrow(fit$data) - model$n
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The odds ratio is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  model
}

#' Assemble the conditional odds-ratio contrast from a clogit fit
#'
#' Turns a fitted matched case-control conditional logistic regression into a
#' `matchatr_result` reporting the exposure's conditional odds ratio with a Wald
#' confidence interval. Reuses the shared `conditional_or_result()` assembly
#' (the same OR layer as the unmatched logistic engine).
#'
#' @details
#' The variance is the inverse partial-likelihood information matrix that
#' `survival::clogit` returns; the Wald interval is on the log-odds scale and
#' exponentiated. The risk difference and risk ratio are rejected as unidentified
#' from a case-control sample (no source-population prevalence q0), shared with
#' the unmatched engine. The conditional fit needs no sandwich; a cluster-robust
#' variance (relevant only when controls are reused, as in a nested case-control
#' sample) is not offered here, so `ci_method = "sandwich"` / `"bootstrap"` abort
#' with `matchatr_unsupported_variance`.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted `survival::clogit`.
#' @param type Character contrast scale; only `"or"` is computed, while
#'   `"difference"` / `"ratio"` abort with `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source; only `"model"` (the
#'   partial-likelihood information matrix) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the conditional odds ratio(s) for the
#'   exposure term.
#' @family estimators
#' @seealso [contrast()], `fit_clogit()`, `conditional_or_result()`
#' @noRd
contrast_clogit <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  # The conditional fit reports the partial-likelihood information-matrix
  # interval. A cluster-robust sandwich (for reused controls) and the bootstrap
  # belong to the risk-set / inclusion-weight designs, not the matched CMLE.
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the conditional logistic estimator."
        ),
        i = "It reports the partial-likelihood information-matrix interval; use `ci_method = \"model\"` (the default)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # With an effect modifier the contrast is one exposure OR per modifier level
  # (the stratum-specific OR), assembled from the joint partial-likelihood
  # variance; otherwise it is the single conditional OR shared with the
  # unmatched logistic engine.
  if (!is.null(fit$effect_modifier)) {
    return(stratum_specific_or_result(
      fit,
      model = fit$model,
      conf_level = conf_level,
      ci_method = ci_method,
      # coxph's nobs() counts events; the analysis size is the rows used (`$n`).
      n = fit$model$n,
      call = call
    ))
  }

  conditional_or_result(
    fit,
    model = fit$model,
    robust = FALSE,
    ci_method = ci_method,
    conf_level = conf_level,
    estimand = "conditional OR",
    # coxph's nobs() counts events; the analysis size is the rows used (`$n`).
    n = fit$model$n,
    call = call
  )
}
