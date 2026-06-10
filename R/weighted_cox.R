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
