#' Fit the Samuelsen IPW-weighted additive-hazards model for NCC data
#'
#' Fits a weighted additive-hazards model with time-constant covariate effects on
#' a nested case-control sample drawn with [sample_ncc()] (`incl_prob = TRUE`).
#' The matching is broken: each unique cohort subject (deduplicated by
#' `.cohort_row`) enters once, upweighted by its Samuelsen inclusion weight 1/π_j
#' (cases at weight 1), and the additive model is fit by the weighted Lin-Ying
#' estimator (`lin_ying_additive()`) with a robust sandwich variance.
#'
#' @details
#' The additive-hazards model writes the hazard as a sum rather than a product,
#'
#'   λ(t | x) = λ₀(t) + γᵀ x,
#'
#' so each γ_j is the **excess hazard** (additive rate difference) for covariate
#' j — the change in the event rate per unit of x, the additive analogue of the
#' Cox hazard ratio. Unlike the hazard ratio it is a difference, so it can be
#' negative (a protective covariate lowers the rate). Effects are time-constant
#' (Lin & Ying 1994), giving one excess hazard per covariate rather than a
#' time-varying cumulative regression function. Estimation reuses controls via
#' the inverse-probability weights, and the validity under nested case-control
#' sampling follows Borgan & Langholz (1997, Biometrics 53(2)).
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"ipw_aalen"`, carrying
#'   the NCC analysis `data` with `.cohort_row` / `ipw_weight` columns plus the
#'   `outcome` / `exposure` / `confounders` slots and the design `time`.
#' @returns A `matchatr_aalen` list with the fitted constant excess hazards
#'   (`gamma`), their robust variance (`robvar`), the exposure term names
#'   (`exposure_terms`), the factor reference level (`reference`, or `NULL`), and
#'   the analysis sample size (`n`).
#' @family estimators
#' @seealso [matcha()], [contrast()], [sample_ncc()], `lin_ying_additive()`
#' @noRd
fit_ipw_aalen <- function(fit) {
  time_col <- require_ipw_ncc_columns(fit, "ipw_aalen")

  # Break the matching: the deduplicated, case-weight-forced analysis sample.
  dt <- ncc_ipw_analysis_data(fit)

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }

  # Drop rows with a missing time / outcome / covariate up front (the estimator
  # needs complete records), keeping the data.table aligned with the design.
  rhs <- stats::reformulate(c(fit$exposure, conf_terms))
  needed <- c(time_col, fit$outcome, all.vars(rhs))
  complete <- stats::complete.cases(
    do.call(cbind.data.frame, lapply(needed, function(v) dt[[v]]))
  )
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the IPW additive fit."
        ),
        i = "The excess hazard is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  dt <- dt[complete, , drop = FALSE]

  # Build the covariate design (no intercept — the additive baseline λ₀(t) is the
  # intercept). model.matrix expands a factor exposure / confounder to dummy
  # columns and records each column's originating term in `assign`, which locates
  # the exposure columns without re-deriving names. The data are complete, so the
  # design rows align with `dt` directly.
  mm <- stats::model.matrix(rhs, data = dt)
  assign <- attr(mm, "assign")
  Z <- mm[, assign > 0L, drop = FALSE] # drop the intercept column

  # The exposure is term 1 in `rhs`; its column(s) carry assign == 1.
  exposure_terms <- colnames(mm)[assign == 1L]

  ly <- lin_ying_additive(
    time = dt[[time_col]],
    status = as.integer(dt[[fit$outcome]]),
    Z = Z,
    w = dt[["ipw_weight"]]
  )

  exposure_col <- fit$data[[fit$exposure]]
  reference <- if (is.factor(exposure_col)) {
    levels(droplevels(exposure_col))[1L]
  } else {
    NULL
  }

  structure(
    list(
      gamma = ly$gamma,
      robvar = ly$robvar,
      exposure_terms = exposure_terms,
      reference = reference,
      n = nrow(mm)
    ),
    class = "matchatr_aalen"
  )
}

#' Assemble the IPW NCC additive excess-hazard contrast
#'
#' Turns a fitted Samuelsen IPW weighted additive-hazards model into a
#' `matchatr_result` reporting the exposure's **excess hazard** γ (additive rate
#' difference) with a robust Wald confidence interval. Only `type = "excess"` is
#' identified; the odds ratio, hazard ratio, time ratio, risk difference, and
#' risk ratio are rejected as off-scale.
#'
#' @details
#' The excess hazard is reported on its natural (linear) scale — it is a rate
#' difference, so the Wald interval γ ± z·SE(γ) is symmetric and the estimate may
#' be negative, unlike the exponentiated odds / hazard / time ratios. The
#' variance is the robust sandwich `lin_ying_additive()` returns, so
#' `ci_method = "model"` reads it directly; `"sandwich"` and `"bootstrap"` are
#' rejected (the robust sandwich is already the appropriate variance and
#' resampling does not account for the NCC weight estimation).
#'
#' @param fit A `matchatr_fit` whose `model` is a `matchatr_aalen` from
#'   `fit_ipw_aalen()`.
#' @param type Character contrast scale; only `"excess"` (excess hazard / rate
#'   difference) is computed.
#' @param ci_method Character variance source; only `"model"` (the robust
#'   sandwich) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the excess hazard(s) for the exposure
#'   term with the robust variance.
#' @family estimators
#' @seealso [contrast()], `fit_ipw_aalen()`
#' @noRd
contrast_ipw_aalen <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  if (type %in% c("or", "hr", "af")) {
    rlang::abort(
      c(
        "The IPW NCC additive-hazards model reports an excess hazard (rate difference), not a ratio.",
        i = 'Use `type = "excess"` (the default for `ipw_aalen`).'
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
          "\"` is not available for the IPW additive estimator."
        ),
        i = paste0(
          "It reports the robust sandwich of the weighted Lin-Ying estimator. ",
          "Use `ci_method = \"model\"` (the default)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  additive_excess_result(
    fit,
    conf_level = conf_level,
    ci_method = ci_method,
    call = call
  )
}

#' Assemble a linear-scale excess-hazard result from a fitted additive model
#'
#' The additive analogue of `conditional_or_result()` for the `ipw_aalen` engine.
#' Unlike the exponentiated ratios, the excess hazard is a rate difference on the
#' linear scale, so the Wald interval is symmetric (γ ± z·SE) and is **not**
#' exponentiated. The constant effects and their robust variance live on the
#' `matchatr_aalen` fit (`gamma`, `robvar`), and the exposure's coefficient(s)
#' are the `exposure_terms` recorded there.
#'
#' @param fit A `matchatr_fit` carrying a `matchatr_aalen` `model` and the
#'   `estimator` / `engine` labels recorded on the result.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param ci_method Character variance source recorded on the result.
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` whose `estimates` and `contrasts` carry the
#'   excess hazard(s) on the linear scale with a symmetric Wald interval.
#' @family estimators
#' @seealso `contrast_ipw_aalen()`, `conditional_or_result()`
#' @noRd
additive_excess_result <- function(
  fit,
  conf_level,
  ci_method,
  call = rlang::caller_env()
) {
  model <- fit$model
  term_names <- model$exposure_terms
  idx <- match(term_names, names(model$gamma))
  if (anyNA(idx)) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` has no estimable excess hazard."),
        i = paste0(
          "It is constant or collinear with the confounders, so its additive ",
          "effect is not identified."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  g <- unname(model$gamma[idx])
  vcov_exp <- model$robvar[idx, idx, drop = FALSE]
  dimnames(vcov_exp) <- list(term_names, term_names)
  s <- unname(sqrt(diag(vcov_exp)))
  lower <- g - z * s
  upper <- g + z * s

  # Excess hazard is a rate difference: the estimate, the CI, and the SE are all
  # on the same linear scale, so `estimates` and `contrasts` coincide (no
  # exp-transform separates a log-scale estimate from an exponentiated contrast).
  effect <- data.table::data.table(
    term = term_names,
    estimate = g,
    se = s,
    ci_lower = lower,
    ci_upper = upper
  )

  new_matchatr_result(
    estimates = effect,
    contrasts = data.table::data.table(
      comparison = term_names,
      estimate = g,
      se = s,
      ci_lower = lower,
      ci_upper = upper
    ),
    type = "excess",
    estimand = "excess hazard",
    ci_method = ci_method,
    reference = model$reference,
    n = model$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_exp,
    call = call
  )
}
