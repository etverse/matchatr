#' Fit an unmatched case-control logistic regression
#'
#' Fits a binomial-family regression to estimate the conditional (adjusted) log
#' odds ratios from an unmatched case-control sample. The fitted model is
#' `outcome ~ exposure + confounders`; the confounder formula's terms, including
#' transforms (`poly(age, 2)`), smooths (`s(age)` under a GAM fitter), and
#' interactions (`age:smoke`), are carried through verbatim.
#'
#' @details
#' Only the slope coefficients carry over from the cohort logistic model: under
#' separate case / control sampling the intercept is offset by
#' log(case sampling fraction / control sampling fraction) (Prentice & Pyke,
#' 1979), so it is never a baseline risk. The outcome is modelled on its native
#' encoding — a numeric 0/1 column, a logical (`TRUE` = case), or a two-level
#' factor (second level = case) all reproduce the 0/1 coding of
#' `resolve_binary_outcome()`. Rows with a missing outcome, exposure, or
#' confounder are dropped by the default `na.action`; a `matchatr_dropped_rows`
#' warning reports how many.
#'
#' The fitter is pluggable via `details$model_fn` (default [stats::glm()]); any
#' function with a `(formula, family, data)` interface works, e.g.
#' `mgcv::gam` for a smooth confounder term such as `s(age)`. The exposure is
#' entered parametrically so its odds ratio stays a single (or per-level) number
#' regardless of the fitter.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"glm_logistic"`,
#'   carrying the analysis `data`, the `outcome` / `exposure` column names, the
#'   `confounders` formula (or `NULL`), and `details$model_fn`.
#' @returns The fitted binomial model object (a `glm`, or whatever `model_fn`
#'   returns).
#' @family estimators
#' @seealso [matcha()], [contrast()]
#' @noRd
fit_logistic_cc <- function(fit) {
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # outcome ~ exposure + confounders. reformulate() preserves transforms,
  # smooths, and interactions from the confounder formula; the exposure enters
  # as a main effect so its coefficient(s) identify the conditional log OR.
  model_formula <- stats::reformulate(
    termlabels = c(fit$exposure, conf_terms),
    response = fit$outcome
  )
  # Pluggable fitter (default stats::glm); a GAM fitter allows smooth confounder
  # adjustment while the parametric exposure keeps an interpretable OR.
  model_fn <- if (is.null(fit$details$model_fn)) {
    stats::glm
  } else {
    fit$details$model_fn
  }
  model <- model_fn(
    model_formula,
    family = stats::binomial(),
    data = fit$data
  )
  # A custom model_fn may ignore `family` and return a non-binomial fit (e.g. an
  # OLS lm), whose exponentiated slope is not an odds ratio. Verify the fit is
  # binomial before any OR is reported. Attributed to the user's matcha() call.
  fitted_family <- tryCatch(
    stats::family(model)$family,
    error = function(e) NA_character_
  )
  if (!identical(fitted_family, "binomial")) {
    rlang::abort(
      c(
        "`model_fn` did not return a binomial fit.",
        i = "It must fit a binomial-family logistic model (e.g. `stats::glm`, `mgcv::gam`)."
      ),
      class = c("matchatr_bad_model_fit", "matchatr_error"),
      call = fit$call
    )
  }
  # glm's default na.action silently drops rows with a missing outcome,
  # exposure, or confounder. Surface that listwise deletion so a missing-data
  # problem is not mistaken for the full sample.
  n_dropped <- nrow(fit$data) - stats::nobs(model)
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

#' Assemble the conditional odds-ratio contrast from a logistic fit
#'
#' Turns a fitted unmatched case-control logistic regression into a
#' `matchatr_result` reporting the exposure's conditional
#' odds ratio with a Wald confidence interval. The risk difference and risk
#' ratio are rejected: under separate case / control sampling the marginal
#' outcome frequency is fixed by design, so absolute risks (and hence RD / RR)
#' are not identified without the source-population prevalence q0 — only the
#' (non-collapsible) conditional OR is.
#'
#' @details
#' The interval is Wald on the log-odds scale, exponentiated to the OR scale:
#' OR = exp(b), CI = exp(b +/- z * SE(b)), with the standard error from the
#' model information matrix (`ci_method = "model"`) or the Huber-White sandwich
#' (`ci_method = "sandwich"`). The intercept is excluded — it is not an
#' interpretable baseline risk in a case-control sample.
#'
#' The reported OR-scale SE is the delta-method value OR * SE(b) (matching the
#' `causatr` result convention). Because the interval is symmetric on the log
#' scale, it is asymmetric on the OR scale, so `estimate +/- z * se` does NOT
#' reproduce `ci_lower` / `ci_upper`; use the reported bounds directly. The SE
#' is provided for reference and downstream composition, not for reconstructing
#' the interval.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted binomial `glm`.
#' Reject an unidentified risk-difference / risk-ratio request
#'
#' Shared by the classical case-control odds-ratio engines (logistic and
#' Mantel-Haenszel): a marginal risk difference or risk ratio is not identified
#' from a case-control sample whose marginal outcome frequency is fixed by
#' design, so `type = "difference"` / `"ratio"` aborts, pointing to the
#' conditional OR or to a case-control-weighted estimator with a prevalence q0.
#'
#' @param type Character contrast scale requested by the caller.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly for `type = "or"`; otherwise aborts with class
#'   `matchatr_unidentified_estimand`.
#' @family estimators
#' @noRd
reject_unidentified_rd_rr <- function(type, call = rlang::caller_env()) {
  if (!type %in% c("difference", "ratio")) {
    return(invisible(NULL))
  }
  estimand_name <- if (identical(type, "difference")) {
    "risk difference"
  } else {
    "risk ratio"
  }
  rlang::abort(
    c(
      paste0(
        "The ",
        estimand_name,
        " is not identified from an unmatched case-control sample without ",
        "the source-population prevalence q0."
      ),
      i = "Report the conditional odds ratio with `type = \"or\"`.",
      i = paste0(
        "For a marginal risk difference / ratio, supply `prevalence =` on ",
        "the design and use a case-control-weighted estimator (e.g. ",
        "`estimator = \"ccw_gformula\"`)."
      )
    ),
    class = c("matchatr_unidentified_estimand", "matchatr_error"),
    call = call
  )
}

#' @param type Character contrast scale; `"or"` is computed, while
#'   `"difference"` / `"ratio"` abort with `matchatr_unidentified_estimand`.
#' @param ci_method Character variance source: `"model"` (information matrix) or
#'   `"sandwich"` (robust). `"bootstrap"` aborts with
#'   `matchatr_unsupported_variance` — the conditional OR is reported with a
#'   Wald interval.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the log-OR estimates and the OR-scale
#'   contrast for the exposure term(s).
#' @family estimators
#' @seealso [contrast()], `fit_logistic_cc()`
#' @noRd
contrast_logistic <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  # Bootstrap resampling of an unmatched case-control sample must respect the
  # fixed case / control counts; a single conditional OR is reported with the
  # exact Wald interval instead.
  if (identical(ci_method, "bootstrap")) {
    rlang::abort(
      c(
        "Bootstrap confidence intervals are not provided for the conditional odds ratio.",
        i = "Use `ci_method = \"model\"` (Wald) or `ci_method = \"sandwich\"` (robust)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # Shared conditional-OR assembly: locate the exposure coefficient(s), form the
  # Wald interval on the log-odds scale, exponentiate to the OR scale. The
  # logistic engine offers the model information matrix or the Huber-White
  # sandwich. The analysis n is the number of rows glm actually used (complete
  # cases), not the full sample, which may differ when confounders carry NAs.
  conditional_or_result(
    fit,
    model = fit$model,
    robust = identical(ci_method, "sandwich"),
    ci_method = ci_method,
    conf_level = conf_level,
    estimand = "conditional OR",
    n = stats::nobs(fit$model),
    call = call
  )
}
