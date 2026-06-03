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
  if (type %in% c("difference", "ratio")) {
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

  model <- fit$model
  beta <- stats::coef(model)
  # Two-sided Wald critical value for the requested confidence level.
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  idx <- exposure_coef_index(model, fit$exposure, fit$data, call = call)
  term_labels <- names(beta)[idx]
  b <- unname(beta[idx])
  # A constant or collinear exposure is aliased to NA by glm: it has no
  # estimable coefficient, so its odds ratio is not identified. Refuse rather
  # than return a silent NA (mirrors the degenerate-outcome rejection in
  # resolve_binary_outcome()).
  if (anyNA(b)) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` has no estimable effect."),
        i = paste0(
          "It is constant or collinear with the confounders / intercept, so ",
          "its odds ratio is not identified."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }

  # Index the variance by coefficient NAME, not position: the model vcov keeps
  # aliased rows while the sandwich drops them, so a positional index into
  # `diag(vcov)` would misalign once any earlier term is aliased.
  # `estimable_vcov()` reduces both sources to the estimable set, and the
  # (non-aliased) exposure terms are guaranteed present after the guard above.
  vcov_exp <- estimable_vcov(
    model,
    robust = identical(ci_method, "sandwich")
  )[term_labels, term_labels, drop = FALSE]
  s <- unname(sqrt(diag(vcov_exp)))
  log_lower <- b - z * s
  log_upper <- b + z * s

  # For a factor exposure, each contrast is a level versus the factor's
  # reference (baseline) level; record it so the OR rows are interpretable.
  exposure_col <- fit$data[[fit$exposure]]
  reference <- if (is.factor(exposure_col)) levels(exposure_col)[1] else NULL

  estimates <- data.table::data.table(
    term = term_labels,
    estimate = b, # log OR (raw coefficient)
    se = s,
    ci_lower = log_lower,
    ci_upper = log_upper
  )
  contrasts <- data.table::data.table(
    comparison = term_labels,
    estimate = exp(b), # OR
    se = exp(b) * s, # delta-method SE on the OR scale
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "or",
    estimand = "conditional OR",
    ci_method = ci_method,
    reference = reference,
    # The analysis n is the number of rows glm actually used (complete cases),
    # not the full sample, which may differ when confounders carry NAs.
    n = stats::nobs(model),
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_exp,
    call = call
  )
}

#' Variance-covariance of the estimable coefficients
#'
#' Returns the model information-matrix or Huber-White sandwich
#' variance-covariance matrix restricted to the *estimable* coefficients — those
#' [stats::coef()] does not set to `NA` for a rank-deficient / aliased fit —
#' indexed by coefficient name. `stats::vcov()` keeps aliased rows (carrying
#' `NA`) while `sandwich::sandwich()` drops them, so a positional index into one
#' would misalign against the other; reducing both to the same estimable set and
#' indexing by name avoids that.
#'
#' @param model A fitted model (e.g. `glm`).
#' @param robust Logical; use the Huber-White sandwich
#'   ([sandwich::sandwich()]) instead of the model information matrix.
#' @returns A named numeric matrix over the estimable coefficients.
#' @family estimators
#' @noRd
estimable_vcov <- function(model, robust = FALSE) {
  beta <- stats::coef(model)
  estimable <- names(beta)[!is.na(beta)]
  v <- if (isTRUE(robust)) {
    sandwich::sandwich(model)
  } else {
    stats::vcov(model)
  }
  v[estimable, estimable, drop = FALSE]
}

#' Locate the coefficient(s) belonging to the exposure term
#'
#' Maps the exposure column to its coefficient position(s) by NAME, derived from
#' the exposure's storage type rather than from a model-internal `assign`
#' attribute, so it works uniformly across fitters (a `glm` and an `mgcv::gam`
#' do not expose the same `model.matrix` `assign`). A numeric / logical exposure
#' contributes the single coefficient named after the column; an (unordered)
#' factor with levels `l1, ..., lk` contributes `paste0(exposure, l2..lk)` (the
#' reference level `l1` has no coefficient). A binary or continuous exposure is
#' the one-coefficient special case.
#'
#' An *ordered* factor is rejected: `glm`/`gam` fit it with polynomial contrasts
#' (`.L`, `.Q`, ...), whose coefficients are not per-level odds ratios. The
#' caller should pass a numeric score (for a trend OR) or an unordered factor
#' (for per-level ORs).
#'
#' @param model A fitted model (e.g. `glm`, `mgcv::gam`).
#' @param exposure Character scalar exposure column name.
#' @param data The analysis data, used to read the exposure's type / levels.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of coefficient indices for the exposure term;
#'   aborts with `matchatr_bad_input` if the exposure is an ordered factor or
#'   contributes no coefficient.
#' @family estimators
#' @noRd
exposure_coef_index <- function(
  model,
  exposure,
  data,
  call = rlang::caller_env()
) {
  exposure_col <- data[[exposure]]
  if (is.factor(exposure_col) && is.ordered(exposure_col)) {
    rlang::abort(
      c(
        paste0(
          "Exposure `",
          exposure,
          "` is an ordered factor, which is fit with polynomial contrasts ",
          "(not per-level odds ratios)."
        ),
        i = "Pass a numeric score for a trend OR, or an unordered factor for per-level ORs."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  coef_names <- names(stats::coef(model))
  # Coefficient name(s) the exposure column contributes: the bare column name
  # for numeric / logical, or `paste0(exposure, level)` for each factor level
  # (the reference level produces no coefficient, so it simply does not match).
  candidates <- if (is.factor(exposure_col)) {
    paste0(exposure, levels(exposure_col))
  } else {
    exposure
  }
  idx <- which(coef_names %in% candidates)
  if (length(idx) == 0L) {
    rlang::abort(
      paste0(
        "Exposure `",
        exposure,
        "` contributes no coefficient to the fitted model."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  idx
}
