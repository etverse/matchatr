#' Fit an unmatched case-control logistic regression
#'
#' Wraps [stats::glm()] with a binomial family to estimate the conditional
#' (adjusted) log odds ratios from an unmatched case-control sample. The fitted
#' model is `outcome ~ exposure + confounders`; the confounder formula's terms,
#' including transforms (`poly(age, 2)`) and interactions (`age:smoke`), are
#' carried through verbatim.
#'
#' @details
#' Only the slope coefficients carry over from the cohort logistic model: under
#' separate case / control sampling the intercept is offset by
#' log(case sampling fraction / control sampling fraction) (Prentice & Pyke,
#' 1979), so it is never a baseline risk. The outcome is modelled on its native
#' encoding — a numeric 0/1 column, a logical (`TRUE` = case), or a two-level
#' factor (second level = case) all reproduce the 0/1 coding of
#' `resolve_binary_outcome()`.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"glm_logistic"`,
#'   carrying the analysis `data`, the `outcome` / `exposure` column names, and
#'   the `confounders` formula (or `NULL`).
#' @returns A `glm` object, the fitted binomial model.
#' @family estimators
#' @seealso [matcha()], [contrast()]
#' @noRd
fit_logistic_cc <- function(fit) {
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  # outcome ~ exposure + confounders. reformulate() preserves transforms and
  # interactions from the confounder formula; the exposure enters as a main
  # effect so its coefficient(s) identify the conditional log OR.
  model_formula <- stats::reformulate(
    termlabels = c(fit$exposure, conf_terms),
    response = fit$outcome
  )
  stats::glm(model_formula, family = stats::binomial(), data = fit$data)
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
#' (`ci_method = "sandwich"`). The reported OR-scale SE is the delta-method
#' value OR * SE(b). The intercept is excluded — it is not an interpretable
#' baseline risk in a case-control sample.
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

  idx <- exposure_coef_index(model, fit$exposure, call = call)
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
    reference = NULL,
    n = nrow(fit$data),
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
#' Maps the exposure column to its coefficient position(s) in a fitted model via
#' the model matrix `assign` attribute, so a binary, continuous, or (for later
#' multi-level support) factor exposure all resolve correctly without parsing
#' coefficient names. The intercept (`assign == 0`) is therefore never returned.
#'
#' @param model A fitted model (e.g. `glm`) with a terms object and a model
#'   matrix carrying an `assign` attribute.
#' @param exposure Character scalar exposure column name; it must enter the
#'   model as a main-effect term.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of coefficient indices for the exposure term;
#'   aborts with `matchatr_bad_input` if the exposure is not a main-effect term.
#' @family estimators
#' @noRd
exposure_coef_index <- function(model, exposure, call = rlang::caller_env()) {
  assign <- attr(stats::model.matrix(model), "assign")
  term_labels <- attr(stats::terms(model), "term.labels")
  pos <- match(exposure, term_labels)
  if (is.na(pos)) {
    rlang::abort(
      paste0(
        "Exposure `",
        exposure,
        "` is not a main-effect term in the fitted model."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  which(assign == pos)
}
