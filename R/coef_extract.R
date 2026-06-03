#' Variance-covariance of the estimable coefficients, by position
#'
#' Returns the model information-matrix or Huber-White sandwich variance over the
#' *estimable* coefficients (those [stats::coef()] does not set to `NA` for a
#' rank-deficient / aliased fit), together with their positions in the full
#' coefficient vector. `stats::vcov()` keeps aliased rows (carrying `NA`) while
#' `sandwich::sandwich()` drops them; both are reduced to the estimable set here,
#' whose rows correspond — *in order* — to the returned `est_pos`. Callers index
#' by position rather than by name, because `glm` permits non-unique coefficient
#' names (a `factor x level` concatenation can collide, e.g. `ses`+`low` equals
#' `se`+`slow`), which makes name indexing ambiguous.
#'
#' @param model A fitted model (e.g. `glm`, `mgcv::gam`).
#' @param robust Logical; use the Huber-White sandwich
#'   ([sandwich::sandwich()]) instead of the model information matrix.
#' @returns A list with `est_pos` (integer positions of the estimable
#'   coefficients in `coef(model)`) and `vcov` (their variance matrix, rows
#'   aligned to `est_pos`).
#' @family estimators
#' @noRd
estimable_vcov <- function(model, robust = FALSE) {
  beta <- stats::coef(model)
  est_pos <- which(!is.na(beta))
  vcov_mat <- if (isTRUE(robust)) {
    # sandwich() already excludes aliased coefficients, in coefficient order.
    sandwich::sandwich(model)
  } else {
    stats::vcov(model)[est_pos, est_pos, drop = FALSE]
  }
  list(est_pos = est_pos, vcov = vcov_mat)
}

#' Parametric term -> coefficient assignment, across fitters
#'
#' Maps the parametric model terms to the positions of the coefficients they
#' contribute, using term *position* (collision-free) rather than reconstructed
#' names. A `glm`/`lm` exposes this via the `model.matrix` `assign` attribute and
#' `terms()`; an `mgcv::gam` keeps its parametric coefficients first and stores
#' their term map in `model$assign` against `model$pterms` (the smooth-basis
#' coefficients have no parametric term and are excluded).
#'
#' @param model A fitted model (e.g. `glm`, `mgcv::gam`).
#' @returns A list with `assign` (integer term index per parametric coefficient,
#'   `0` for the intercept) and `labels` (the parametric `term.labels`). The
#'   `assign` vector aligns with the leading coefficients of `coef(model)`.
#' @family estimators
#' @noRd
term_assign <- function(model) {
  if (!is.null(model$assign) && !is.null(model$pterms)) {
    # mgcv::gam: parametric coefficients lead coef(model); $assign maps them to
    # $pterms, and the trailing smooth-basis coefficients are left unmapped.
    list(
      assign = model$assign,
      labels = attr(model$pterms, "term.labels")
    )
  } else {
    list(
      assign = attr(stats::model.matrix(model), "assign"),
      labels = attr(stats::terms(model), "term.labels")
    )
  }
}

#' Positions of the parametric (fixed-effect) coefficients
#'
#' The coefficient positions a tidy / summary table should report. An
#' `mgcv::gam` keeps its `nsdf` parametric coefficients first, followed by
#' smooth-basis coefficients (`s(age).1`, ...) that are penalized basis weights,
#' not odds ratios; those are excluded. Any other fit (e.g. `glm`) is fully
#' parametric.
#'
#' @param model A fitted model (e.g. `glm`, `mgcv::gam`).
#' @returns An integer vector of positions in `coef(model)`.
#' @family estimators
#' @noRd
parametric_positions <- function(model) {
  if (!is.null(model$nsdf)) {
    seq_len(model$nsdf)
  } else {
    seq_along(stats::coef(model))
  }
}

#' Assemble a conditional odds-ratio result from a fitted model
#'
#' Shared odds-ratio assembly for the conditional-likelihood OR engines — the
#' unmatched-CC logistic `glm` and the matched-CC conditional
#' `survival::clogit`. It locates the exposure coefficient(s) by term position
#' (`exposure_coef_index()`), reads their variance over the estimable
#' coefficients (`estimable_vcov()`), forms the Wald interval on the log-odds
#' scale, and exponentiates to the OR scale.
#'
#' @details
#' The interval is symmetric on the log scale and therefore asymmetric on the OR
#' scale, so `estimate +/- z * se` does not reproduce `ci_lower` / `ci_upper`;
#' the OR-scale `se` is the delta-method value OR * SE(log OR), kept for
#' reference and downstream composition, not to reconstruct the interval. The
#' reconstructable log-scale estimate and SE live in the result's `estimates`.
#' A constant or collinear exposure is aliased to `NA` by the fitter and aborts
#' with `matchatr_unestimable_exposure` rather than returning a silent `NA` OR.
#'
#' @param fit A `matchatr_fit` supplying the `exposure` / `data` columns and the
#'   `estimator` / `engine` labels recorded on the result.
#' @param model The fitted model (a `glm` or a `survival::clogit`) carrying the
#'   coefficients and their variance.
#' @param robust Logical; use the Huber-White sandwich
#'   ([sandwich::sandwich()]) instead of the model information matrix.
#' @param ci_method Character variance source recorded on the result.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param estimand Character estimand label (e.g. `"conditional OR"`).
#' @param n Integer analysis sample size to record (the fitter's complete-case
#'   count; `stats::nobs()` for a `glm`, `model$n` for a `clogit`, whose
#'   `nobs()` counts events rather than rows).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the log-OR estimates and the OR-scale
#'   contrasts for the exposure term(s).
#' @family estimators
#' @seealso `contrast_logistic()`, `contrast_clogit()`
#' @noRd
conditional_or_result <- function(
  fit,
  model,
  robust,
  ci_method,
  conf_level,
  estimand,
  n,
  call = rlang::caller_env()
) {
  beta <- stats::coef(model)
  # Two-sided Wald critical value for the requested confidence level.
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  idx <- exposure_coef_index(model, fit$exposure, call = call)
  term_labels <- names(beta)[idx]
  b <- unname(beta[idx])
  # A constant or collinear exposure is aliased to NA by the fitter: it has no
  # estimable coefficient, so its odds ratio is not identified. Refuse rather
  # than return a silent NA (mirrors the degenerate-outcome rejection in
  # resolve_binary_outcome()).
  if (anyNA(b)) {
    rlang::abort(
      c(
        paste0("Exposure `", fit$exposure, "` has no estimable effect."),
        i = paste0(
          "It is constant or collinear with the confounders / strata, so ",
          "its odds ratio is not identified."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }

  # Index the variance by POSITION, not coefficient name: names can collide
  # (factor x level concatenation) and the model vcov keeps aliased rows while
  # the sandwich drops them. `estimable_vcov()` returns the variance over the
  # estimable coefficients with their positions (`est_pos`); the exposure
  # coefficients are non-aliased (guard above), so each maps into that set.
  ev <- estimable_vcov(model, robust = robust)
  sel <- match(idx, ev$est_pos)
  vcov_exp <- ev$vcov[sel, sel, drop = FALSE]
  dimnames(vcov_exp) <- list(term_labels, term_labels)
  s <- unname(sqrt(diag(vcov_exp)))
  log_lower <- b - z * s
  log_upper <- b + z * s

  # For a factor exposure, each contrast is a level versus the factor's
  # reference (baseline) level; record it so the OR rows are interpretable. Read
  # it from the model's `xlevels` (the levels actually used in fitting, with
  # unused declared levels dropped), not from the raw column whose first declared
  # level may never occur. Fall back to the present levels if `xlevels` is absent.
  exposure_col <- fit$data[[fit$exposure]]
  reference <- if (is.factor(exposure_col)) {
    xl <- model$xlevels[[fit$exposure]]
    if (is.null(xl)) levels(droplevels(exposure_col))[1] else xl[1]
  } else {
    NULL
  }

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
    estimand = estimand,
    ci_method = ci_method,
    reference = reference,
    n = n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_exp,
    call = call
  )
}

#' Locate the coefficient(s) belonging to the exposure term
#'
#' Maps the exposure to its coefficient position(s) via the parametric term
#' assignment (`term_assign()`), i.e. by term *position*, so it is collision-free
#' even when `glm` produces non-unique coefficient names from `factor x level`
#' concatenation (e.g. exposure `ses`+`low` shares the name `"seslow"` with a
#' confounder `se`+`slow`). The intercept is never returned, and the
#' approach works uniformly for a binary, continuous, or (unordered) factor
#' exposure across `glm` and `mgcv::gam`.
#'
#' @param model A fitted model (e.g. `glm`, `mgcv::gam`).
#' @param exposure Character scalar exposure column name; it must enter the
#'   model as a parametric main-effect term.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of coefficient positions (in `coef(model)`) for
#'   the exposure term; aborts with `matchatr_bad_input` if the exposure is not
#'   a parametric main-effect term.
#' @family estimators
#' @noRd
exposure_coef_index <- function(model, exposure, call = rlang::caller_env()) {
  ta <- term_assign(model)
  pos <- match(exposure, ta$labels)
  if (is.na(pos)) {
    rlang::abort(
      paste0(
        "Exposure `",
        exposure,
        "` is not a parametric main-effect term in the fitted model."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  which(ta$assign == pos)
}
