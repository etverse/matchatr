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
