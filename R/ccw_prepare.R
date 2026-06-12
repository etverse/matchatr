#' Shared preprocessing for the case-control-weighted estimators
#'
#' Validates that an adjustment set is supplied, coerces the outcome and exposure
#' to 0/1, drops rows with missing outcome / exposure / confounders (listwise
#' deletion, warning `matchatr_dropped_rows`), and builds the Rose & van der Laan
#' case-control weight vector on the complete-case sample — the common front end
#' for every CCW engine (the causatr-delegated g-formula / IPW / AIPW and the
#' hand-rolled TMLE), which differ only in what they do with the weighted,
#' complete-case, 0/1-coded sample.
#'
#' @param fit A `matchatr_fit` whose `engine` is a `ccw_*` estimator, carrying the
#'   case-control `data`, the `outcome` / `exposure` / `confounders` roles, and a
#'   `design` whose `prevalence` (q0) is set.
#' @returns A list with `dt` (a complete-case `data.table` copy of `data` with the
#'   outcome and exposure recoded to 0/1) and `weights` (the numeric case-control
#'   weights, one per row of `dt`, computed on the complete-case outcome so the
#'   weighted case fraction equals q0). Warns `matchatr_dropped_rows` when rows are
#'   dropped; aborts with `matchatr_bad_input` when no `confounders` are supplied or
#'   the exposure is non-binary.
#' @family estimators
#' @seealso `fit_ccw()`, `fit_ccw_tmle()`, `cc_weights()`
#' @noRd
ccw_prepare <- function(fit) {
  # Every CCW estimator needs an adjustment set: g-formula standardizes a
  # confounder-adjusted outcome model, IPW needs a propensity model, AIPW and
  # TMLE both. With no confounders there is nothing to adjust for; reject early.
  if (is.null(fit$confounders)) {
    rlang::abort(
      c(
        "The case-control-weighted estimators require `confounders` for the adjustment model(s).",
        i = paste0(
          "Supply an adjustment set, e.g. `confounders = ~ age + smoke`, on ",
          "`matcha()`."
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }

  # Coerce both roles to 0/1: the outcome so the weighted GLM reads a proper 0/1
  # response, the exposure so the treat-all / treat-none interventions match the
  # treatment coding. A non-binary exposure has no binary average-treatment-effect
  # contrast and is rejected here.
  y01 <- resolve_binary_outcome(fit$data, fit$outcome)
  x01 <- resolve_binary_exposure(
    fit$data,
    fit$exposure,
    estimator_label = ccw_estimator_label(fit$estimator),
    alternative = "a conditional estimator (e.g. estimator = \"logistic\")"
  )

  dt <- data.table::copy(data.table::as.data.table(fit$data))
  dt[[fit$outcome]] <- y01
  dt[[fit$exposure]] <- x01

  # Complete-case the analysis sample here, once, for the whole CCW family. The
  # hand-rolled TMLE engine cannot tolerate NA (its prediction / clever-covariate
  # vectors would misalign), and causatr's IPW / AIPW reject a confounder NA the
  # outcome mask does not cover, so doing the listwise deletion up front gives
  # every engine an aligned complete-case sample and one consistent behaviour.
  # The weights are computed AFTER dropping so the weighted case fraction still
  # equals q0 over the analysed (complete-case) sample. Multiple imputation is the
  # principled alternative for missing confounders (a later phase).
  conf_vars <- all.vars(fit$confounders)
  complete <- !is.na(y01) &
    !is.na(x01) &
    stats::complete.cases(as.data.frame(dt)[, conf_vars, drop = FALSE])
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the case-control-weighted fit."
        ),
        i = "The marginal effect is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
    dt <- dt[complete, , drop = FALSE]
    y01 <- y01[complete]
  }

  list(dt = dt, weights = cc_weights(fit$design$prevalence, y01))
}

#' Map a CCW estimator name to its causatr engine
#'
#' @param estimator Character scalar matchatr CCW estimator (`"ccw_gformula"`,
#'   `"ccw_ipw"`, `"ccw_aipw"`).
#' @param call Caller environment surfaced in the defensive error.
#' @returns The causatr `estimator` string (`"gcomp"`, `"ipw"`, `"aipw"`).
#' @family estimators
#' @noRd
ccw_causat_estimator <- function(estimator, call = rlang::caller_env()) {
  switch(
    estimator,
    ccw_gformula = "gcomp",
    ccw_ipw = "ipw",
    ccw_aipw = "aipw",
    # Defensive: the dispatch table only routes the three names above here.
    rlang::abort(
      paste0("Unknown case-control-weighted estimator `", estimator, "`."),
      class = c("matchatr_bad_estimator", "matchatr_error"),
      call = call
    )
  )
}

#' Human-readable label for a CCW estimator
#'
#' @param estimator Character scalar matchatr CCW estimator name.
#' @returns A short label used in error messages (e.g. `"CCW AIPW"`).
#' @family estimators
#' @noRd
ccw_estimator_label <- function(estimator) {
  switch(
    estimator,
    ccw_gformula = "CCW g-formula",
    ccw_ipw = "CCW IPW",
    ccw_aipw = "CCW AIPW",
    ccw_tmle = "CCW TMLE",
    "case-control-weighted"
  )
}
