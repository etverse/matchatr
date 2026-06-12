#' Fit a case-control-weighted causal model (g-formula / IPW / AIPW)
#'
#' Implements the Rose & van der Laan case-control-weighted (CCW) family: the
#' case-control sample is reweighted to the source population with the
#' marginal-prevalence weights from `cc_weights()`, a cohort causal estimator is
#' fitted on the weighted sample, and the contrast step standardizes it to a
#' **marginal** causal effect (risk difference, risk ratio, or marginal odds
#' ratio). The estimator is chosen by `fit$estimator`: `"ccw_gformula"` fits a
#' weighted outcome model (g-computation), `"ccw_ipw"` a weighted propensity
#' (inverse-probability) model, and `"ccw_aipw"` the doubly-robust augmented IPW
#' (consistent if **either** the outcome or the propensity model is correct).
#' Point estimation and variance are delegated to [causatr::causat()] with the
#' case-control weights passed as observation weights; matchatr owns only the
#' weighting layer.
#'
#' @details
#' Under separate case / control sampling the marginal outcome frequency is
#' fixed by design, so a case-control fit identifies only the conditional odds
#' ratio. The case-control weights (Rose & van der Laan 2008) restore the source
#' population's outcome margin q0, so the weighted empirical distribution mimics
#' the cohort and a cohort estimator on it targets the marginal estimand
#' (Rose & van der Laan 2008, *Int. J. Biostat.* 4(1); the doubly-robust CCW-AIPW
#' is Rose & van der Laan 2014, *Biometrics* 70(1)).
#'
#' The outcome and exposure are coerced to 0/1 so the marginal contrast's
#' interventions (treat-all versus treat-none) align with the fitted treatment
#' coding; a non-binary exposure is rejected (the binary average treatment effect
#' only). The model is fitted with `family = "quasibinomial"`: the fractional
#' case-control weights make the Bernoulli "successes" non-integer, and
#' quasibinomial fits the same mean model as binomial without the spurious
#' `non-integer #successes` warning — and the dispersion it adds is unused by the
#' marginal contrast and causatr's sandwich. The IPW and AIPW estimators also fit
#' a propensity model, whose fitter (`propensity_model_fn`) is named explicitly so
#' causatr does not warn about defaulting it.
#'
#' @param fit A `matchatr_fit` whose `engine` is one of `"ccw_gformula"` /
#'   `"ccw_ipw"` / `"ccw_aipw"`, carrying the case-control `data`, the `outcome` /
#'   `exposure` / `confounders` roles, and a `design` whose `prevalence` (q0) is
#'   set. The `matchatr_missing_prevalence` guard in [matcha()] ensures q0 is
#'   present before this runs.
#' @returns A `causatr_fit` object (the weighted causal fit) stored in the
#'   `matchatr_fit`'s `model` slot; [contrast()] turns it into the marginal
#'   effect.
#' @family estimators
#' @seealso [matcha()], [contrast()], `cc_weights()`, [causatr::causat()]
#' @noRd
fit_ccw <- function(fit) {
  causat_estimator <- ccw_causat_estimator(fit$estimator)
  prep <- ccw_prepare(fit)

  # The outcome model fitter (NULL -> stats::glm); a user may pass mgcv::gam for
  # smooth confounder adjustment.
  model_fn <- fit$details$model_fn
  if (is.null(model_fn)) {
    model_fn <- stats::glm
  }

  # `family` governs the weighted outcome / marginal-mean fit (g-formula's outcome
  # model, AIPW's outcome model, IPW's weighted marginal mean) — quasibinomial
  # fits the same mean model as binomial but is silent on the fractional
  # "successes" the case-control weights produce. The propensity model's family is
  # auto-detected by causatr from the (binary) treatment, so this is the right
  # family for all three estimators.
  args <- list(
    data = prep$dt,
    outcome = fit$outcome,
    treatment = fit$exposure,
    confounders = fit$confounders,
    estimator = causat_estimator,
    family = "quasibinomial",
    weights = as.numeric(prep$weights),
    model_fn = model_fn
  )
  # IPW / AIPW fit a propensity model; name its fitter so causatr does not warn
  # about defaulting it (the outcome `model_fn` is unused by the pure IPW path).
  if (causat_estimator %in% c("ipw", "aipw")) {
    args$propensity_model_fn <- stats::glm
  }
  do.call(causatr::causat, args)
}

#' Shared preprocessing for the case-control-weighted estimators
#'
#' Validates that an adjustment set is supplied, coerces the outcome and exposure
#' to 0/1, and builds the Rose & van der Laan case-control weight vector — the
#' common front end for every CCW engine (the causatr-delegated g-formula / IPW /
#' AIPW and the hand-rolled TMLE), which differ only in what they do with the
#' weighted, 0/1-coded sample.
#'
#' @param fit A `matchatr_fit` whose `engine` is a `ccw_*` estimator, carrying the
#'   case-control `data`, the `outcome` / `exposure` / `confounders` roles, and a
#'   `design` whose `prevalence` (q0) is set.
#' @returns A list with `dt` (a `data.table` copy of `data` with the outcome and
#'   exposure recoded to 0/1) and `weights` (the numeric case-control weights, one
#'   per row of `dt`). Aborts with `matchatr_bad_input` when no `confounders` are
#'   supplied or the exposure is non-binary.
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

#' Assemble the case-control-weighted marginal contrast
#'
#' Turns a fitted CCW g-computation model into a `matchatr_result` reporting the
#' marginal causal effect on the requested scale: the risk difference
#' (`type = "difference"`), risk ratio (`type = "ratio"`), or marginal odds
#' ratio (`type = "or"`). The standardization and its variance are delegated to
#' [causatr::contrast()] on the weighted fit, comparing the treat-all
#' (`static(1)`) and treat-none (`static(0)`) interventions.
#'
#' @details
#' The two interventions are the static "everyone exposed" and "everyone
#' unexposed" regimes, so the contrast is the average treatment effect on the
#' chosen scale, standardized over the (weight-restored) population covariate
#' distribution. The conditional odds ratio, hazard ratio, and the IPW NCC
#' scales are not marginal effects and are rejected.
#'
#' causatr returns the influence-function / sandwich variance on the weighted
#' fit, so `ci_method = "model"` and `ci_method = "sandwich"` both forward to it.
#' A bootstrap interval must resample within the case and control strata and
#' recompute the q0 weights each replicate (a design-specific correction not yet
#' wired), so `ci_method = "bootstrap"` is rejected here.
#'
#' @param fit A `matchatr_fit` whose `model` is a `causatr_fit` fitted by
#'   `fit_ccw()`.
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"or"`.
#' @param ci_method Character variance source. A marginal g-formula contrast has
#'   no information-matrix variance distinct from the influence-function /
#'   sandwich one, so `"model"` and `"sandwich"` both yield causatr's sandwich
#'   interval and the result records `ci_method = "sandwich"`. `"bootstrap"` is
#'   rejected (it must resample within the case / control strata and recompute
#'   the q0 weights, a later addition).
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the intervention means and the marginal
#'   contrast with causatr's variance.
#' @family estimators
#' @seealso [contrast()], `fit_ccw()`
#' @noRd
contrast_ccw <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  # Case-control weighting targets a marginal effect on the RD / RR / marginal-OR
  # scales; the conditional OR, hazard ratio, and AFT / additive scales belong to
  # other engines and are not relabelings of this fit.
  if (!type %in% c("difference", "ratio", "or")) {
    rlang::abort(
      c(
        paste0(
          "A case-control-weighted estimator reports a marginal effect, not `type = \"",
          type,
          "\"`."
        ),
        i = paste0(
          'Use `type = "difference"` (risk difference), `"ratio"` (risk ',
          'ratio), or `"or"` (marginal odds ratio).'
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (identical(ci_method, "bootstrap")) {
    rlang::abort(
      c(
        '`ci_method = "bootstrap"` is not available for the case-control-weighted estimators.',
        i = paste0(
          'Use `ci_method = "model"` or `ci_method = "sandwich"` (causatr\'s ',
          "influence-function variance on the weighted fit)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # The treat-all / treat-none static regimes give the average treatment effect;
  # the reference is "everyone unexposed".
  interventions <- list(
    treated = causatr::static(1),
    control = causatr::static(0)
  )

  res <- causatr::contrast(
    fit$model,
    interventions = interventions,
    type = type,
    reference = "control",
    ci_method = "sandwich",
    conf_level = conf_level
  )

  estimand <- switch(
    type,
    difference = "marginal risk difference",
    ratio = "marginal risk ratio",
    or = "marginal odds ratio"
  )

  new_matchatr_result(
    estimates = data.table::as.data.table(res$estimates),
    contrasts = data.table::as.data.table(res$contrasts),
    type = type,
    estimand = estimand,
    # A marginal g-formula contrast has no information-matrix variance distinct
    # from the influence-function / sandwich one causatr computes, so `"model"`
    # and `"sandwich"` both yield it. Record what was actually used (causatr's
    # `"sandwich"`) rather than the requested label, so the result does not
    # claim a model-based interval it did not compute.
    ci_method = res$ci_method,
    reference = "control",
    n = res$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = res$vcov,
    call = call
  )
}
