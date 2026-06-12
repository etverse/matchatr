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
#'   interval and the result records `ci_method = "sandwich"`. `"bootstrap"`
#'   resamples within the case / control strata (the design-preserving
#'   `ccw_bootstrap_ci()`) and reports the percentile interval.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @param n_boot Integer number of bootstrap replicates when
#'   `ci_method = "bootstrap"`.
#' @returns A `matchatr_result` carrying the intervention means and the marginal
#'   contrast with causatr's sandwich variance (or the within-stratum bootstrap
#'   percentile interval).
#' @family estimators
#' @seealso [contrast()], `fit_ccw()`, `ccw_bootstrap_ci()`
#' @noRd
contrast_ccw <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env(),
  n_boot = 1000L
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

  contrasts <- data.table::as.data.table(res$contrasts)
  # A marginal g-formula contrast has no information-matrix variance distinct from
  # the influence-function / sandwich one causatr computes, so `"model"` and
  # `"sandwich"` both yield it; record causatr's `"sandwich"` rather than the
  # requested label. `"bootstrap"` overrides the interval with the within-stratum
  # percentile bootstrap (the point estimate stays the sandwich plug-in).
  # An estimated q0 (the design's `prevalence_n` set) widens the interval: the
  # bootstrap redraws q0* per replicate, the analytic path adds the delta-method
  # term to the sandwich SE.
  recorded_ci <- res$ci_method
  if (identical(ci_method, "bootstrap")) {
    boot <- ccw_bootstrap_ci(fit, type, conf_level, n_boot)
    contrasts$se <- boot$se
    contrasts$ci_lower <- boot$lower
    contrasts$ci_upper <- boot$upper
    recorded_ci <- "bootstrap"
  } else if (!is.null(fit$design$prevalence_n)) {
    adj <- ccw_apply_estimated_q0(
      fit,
      type,
      conf_level,
      contrasts$estimate,
      contrasts$ci_lower,
      contrasts$ci_upper
    )
    contrasts$se <- adj$se
    contrasts$ci_lower <- adj$lower
    contrasts$ci_upper <- adj$upper
  }

  new_matchatr_result(
    estimates = data.table::as.data.table(res$estimates),
    contrasts = contrasts,
    type = type,
    estimand = estimand,
    ci_method = recorded_ci,
    reference = "control",
    n = res$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = res$vcov,
    call = call
  )
}
