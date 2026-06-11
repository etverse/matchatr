#' Fit a case-control-weighted g-computation model
#'
#' Implements the Rose & van der Laan case-control-weighted (CCW) g-formula:
#' the case-control sample is reweighted to the source population with the
#' marginal-prevalence weights from `cc_weights()`, a weighted outcome model is
#' fitted, and the contrast step standardizes it to a **marginal** causal effect
#' (risk difference, risk ratio, or marginal odds ratio). Point estimation and
#' variance are delegated to [causatr::causat()] with the case-control weights
#' passed as observation weights; matchatr owns only the weighting layer.
#'
#' @details
#' Under separate case / control sampling the marginal outcome frequency is
#' fixed by design, so a case-control fit identifies only the conditional odds
#' ratio. The case-control weights (Rose & van der Laan 2008) restore the source
#' population's outcome margin q0, so the weighted empirical distribution mimics
#' the cohort and a cohort g-computation on it targets the marginal estimand
#' (Rose & van der Laan 2008, *Int. J. Biostat.* 4(1)).
#'
#' The outcome and exposure are coerced to 0/1 so the marginal contrast's
#' interventions (treat-all versus treat-none) align with the fitted treatment
#' coding; a non-binary exposure is rejected (this chunk supports the binary
#' average treatment effect only). The weighted binomial GLM reports fractional
#' "successes" — an expected consequence of non-integer weights — so its
#' `non-integer #successes` warning is muffled here.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"ccw_gformula"`,
#'   carrying the case-control `data`, the `outcome` / `exposure` / `confounders`
#'   roles, and a `design` whose `prevalence` (q0) is set. The
#'   `matchatr_missing_prevalence` guard in [matcha()] ensures q0 is present
#'   before this runs.
#' @returns A `causatr_fit` object (the weighted g-computation fit) stored in the
#'   `matchatr_fit`'s `model` slot; [contrast()] turns it into the marginal
#'   effect.
#' @family estimators
#' @seealso [matcha()], [contrast()], `cc_weights()`, [causatr::causat()]
#' @noRd
fit_ccw <- function(fit) {
  prevalence <- fit$design$prevalence

  # The g-formula standardizes a confounder-adjusted outcome model; with no
  # confounders there is nothing to standardize over (and causatr's g-computation
  # requires an outcome-model adjustment set). Reject early with a matchatr error
  # rather than leaking causatr's.
  if (is.null(fit$confounders)) {
    rlang::abort(
      c(
        "The CCW g-formula requires `confounders` to standardize over.",
        i = paste0(
          "Supply an adjustment set, e.g. `confounders = ~ age + smoke`, on ",
          "`matcha()`."
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }

  # Coerce both roles to 0/1: the outcome so the weighted binomial GLM reads a
  # proper response, the exposure so the marginal contrast's static(1)/static(0)
  # interventions match the treatment coding. A non-binary exposure has no
  # binary average-treatment-effect contrast and is rejected here.
  y01 <- resolve_binary_outcome(fit$data, fit$outcome)
  x01 <- resolve_binary_exposure(
    fit$data,
    fit$exposure,
    estimator_label = "CCW g-formula",
    alternative = "a conditional estimator (e.g. estimator = \"logistic\")"
  )

  dt <- data.table::copy(data.table::as.data.table(fit$data))
  dt[[fit$outcome]] <- y01
  dt[[fit$exposure]] <- x01

  weights <- cc_weights(prevalence, y01)

  # The outcome model fitter (NULL -> stats::glm); a user may pass mgcv::gam for
  # smooth confounder adjustment, mirroring the glm_logistic engine.
  model_fn <- fit$details$model_fn
  if (is.null(model_fn)) {
    model_fn <- stats::glm
  }

  # Fractional case-control weights make the binomial "successes" non-integer;
  # the resulting GLM warning is expected and benign, so muffle just that one.
  withCallingHandlers(
    causatr::causat(
      dt,
      outcome = fit$outcome,
      treatment = fit$exposure,
      confounders = fit$confounders,
      estimator = "gcomp",
      family = "binomial",
      weights = as.numeric(weights),
      model_fn = model_fn
    ),
    warning = function(w) {
      if (grepl("non-integer #successes", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
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
#' @param ci_method Character variance source; `"model"` / `"sandwich"` forward
#'   to causatr's sandwich; `"bootstrap"` is rejected.
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
          "The CCW g-formula reports a marginal effect, not `type = \"",
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
        '`ci_method = "bootstrap"` is not available for the CCW g-formula.',
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
    ci_method = ci_method,
    reference = "control",
    n = res$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = res$vcov,
    call = call
  )
}
