#' Classical-estimator dispatch table
#'
#' Maps each sampling design to the classical (non-causal) estimators it
#' admits and the engine key each resolves to. The case-control-weighted
#' (`"ccw_*"`) estimators are handled separately in `resolve_engine()` because
#' they apply to *any* design (they reweight a sample back to the source
#' population), so they are not enumerated per design here.
#'
#' @returns A named list keyed by design type; each element is a named
#'   character vector mapping an estimator string to its engine key.
#' @family dispatch
#' @noRd
dispatch_table <- function() {
  list(
    unmatched_cc = c(
      logistic = "glm_logistic",
      mh = "mantel_haenszel",
      # Multiple case / control groups: a multinomial logistic with a shared
      # reference, each non-reference equation a subtype-vs-reference log OR.
      polytomous = "multinom"
    ),
    # Matched CC and NCC share the conditional partial-likelihood engine: a
    # matched set and a sampled risk set are the same stratum construction. The
    # McNemar closed form is the 1:1 binary-exposure reduction of that
    # likelihood, offered as a faster, formula-level alternative for paired data.
    matched_cc = c(clogit = "clogit", mcnemar = "mcnemar"),
    # The risk-set conditional likelihood (clogit) plus the three Samuelsen-IPW
    # estimators that break the matching: weighted Cox (hazard ratio), weighted
    # AFT (time ratio), and weighted Aalen additive hazards (excess hazard).
    nested_cc = c(
      clogit = "clogit",
      ipw_cox = "ipw_cox",
      ipw_aft = "ipw_aft",
      ipw_aalen = "ipw_aalen"
    ),
    case_cohort = c(cch = "cch"),
    two_phase = c(survey = "survey_twophase"),
    counter_matched = c(weighted_cox = "weighted_cox")
  )
}

#' The case-control-weighted estimator names
#'
#' These reuse the causatr g-computation / IPW / AIPW engines (and the new
#' targeting step for TMLE) by passing q0-based case-control weights as
#' observation weights; they are valid on any design but require a prevalence.
#'
#' @returns A character vector of the `"ccw_*"` estimator strings.
#' @family dispatch
#' @noRd
ccw_estimators <- function() {
  c("ccw_gformula", "ccw_ipw", "ccw_aipw", "ccw_tmle")
}

#' Default estimator for a design
#'
#' The canonical, design-faithful analysis used when the caller does not name
#' an `estimator` in [matcha()].
#'
#' @param design_type Character scalar design type.
#' @returns A character scalar estimator name.
#' @family dispatch
#' @noRd
default_estimator <- function(design_type) {
  switch(
    design_type,
    unmatched_cc = "logistic",
    matched_cc = "clogit",
    nested_cc = "clogit",
    case_cohort = "cch",
    two_phase = "survey",
    counter_matched = "weighted_cox",
    # Defensive: a design type with no default should never reach here because
    # the constructors are the only source of `type`.
    rlang::abort(
      paste0("No default estimator for design type `", design_type, "`."),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  )
}

#' Default contrast scale for an engine
#'
#' The estimand an engine identifies, used by [contrast()] when the caller does
#' not name a `type`. The classical odds-ratio engines identify only the
#' conditional OR, so they default to `"or"`; engines that target a marginal
#' effect default to the risk difference (the etverse convention). The
#' conditional partial-likelihood engine serves two designs through the same
#' `"clogit"` engine key, and the design fixes the scale: a matched case-control
#' design reports the odds ratio, a risk-set-sampled nested case-control design
#' the hazard ratio (OR = HR exactly under risk-set sampling).
#'
#' @param engine Character scalar engine key.
#' @param design_type Character scalar design type, or `NULL`. Distinguishes the
#'   matched (`"or"`) and nested (`"hr"`) case-control designs, which share the
#'   `"clogit"` engine; ignored by every other engine.
#' @returns A character scalar contrast type (`"or"`, `"hr"`, or `"difference"`).
#' @family dispatch
#' @noRd
default_contrast_type <- function(engine, design_type = NULL) {
  # The nested case-control design reports the hazard ratio; the matched design
  # the odds ratio. Both resolve to the `"clogit"` engine, so the design type
  # disambiguates the default scale.
  if (identical(engine, "clogit") && identical(design_type, "nested_cc")) {
    return("hr")
  }
  switch(
    engine,
    glm_logistic = "or",
    mantel_haenszel = "or",
    clogit = "or",
    mcnemar = "or",
    multinom = "or",
    # Counter-matched partial likelihood identifies the hazard ratio (Langholz &
    # Borgan 1995), not the OR. The weighted coxph returns exp(beta); the design
    # fixes the scale. IPW NCC weighted Cox similarly reports the Cox HR (the
    # IPW weighted partial likelihood is consistent for the same estimand).
    weighted_cox = "hr",
    ipw_cox = "hr",
    cch = "hr",
    # The IPW NCC alternative-model engines each identify one non-Cox scale: the
    # AFT a time ratio (acceleration factor), the additive Aalen model an excess
    # hazard (a rate difference).
    ipw_aft = "af",
    ipw_aalen = "excess",
    # Case-control weighting targets a marginal effect; the etverse convention is
    # to default to the risk difference.
    ccw_gformula = "difference",
    ccw_ipw = "difference",
    ccw_aipw = "difference",
    "difference"
  )
}

#' Outcome encoding an engine expects
#'
#' Most engines analyse a binary case indicator; the polytomous (multinomial)
#' engine analyses a multi-group outcome instead. [matcha()] reads this to pick
#' the outcome resolver — `resolve_binary_outcome()` versus
#' `resolve_polytomous_outcome()`.
#'
#' @param engine Character scalar engine key.
#' @returns `"polytomous"` for the multinomial engine, otherwise `"binary"`.
#' @family dispatch
#' @noRd
engine_outcome_kind <- function(engine) {
  if (identical(engine, "multinom")) "polytomous" else "binary"
}

#' Resolve a (design, estimator) pair to an engine
#'
#' Implements the two orthogonal axes: the `design_type` selects the sampling
#' structure and the `estimator` selects the analysis. Case-control-weighted
#' estimators short-circuit (valid on any design); otherwise the classical
#' dispatch table is consulted and an unrecognised estimator is rejected with a
#' classed error that lists the admissible choices.
#'
#' @param design_type Character scalar design type.
#' @param estimator Character scalar estimator name.
#' @param call Caller environment surfaced in the error.
#' @returns A named list with `engine` (the engine key), `kind`
#'   (`"classical"` or `"ccw"`), and `conditional` (`TRUE` for the
#'   conditional-likelihood engine, which triggers stratum-informativeness
#'   checks).
#' @family dispatch
#' @noRd
resolve_engine <- function(design_type, estimator, call = rlang::caller_env()) {
  ccw <- ccw_estimators()
  if (estimator %in% ccw) {
    return(list(
      engine = estimator,
      kind = "ccw",
      conditional = FALSE,
      outcome_kind = "binary"
    ))
  }

  allowed <- dispatch_table()[[design_type]]
  if (is.null(allowed) || !estimator %in% names(allowed)) {
    valid <- c(names(allowed), ccw)
    rlang::abort(
      c(
        paste0(
          "Estimator `",
          estimator,
          "` is not available for design `",
          design_type,
          "`."
        ),
        i = paste0(
          "Supported estimators for this design: ",
          paste0("\"", valid, "\"", collapse = ", "),
          "."
        )
      ),
      class = c("matchatr_bad_estimator", "matchatr_error"),
      call = call
    )
  }

  engine <- unname(allowed[[estimator]])
  list(
    engine = engine,
    kind = "classical",
    # Both the conditional logistic and its 1:1 McNemar reduction condition on
    # the matched sets, so both want the uninformative-stratum check that the
    # `conditional` flag triggers in `matcha()`.
    conditional = estimator %in% c("clogit", "mcnemar"),
    # The polytomous engine analyses a multi-group outcome; everything else a
    # binary case indicator. `matcha()` resolves the outcome accordingly.
    outcome_kind = engine_outcome_kind(engine)
  )
}

#' Run the resolved estimation engine on a fit
#'
#' Bridges [matcha()] to the estimator implementations: it switches on the
#' fit's resolved `engine` key and returns the fitted model object. Engines
#' with no wired estimator return `NULL`, leaving the fit's `model` slot empty
#' so [contrast()] and the summary methods guard on it with the classed
#' `matchatr_not_estimated` condition.
#'
#' @param fit A `matchatr_fit` with `model = NULL`, carrying the validated
#'   analysis specification (data, outcome / exposure roles, confounders).
#' @returns The fitted model object (e.g. a `glm`), or `NULL` for an engine
#'   without a wired estimator.
#' @family dispatch
#' @noRd
run_engine <- function(fit) {
  switch(
    fit$engine,
    glm_logistic = fit_logistic_cc(fit),
    mantel_haenszel = fit_mh(fit),
    clogit = fit_clogit(fit),
    mcnemar = fit_mcnemar(fit),
    multinom = fit_polytomous(fit),
    weighted_cox = fit_weighted_cox(fit),
    ipw_cox = fit_ipw_cox(fit),
    ipw_aft = fit_ipw_aft(fit),
    ipw_aalen = fit_ipw_aalen(fit),
    cch = fit_cch(fit),
    # Case-control-weighted causal estimators: reweight to the source population
    # and delegate the marginal estimate to causatr (g-computation / IPW / AIPW).
    # All three route through the same fit_ccw(), which reads fit$estimator. The
    # CCW-TMLE engine is not yet wired and falls through to NULL.
    ccw_gformula = fit_ccw(fit),
    ccw_ipw = fit_ccw(fit),
    ccw_aipw = fit_ccw(fit),
    NULL
  )
}

#' Data columns a design references
#'
#' Collects every `data` column named by a design (strata, time, subcohort,
#' phase-1 strata, phase-2 selection) so [matcha()] can check they all exist
#' before resolving the engine. `NULL` slots drop out of the concatenation.
#'
#' @param design A `matchatr_design` object.
#' @returns A character vector of unique column names (possibly length 0).
#' @family dispatch
#' @noRd
design_columns <- function(design) {
  cols <- c(
    design$strata,
    design$time,
    design$subcohort,
    design$weights,
    design$id,
    design$stratum,
    design$phase1,
    design$phase2
  )
  # `c()` of all-NULL slots is NULL; coerce so the return is always a character
  # vector (possibly length 0) for predictable downstream length() / setdiff().
  as.character(unique(cols))
}
