#' Classical-estimator dispatch table
#'
#' Maps each sampling design to the classical (non-causal) estimators it
#' admits and the engine key each resolves to. The case-control-weighted
#' (`"ccw_*"`) estimators are handled separately in [resolve_engine()] because
#' they apply to *any* design (they reweight a sample back to the source
#' population), so they are not enumerated per design here.
#'
#' @returns A named list keyed by design type; each element is a named
#'   character vector mapping an estimator string to its engine key.
#' @family dispatch
#' @noRd
dispatch_table <- function() {
  list(
    unmatched_cc = c(logistic = "glm_logistic", mh = "mantel_haenszel"),
    # Matched CC and NCC share the conditional partial-likelihood engine: a
    # matched set and a sampled risk set are the same stratum construction.
    matched_cc = c(clogit = "clogit"),
    nested_cc = c(clogit = "clogit"),
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
    return(list(engine = estimator, kind = "ccw", conditional = FALSE))
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

  list(
    engine = unname(allowed[[estimator]]),
    kind = "classical",
    conditional = identical(estimator, "clogit")
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
    design$phase1,
    design$phase2
  )
  # `c()` of all-NULL slots is NULL; coerce so the return is always a character
  # vector (possibly length 0) for predictable downstream length() / setdiff().
  as.character(unique(cols))
}
