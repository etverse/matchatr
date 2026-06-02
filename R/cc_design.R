#' Unmatched case-control design
#'
#' @description
#' Declares an independent (unmatched) case-control sample: cases and controls
#' are drawn separately from the source population with no individual or
#' frequency matching. This is the design behind the classical logistic /
#' Mantel-Haenszel odds ratio and, when the marginal outcome prevalence q0 is
#' supplied, the case-control-weighted marginal causal contrasts of Rose &
#' van der Laan (2009).
#'
#' @param prevalence `NULL` or a single number in (0, 1). The marginal outcome
#'   prevalence q0 in the source population. Required for case-control-weighted
#'   (`"ccw_*"`) estimators, which reweight the sample back to the population;
#'   optional for the conditional / classical odds-ratio estimators, which do
#'   not need it.
#'
#' @returns A `matchatr_design` object of `type` `"unmatched_cc"` carrying the
#'   prevalence and a `weight_spec` describing the case-control weighting
#'   scheme.
#'
#' @details
#' Weights are never stored as a data column. The `weight_spec` records the
#' intended scheme (case-control weighting when q0 is present, otherwise none);
#' the weighting layer realises it as observation weights into the estimation
#' engine.
#'
#' @examples
#' unmatched_cc()
#' unmatched_cc(prevalence = 0.02)
#'
#' @family design constructors
#' @seealso [matched_cc()], [nested_cc()], [case_cohort()], [matcha()]
#' @export
unmatched_cc <- function(prevalence = NULL) {
  check_prevalence(prevalence)
  # q0 present -> the sample can be reweighted to the population (CCW);
  # absent -> only the conditional OR / MH estimators apply, no weights.
  weight_spec <- list(
    kind = if (is.null(prevalence)) "none" else "case_control",
    prevalence = prevalence
  )
  new_matchatr_design(
    type = "unmatched_cc",
    prevalence = prevalence,
    weight_spec = weight_spec,
    call = match.call()
  )
}

#' Matched case-control design
#'
#' @description
#' Declares an individually or frequency matched case-control sample: each case
#' is matched to one or more controls sharing the values of the matching
#' variable(s). The design-faithful analysis is conditional logistic
#' regression (conditional maximum likelihood), which conditions on the
#' matched-set totals and so removes the matching-variable nuisance
#' parameters; fitting matched-set indicators by ordinary logistic regression
#' biases the odds ratio and is never used.
#'
#' @param strata A non-empty character vector naming the column(s) that
#'   identify the matched sets (e.g. `"set"`, or `c("age_grp", "sex")` for
#'   frequency matching).
#' @param ratio `NULL` or a single whole number >= 1. The number of controls
#'   matched per case (m:1). Optional metadata; the conditional likelihood does
#'   not require a fixed ratio.
#'
#' @returns A `matchatr_design` object of `type` `"matched_cc"` carrying the
#'   strata and ratio.
#'
#' @examples
#' matched_cc(strata = "set")
#' matched_cc(strata = c("age_grp", "sex"), ratio = 2)
#'
#' @family design constructors
#' @seealso [unmatched_cc()], [nested_cc()], [matcha()]
#' @export
matched_cc <- function(strata, ratio = NULL) {
  check_character(strata, class = "matchatr_bad_strata")
  check_ratio(ratio)
  new_matchatr_design(
    type = "matched_cc",
    strata = strata,
    ratio = ratio,
    # Matched CC is fit by conditional likelihood, so no observation weights
    # enter the engine.
    weight_spec = list(kind = "none"),
    call = match.call()
  )
}

#' Nested case-control design
#'
#' @description
#' Declares a nested case-control (NCC) sample drawn from a cohort by
#' risk-set (incidence-density) sampling: at each event time, the case is
#' matched to `ratio` controls sampled from those still at risk. The classical
#' analysis treats each sampled risk set as a stratum in a conditional partial
#' likelihood (`survival::clogit`/`coxph`); because the controls are sampled
#' from the risk set, the conditional odds ratio equals the hazard ratio with
#' no rare-disease assumption.
#'
#' @param strata A non-empty character vector naming the column(s) that
#'   identify the sampled risk sets (the matched-set id).
#' @param time A single character string naming the event/entry time column
#'   that defines the risk sets. Carried for the inclusion-probability weight
#'   and weighted-Cox estimators built on this design.
#' @param ratio `NULL` or a single whole number >= 1. Controls sampled per
#'   case.
#'
#' @returns A `matchatr_design` object of `type` `"nested_cc"` carrying the
#'   strata, time, and ratio, with a `weight_spec` flagged for
#'   inclusion-probability weighting.
#'
#' @examples
#' nested_cc(strata = "set", time = "t")
#' nested_cc(strata = "set", time = "t", ratio = 3)
#'
#' @family design constructors
#' @seealso [case_cohort()], [counter_matched()], [matcha()]
#' @export
nested_cc <- function(strata, time, ratio = NULL) {
  check_character(strata, class = "matchatr_bad_strata")
  check_string(time)
  check_ratio(ratio)
  new_matchatr_design(
    type = "nested_cc",
    strata = strata,
    time = time,
    ratio = ratio,
    # Inclusion-probability (Samuelsen/Borgan) weights are available for the
    # full-cohort weighted-Cox analyses built on an NCC sample.
    weight_spec = list(kind = "inclusion"),
    call = match.call()
  )
}

#' Case-cohort design
#'
#' @description
#' Declares a case-cohort sample: a random subcohort drawn from the full cohort
#' at baseline, augmented by every case occurring in the cohort (including
#' cases outside the subcohort). The design-faithful analyses are the Prentice,
#' Self-Prentice, and Borgan pseudo-likelihood Cox estimators
#' (`survival::cch`), whose controls are reused across failure times so the
#' variance comes from a robust / asymptotic correction, not the information
#' matrix.
#'
#' @param subcohort A single character string naming the 0/1 (or logical)
#'   column flagging subcohort membership.
#' @param time A single character string naming the event/follow-up time column.
#'
#' @returns A `matchatr_design` object of `type` `"case_cohort"` carrying the
#'   subcohort and time columns, with a `weight_spec` flagged for
#'   inclusion-probability weighting.
#'
#' @examples
#' case_cohort(subcohort = "in_subcohort", time = "t")
#'
#' @family design constructors
#' @seealso [nested_cc()], [matcha()]
#' @export
case_cohort <- function(subcohort, time) {
  check_string(subcohort)
  check_string(time)
  new_matchatr_design(
    type = "case_cohort",
    subcohort = subcohort,
    time = time,
    weight_spec = list(kind = "inclusion"),
    call = match.call()
  )
}

#' Two-phase sampling design
#'
#' @description
#' Declares a two-phase (double) sampling design: phase-1 variables are
#' measured on the whole cohort and used to stratify a phase-2 subsample on
#' which the expensive covariate is then measured. The analysis weights
#' phase-2 records by the inverse of their phase-2 selection probability
#' (survey / calibration estimators).
#'
#' @param phase1 A non-empty character vector naming the phase-1 stratification
#'   column(s) measured on the full cohort.
#' @param phase2 A single character string naming the 0/1 (or logical) column
#'   flagging phase-2 selection.
#'
#' @returns A `matchatr_design` object of `type` `"two_phase"` carrying the
#'   phase-1 strata and phase-2 selection columns.
#'
#' @examples
#' two_phase(phase1 = "stratum", phase2 = "in_phase2")
#'
#' @family design constructors
#' @seealso [matcha()]
#' @export
two_phase <- function(phase1, phase2) {
  check_character(phase1)
  check_string(phase2)
  new_matchatr_design(
    type = "two_phase",
    phase1 = phase1,
    phase2 = phase2,
    weight_spec = list(kind = "design"),
    call = match.call()
  )
}

#' Counter-matched design
#'
#' @description
#' Declares a counter-matched nested case-control sample: a stratified
#' risk-set design in which controls are sampled to *differ* from the case on
#' a surrogate of exposure, balancing each sampled risk set across exposure
#' strata. The analysis is a weighted Cox partial likelihood that carries the
#' counter-matching inclusion weights.
#'
#' @param strata A non-empty character vector naming the counter-matching
#'   surrogate stratum column(s).
#' @param time A single character string naming the event/entry time column
#'   that defines the risk sets.
#' @param ratio `NULL` or a single whole number >= 1. Controls sampled per
#'   case within the counter-matched strata.
#'
#' @returns A `matchatr_design` object of `type` `"counter_matched"` carrying
#'   the strata, time, and ratio, with a `weight_spec` flagged for
#'   counter-matching weights.
#'
#' @examples
#' counter_matched(strata = "exposure_surrogate", time = "t")
#'
#' @family design constructors
#' @seealso [nested_cc()], [matcha()]
#' @export
counter_matched <- function(strata, time, ratio = NULL) {
  check_character(strata, class = "matchatr_bad_strata")
  check_string(time)
  check_ratio(ratio)
  new_matchatr_design(
    type = "counter_matched",
    strata = strata,
    time = time,
    ratio = ratio,
    weight_spec = list(kind = "counter_match"),
    call = match.call()
  )
}

#' Construct a `matchatr_design` object
#'
#' Low-level constructor shared by every design constructor. Stores the
#' sampling structure in a fixed slot layout so downstream code (validation,
#' dispatch, print) can read any design uniformly; unused slots are `NULL`.
#'
#' @param type Character scalar design type (one of `"unmatched_cc"`,
#'   `"matched_cc"`, `"nested_cc"`, `"case_cohort"`, `"two_phase"`,
#'   `"counter_matched"`).
#' @param strata Character vector of stratum / matched-set columns, or `NULL`.
#' @param time Character scalar time column, or `NULL`.
#' @param ratio Whole-number controls-per-case, or `NULL`.
#' @param prevalence Marginal prevalence q0, or `NULL`.
#' @param subcohort Character scalar subcohort-membership column, or `NULL`.
#' @param phase1 Character vector of phase-1 strata, or `NULL`.
#' @param phase2 Character scalar phase-2 selection column, or `NULL`.
#' @param weight_spec Named list describing the intended weighting scheme.
#' @param call The originating constructor call, for printing.
#' @returns A list with class `"matchatr_design"`.
#' @family design constructors
#' @noRd
new_matchatr_design <- function(
  type,
  strata = NULL,
  time = NULL,
  ratio = NULL,
  prevalence = NULL,
  subcohort = NULL,
  phase1 = NULL,
  phase2 = NULL,
  weight_spec = list(kind = "none"),
  call = NULL
) {
  structure(
    list(
      type = type,
      strata = strata,
      time = time,
      ratio = ratio,
      prevalence = prevalence,
      subcohort = subcohort,
      phase1 = phase1,
      phase2 = phase2,
      weight_spec = weight_spec,
      call = call
    ),
    class = "matchatr_design"
  )
}
