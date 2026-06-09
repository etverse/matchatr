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
#' @param strata `NULL` or a non-empty character vector naming the column(s) to
#'   stratify on for the Mantel-Haenszel estimator (`estimator = "mh"`), e.g.
#'   `"agegrp"` or `c("agegrp", "sex")`. Several columns are crossed into a
#'   single stratifying factor. Used only by `"mh"`; the `"logistic"` estimator
#'   adjusts for covariates via `confounders` instead, and `"mh"` with no
#'   `strata` reduces to the crude single-table odds ratio.
#'
#' @returns A `matchatr_design` object of `type` `"unmatched_cc"` carrying the
#'   prevalence, the Mantel-Haenszel strata, and a `weight_spec` describing the
#'   case-control weighting scheme.
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
#' unmatched_cc(strata = "agegrp")
#'
#' @family design constructors
#' @seealso [matched_cc()], [nested_cc()], [case_cohort()], [matcha()]
#' @export
unmatched_cc <- function(prevalence = NULL, strata = NULL) {
  check_prevalence(prevalence)
  if (!is.null(strata)) {
    check_character(strata, class = "matchatr_bad_strata")
  }
  # q0 present -> the sample can be reweighted to the population (CCW);
  # absent -> only the conditional OR / MH estimators apply, no weights.
  weight_spec <- list(
    kind = if (is.null(prevalence)) "none" else "case_control",
    prevalence = prevalence
  )
  new_matchatr_design(
    type = "unmatched_cc",
    strata = strata,
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
#' @param method Character scalar naming the pseudo-likelihood method passed to
#'   [survival::cch()]. One of `"Prentice"` (the default), `"SelfPrentice"`,
#'   `"LinYing"`, `"I.Borgan"`, or `"II.Borgan"`. The first three are for
#'   simple (unstratified) subcohorts; the Borgan variants require stratified
#'   subcohort sampling and need `stratum` to be supplied. `"Prentice"` and
#'   `"SelfPrentice"` share the same asymptotic variance; `"LinYing"` uses an
#'   independent variance estimator; `"I.Borgan"` and `"II.Borgan"` are IPW
#'   estimators for stratified subcohorts with plug-in asymptotic variance.
#' @param id `NULL` or a single character string naming the subject-identifier
#'   column. When `NULL` (the default) the original row positions in the data
#'   are used as synthetic IDs. Supply an ID column when subjects can appear
#'   both as subcohort members and as cases (the common case-cohort situation):
#'   [survival::cch()] uses the ID to correctly pair each subject's two
#'   appearances.
#' @param stratum `NULL` or a non-empty character vector naming the column(s)
#'   defining the subcohort sampling strata. Required when `method` is
#'   `"I.Borgan"` or `"II.Borgan"`: both IPW estimators weight each subject by
#'   the inverse of its stratum-specific subcohort sampling fraction, so the
#'   stratum boundaries must be known. When `method` is `"Prentice"`,
#'   `"SelfPrentice"`, or `"LinYing"` the `stratum` argument is ignored (those
#'   estimators assume a simple random subcohort).
#'
#' @returns A `matchatr_design` object of `type` `"case_cohort"` carrying the
#'   subcohort, time, method, stratum, and (optionally) id columns, with a
#'   `weight_spec` flagged for inclusion-probability weighting.
#'
#' @examples
#' case_cohort(subcohort = "in_subcohort", time = "t")
#' case_cohort(subcohort = "in_subcohort", time = "t", method = "LinYing")
#' case_cohort(subcohort = "in_subcohort", time = "t", id = "subject_id")
#' case_cohort(subcohort = "in_subcohort", time = "t",
#'             method = "I.Borgan", stratum = "region")
#'
#' @family design constructors
#' @seealso [nested_cc()], [matcha()]
#' @export
case_cohort <- function(
  subcohort,
  time,
  method = "Prentice",
  id = NULL,
  stratum = NULL
) {
  check_string(subcohort)
  check_string(time)
  valid_methods <- c(
    "Prentice",
    "SelfPrentice",
    "LinYing",
    "I.Borgan",
    "II.Borgan"
  )
  if (
    !is.character(method) || length(method) != 1L || !method %in% valid_methods
  ) {
    rlang::abort(
      c(
        paste0(
          "`method` must be one of: ",
          paste0('"', valid_methods, '"', collapse = ", "),
          "."
        ),
        i = 'Default is `"Prentice"`.'
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  if (!is.null(id)) {
    check_string(id)
  }
  if (!is.null(stratum)) {
    check_character(stratum, class = "matchatr_bad_input")
  }
  new_matchatr_design(
    type = "case_cohort",
    subcohort = subcohort,
    time = time,
    method = method,
    id = id,
    stratum = stratum,
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
#' Declares a counter-matched nested case-control (NCC) sample: a stratified
#' risk-set design in which controls are drawn from the *opposite* surrogate
#' stratum to the case, so each sampled risk set is balanced across surrogate
#' exposure values. The analysis is a weighted Cox partial likelihood using
#' log-sampling-weights as an offset (Langholz & Borgan 1995); the counter-
#' matching weights enter as `offset(weights)` in `survival::coxph`, and the
#' result is a hazard ratio.
#'
#' @param strata A non-empty character vector naming the column(s) that
#'   identify the sampled risk sets (the matched-set id, e.g. `"set"`). Same
#'   role as the `strata` argument of [nested_cc()].
#' @param time A single character string naming the risk-time column that
#'   defines the failure time for each set (e.g. `"risk_time"`).
#' @param weights `NULL` or a single character string naming the log-weight
#'   column in the data. `sample_ncc_counter_matched()` appends a `"log_w"`
#'   column carrying `log(n_stratum / m_stratum)` for each sampled subject;
#'   supply its name here so the analysis engine can enter it as a Cox offset.
#'   Required for `matcha()` to fit the weighted partial likelihood.
#' @param ratio `NULL` or a single whole number >= 1. Controls sampled per
#'   case (from the opposite surrogate stratum).
#'
#' @returns A `matchatr_design` object of `type` `"counter_matched"` carrying
#'   the strata, time, weights, and ratio, with a `weight_spec` flagged for
#'   counter-matching weights.
#'
#' @examples
#' counter_matched(strata = "set", time = "risk_time", weights = "log_w")
#' counter_matched(strata = "set", time = "risk_time", weights = "log_w",
#'                 ratio = 2L)
#'
#' @family design constructors
#' @seealso [nested_cc()], [sample_ncc_counter_matched()], [matcha()]
#' @export
counter_matched <- function(strata, time, weights = NULL, ratio = NULL) {
  check_character(strata, class = "matchatr_bad_strata")
  check_string(time)
  if (!is.null(weights)) {
    check_string(weights)
  }
  check_ratio(ratio)
  new_matchatr_design(
    type = "counter_matched",
    strata = strata,
    time = time,
    weights = weights,
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
#' @param weights Character scalar log-weight column for counter-matched
#'   designs, or `NULL`.
#' @param method Character scalar estimation method for case-cohort designs
#'   (one of `"Prentice"`, `"SelfPrentice"`, `"LinYing"`, `"I.Borgan"`,
#'   `"II.Borgan"`), or `NULL`.
#' @param id Character scalar subject-identifier column for case-cohort
#'   designs, or `NULL` (row positions used as synthetic IDs).
#' @param stratum Character vector of subcohort sampling stratum column(s) for
#'   Borgan IPW estimators, or `NULL` for simple (unstratified) subcohorts.
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
  weights = NULL,
  method = NULL,
  id = NULL,
  stratum = NULL,
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
      weights = weights,
      method = method,
      id = id,
      stratum = stratum,
      phase1 = phase1,
      phase2 = phase2,
      weight_spec = weight_spec,
      call = call
    ),
    class = "matchatr_design"
  )
}
