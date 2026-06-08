#' Coerce a case-status column to a 0/1 integer vector
#'
#' Case-control, nested case-control, and case-cohort designs all have a
#' binary case / event indicator. This resolver accepts the three sane
#' encodings — logical, a two-level factor, or numeric 0/1 — and returns the
#' canonical 0/1 integer vector used by the strata-informativeness checks.
#' Anything else (continuous outcome, 3+ levels, all-equal) is rejected so a
#' mis-specified outcome cannot silently flow into a binary-only estimator.
#'
#' @param data A data.frame or data.table.
#' @param outcome Character scalar naming the case-status column.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of 0/1 the same length as `nrow(data)` (NA
#'   preserved); aborts with class `matchatr_bad_outcome` on a non-binary
#'   column.
#' @family validators
#' @noRd
resolve_binary_outcome <- function(data, outcome, call = rlang::caller_env()) {
  y <- data[[outcome]]

  bad_outcome <- function() {
    rlang::abort(
      c(
        paste0(
          "Outcome `",
          outcome,
          "` must be a binary case indicator (logical, two-level factor, ",
          "or numeric 0/1)."
        ),
        i = "Multiple case / control groups are handled by a polytomous estimator."
      ),
      class = c("matchatr_bad_outcome", "matchatr_error"),
      call = call
    )
  }

  # Coerce each accepted encoding to a 0/1 integer first; the contrast check
  # below is then applied uniformly so logical / factor / numeric cannot
  # diverge in what they let through.
  y01 <- if (is.logical(y)) {
    as.integer(y)
  } else if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      bad_outcome()
    }
    # First level -> 0, second level -> 1, matching the order the user set.
    as.integer(y) - 1L
  } else if (is.numeric(y)) {
    vals <- unique(stats::na.omit(y))
    if (!all(vals %in% c(0, 1))) {
      bad_outcome()
    }
    as.integer(y)
  } else {
    bad_outcome()
  }

  # A case-control analysis needs a contrast: both a case (1) and a control (0)
  # must actually occur. Enforced here for every encoding so an all-cases /
  # all-controls (or all-NA) column cannot reach the estimation layer with a
  # silent n_cases = 0 -- previously only the numeric branch caught this
  # (review Issue R2).
  if (length(unique(stats::na.omit(y01))) < 2L) {
    bad_outcome()
  }
  y01
}

#' Resolve a multi-group outcome for the polytomous estimator
#'
#' The polytomous (multinomial) estimator needs an outcome with three or more
#' groups — multiple disease subtypes, or several control groups. This resolver
#' accepts a factor or character outcome, drops unused factor levels and
#' NA-omits, requires at least three observed groups, and relevels the column so
#' the chosen reference group is the baseline (first) level that
#' [nnet::multinom()] contrasts every other level against.
#'
#' @details
#' A binary outcome (two groups, or a logical / numeric 0/1 column) is rejected
#' with `matchatr_bad_outcome`, pointing back to the binary case-control
#' estimators (`logistic` / `mh` / `clogit`) — the polytomous engine is only for
#' k >= 3. The reference must name one of the observed levels; an absent
#' reference is `matchatr_bad_input`. When `reference` is `NULL` the first level
#' (factor) or the first level in sorted order (character) is used, and the
#' choice is recorded so it can be echoed.
#'
#' @param data A data.frame or data.table.
#' @param outcome Character scalar naming the multi-group outcome column.
#' @param reference `NULL` (use the default baseline) or a character scalar
#'   naming the reference group, which must be one of the observed levels.
#' @param call Caller environment surfaced in the error.
#' @returns A list with `factor` (the reference-first factor, same length as
#'   `nrow(data)`, NA preserved), `levels` (its levels, reference first),
#'   `reference` (the resolved reference group), and `counts` (a named integer
#'   vector of per-group counts). Aborts with `matchatr_bad_outcome` for a
#'   non-categorical or fewer-than-three-group outcome, or `matchatr_bad_input`
#'   for a reference not among the levels.
#' @family validators
#' @noRd
resolve_polytomous_outcome <- function(
  data,
  outcome,
  reference = NULL,
  call = rlang::caller_env()
) {
  y <- data[[outcome]]

  bad_outcome <- function(msg) {
    rlang::abort(
      c(
        msg,
        i = paste0(
          "A two-group outcome is a binary case-control analysis ",
          "(`estimator = \"logistic\"`, `\"mh\"`, or `\"clogit\"`)."
        )
      ),
      class = c("matchatr_bad_outcome", "matchatr_error"),
      call = call
    )
  }

  # A polytomous outcome is categorical: a factor or character column. A logical
  # or numeric 0/1 column is a binary outcome routed to the wrong estimator.
  yf <- if (is.factor(y)) {
    droplevels(y)
  } else if (is.character(y)) {
    factor(y)
  } else {
    bad_outcome(paste0(
      "Outcome `",
      outcome,
      "` must be a factor or character column with three or more groups."
    ))
  }

  levs <- levels(yf)
  if (length(levs) < 3L) {
    bad_outcome(paste0(
      "Outcome `",
      outcome,
      "` has ",
      length(levs),
      " group(s); the polytomous estimator needs at least three."
    ))
  }

  # Resolve the reference (baseline) group; multinom contrasts every other level
  # against the first level, so it is releveled to the front below.
  if (is.null(reference)) {
    reference <- levs[1]
  } else {
    check_string(reference, arg = "reference", call = call)
    if (!reference %in% levs) {
      rlang::abort(
        c(
          paste0(
            "`reference` `",
            reference,
            "` is not one of the outcome groups."
          ),
          i = paste0(
            "Observed groups: ",
            paste0("\"", levs, "\"", collapse = ", "),
            "."
          )
        ),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  yf <- stats::relevel(yf, ref = reference)

  list(
    factor = yf,
    levels = levels(yf),
    reference = reference,
    counts = table(yf)
  )
}

#' Coerce a binary exposure to a 0/1 integer vector
#'
#' The Mantel-Haenszel summary odds ratio and the 1:1 McNemar odds ratio are
#' both defined for a binary exposure (a 2x2 table per stratum / pair). This
#' resolver accepts the same three encodings as `resolve_binary_outcome()` —
#' logical, a two-level factor (second level = exposed), or numeric 0/1 — and
#' returns the 0/1 integer vector. A multi-level or continuous exposure is
#' rejected, pointing to the estimator that handles categorical / continuous
#' exposures (named by `alternative`). Degenerate (single-value) exposures are
#' left to the estimator's zero-margin guard.
#'
#' @param data A data.frame or data.table.
#' @param exposure Character scalar naming the exposure column.
#' @param estimator_label Character scalar naming the estimator in the error
#'   message (e.g. `"Mantel-Haenszel"`, `"McNemar"`).
#' @param alternative Character scalar naming the estimator to use instead for a
#'   non-binary exposure, inserted verbatim into the hint (e.g.
#'   `"estimator = \"logistic\""`).
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of 0/1 (NA preserved); aborts with class
#'   `matchatr_bad_input` on a non-binary exposure.
#' @family validators
#' @noRd
resolve_binary_exposure <- function(
  data,
  exposure,
  estimator_label = "Mantel-Haenszel",
  alternative = "estimator = \"logistic\"",
  call = rlang::caller_env()
) {
  x <- data[[exposure]]
  bad_exposure <- function() {
    rlang::abort(
      c(
        paste0(
          "The ",
          estimator_label,
          " estimator requires a binary exposure; `",
          exposure,
          "` is not binary (logical, two-level factor, or numeric 0/1)."
        ),
        i = paste0(
          "For a categorical (k>2) or continuous exposure use `",
          alternative,
          "`."
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (is.logical(x)) {
    as.integer(x)
  } else if (is.factor(x)) {
    if (nlevels(x) != 2L) {
      bad_exposure()
    }
    as.integer(x) - 1L
  } else if (is.numeric(x)) {
    vals <- unique(stats::na.omit(x))
    if (!all(vals %in% c(0, 1))) {
      bad_exposure()
    }
    as.integer(x)
  } else {
    bad_exposure()
  }
}

#' Coerce a cohort event indicator to a 0/1 integer vector
#'
#' Risk-set sampling ([sample_ncc()]) needs the cohort's event column as a 0/1
#' indicator to locate the cases that anchor each risk set. Unlike
#' `resolve_binary_outcome()` (which requires BOTH a case and a control to occur,
#' because the case-control analysis needs a contrast), this resolver only
#' requires at least one event: controls are sampled by being at risk, not by
#' their eventual event status, so an all-events cohort is still samplable (a
#' future case can serve as an earlier case's control).
#'
#' @param data A data.frame or data.table.
#' @param event Character scalar naming the event column.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of 0/1 the same length as `nrow(data)` (NA
#'   preserved, treated as non-event when locating cases). Aborts with
#'   `matchatr_bad_outcome` on a non-0/1 column or a cohort with no events.
#' @family validators
#' @noRd
resolve_event_indicator <- function(data, event, call = rlang::caller_env()) {
  y <- data[[event]]
  bad_event <- function(msg) {
    rlang::abort(
      c(
        msg,
        i = "Supply a logical, two-level factor, or numeric 0/1 event indicator with at least one event."
      ),
      class = c("matchatr_bad_outcome", "matchatr_error"),
      call = call
    )
  }
  y01 <- if (is.logical(y)) {
    as.integer(y)
  } else if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      bad_event(paste0("Event `", event, "` must be a two-level factor."))
    }
    as.integer(y) - 1L
  } else if (is.numeric(y)) {
    vals <- unique(stats::na.omit(y))
    if (!all(vals %in% c(0, 1))) {
      bad_event(paste0("Event `", event, "` must be numeric 0/1."))
    }
    as.integer(y)
  } else {
    bad_event(paste0(
      "Event `",
      event,
      "` must be a binary indicator (logical, two-level factor, or numeric 0/1)."
    ))
  }
  if (sum(y01, na.rm = TRUE) < 1L) {
    bad_event(paste0(
      "Event `",
      event,
      "` has no events to sample risk sets for."
    ))
  }
  y01
}
