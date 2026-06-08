#' Draw a nested case-control sample from a cohort by risk-set sampling
#'
#' @description
#' Generates a nested case-control (NCC) dataset from a full cohort by risk-set
#' (incidence-density) sampling: at each event time the failing subject (the
#' case) is matched to `m` controls drawn at random from the subjects still at
#' risk at that instant. Each case and its sampled controls form one matched set
#' (a sampled risk set), which is the stratum the conditional partial likelihood
#' conditions on. The result feeds straight into
#' `matcha(design = nested_cc(strata = "set", time = "risk_time"))`, whose
#' conditional logistic fit reports the hazard ratio (OR = HR exactly under
#' risk-set sampling; Prentice & Breslow 1978).
#'
#' @details
#' The risk set at a failure time `tc` is the set of subjects under observation
#' then: those who have entered follow-up and not yet left, i.e.
#' `entry < tc <= time` (with `entry` taken as the origin when not supplied, so
#' the condition reduces to `time >= tc`). The case is excluded from its own
#' control pool, and `min(m, n_eligible)` controls are sampled without
#' replacement; a late failure time with fewer than `m` eligible controls yields
#' a smaller set rather than an error, mirroring real NCC sampling. A subject
#' sampled as a control may itself fail later in the cohort — it serves as a
#' control before its own event — so the per-set `case` indicator is distinct
#' from the cohort-wide `event` column.
#'
#' Additional matching (`match = ~ s1 + s2`) restricts each case's control pool
#' to subjects sharing the case's values on the named population-stratum
#' variables (e.g. sex, birth cohort). Sampling uses the ambient random-number
#' stream, so wrap the call in [withr::with_seed()] or precede it with
#' [set.seed()] for reproducibility.
#'
#' The sampler is implemented natively (base R + data.table) so it is always
#' available and deterministically seedable; [Epi::ccwc()] is an equivalent
#' external implementation used as a cross-check in the package's tests.
#'
#' A case whose risk set contains no eligible control carries no
#' conditional-likelihood information and signals a sampling failure (a
#' misspecified time origin/scale, an entry/exit mismatch, or over-fine `match`
#' strata), so it aborts with `matchatr_empty_risk_set` rather than silently
#' producing a singleton set.
#'
#' @param cohort A data.frame or data.table with one row per subject. Copied,
#'   never mutated.
#' @param time A single character string naming the exit / event-time column (the
#'   time scale on which risk sets are formed). Must be numeric.
#' @param event A single character string naming the cohort event indicator
#'   (logical, a two-level factor, or numeric 0/1); at least one event must
#'   occur.
#' @param m A single whole number >= 1, the number of controls sampled per case
#'   (default 1).
#' @param match `NULL` (no additional matching) or a one-sided formula naming
#'   population-stratum column(s) controls must share with the case
#'   (e.g. `~ sex + birth_cohort`).
#' @param entry `NULL` (everyone enters at the time origin) or a single character
#'   string naming a delayed-entry / left-truncation column. Must be numeric when
#'   supplied; a subject is at risk at `tc` only if `entry < tc`.
#'
#' @returns A `data.table` with one row per sampled subject: the selected rows of
#'   `cohort` (all original columns) plus `set` (integer matched-set id), `case`
#'   (per-set 0/1 indicator, 1 for the case), and `risk_time` (the set's failure
#'   time). Aborts with `matchatr_empty_risk_set` when any case has no eligible
#'   control.
#'
#' @examples
#' # A small cohort with event times; draw 2 controls per case.
#' cohort <- data.frame(
#'   id = 1:8,
#'   t  = c(2, 5, 1, 8, 3, 9, 4, 7),
#'   d  = c(1, 0, 1, 0, 1, 0, 0, 0),
#'   x  = c(1, 0, 1, 0, 0, 1, 0, 1)
#' )
#' set.seed(1)
#' ncc <- sample_ncc(cohort, time = "t", event = "d", m = 2)
#' ncc
#'
#' # Analyse it: each sampled risk set is a stratum -> hazard ratio.
#' fit <- matcha(ncc, outcome = "case", exposure = "x",
#'               design = nested_cc(strata = "set", time = "risk_time"),
#'               estimator = "clogit")
#' contrast(fit)
#'
#' @seealso [nested_cc()], [matcha()], [Epi::ccwc()]
#' @export
sample_ncc <- function(cohort, time, event, m = 1, match = NULL, entry = NULL) {
  call <- rlang::current_env()
  if (!is.data.frame(cohort)) {
    rlang::abort(
      "`cohort` must be a data.frame or data.table.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  check_string(time, call = call)
  check_string(event, call = call)
  check_sample_m(m, call = call)
  if (!is.null(entry)) {
    check_string(entry, call = call)
  }
  match_vars <- character(0)
  if (!is.null(match)) {
    check_formula(match, call = call)
    match_vars <- all.vars(match)
    if (length(match_vars) == 0L) {
      rlang::abort(
        "`match` must name at least one column (e.g. `~ sex`).",
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  check_unique_colnames(cohort, call = call)
  # Every referenced column must exist before sampling reads it.
  check_cols_exist(cohort, time, arg = "time", call = call)
  check_cols_exist(cohort, event, arg = "event", call = call)
  if (!is.null(entry)) {
    check_cols_exist(cohort, entry, arg = "entry", call = call)
  }
  if (length(match_vars) > 0L) {
    check_cols_exist(cohort, match_vars, arg = "match", call = call)
  }
  # The output columns are appended; a cohort that already carries them would
  # collide (and `[[` would silently resolve to the wrong one downstream).
  added <- c("set", "case", "risk_time")
  clash <- intersect(added, names(cohort))
  if (length(clash) > 0L) {
    rlang::abort(
      c(
        paste0(
          "`cohort` already has column(s) ",
          paste0("`", clash, "`", collapse = ", "),
          ", which `sample_ncc()` appends to the output."
        ),
        i = "Rename them so the sampled set id / case indicator / risk time are unambiguous."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  # A base data.frame copy drives the row subsetting: data.table's `[`
  # row-selection is gated on the calling namespace being "data.table-aware"
  # (cedta), which is brittle under load_all, so the sampler reads columns with
  # `[[`, assembles the result by base row indexing, and converts to data.table
  # only at the return. Adding the output columns to the fresh subset never
  # touches `cohort`, so the user's data is not mutated.
  df <- as.data.frame(cohort)
  tvec <- df[[time]]
  if (!is.numeric(tvec)) {
    rlang::abort(
      paste0("`time` column `", time, "` must be numeric."),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  entryvec <- NULL
  if (!is.null(entry)) {
    entryvec <- df[[entry]]
    if (!is.numeric(entryvec)) {
      rlang::abort(
        paste0("`entry` column `", entry, "` must be numeric."),
        class = c("matchatr_bad_input", "matchatr_error"),
        call = call
      )
    }
  }
  evec <- resolve_event_indicator(df, event, call = call)

  n <- nrow(df)
  rid <- seq_len(n)
  # Cases anchor the risk sets; a case needs a finite failure time to define one.
  case_rows <- rid[which(evec == 1L)]
  if (anyNA(tvec[case_rows]) || any(!is.finite(tvec[case_rows]))) {
    rlang::abort(
      paste0(
        "Some case rows have a missing or non-finite `",
        time,
        "`, so their risk set cannot be formed."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  # Order cases by failure time (ties broken by row) so the set index is stable
  # and independent of the input row order.
  case_rows <- case_rows[order(tvec[case_rows], case_rows)]

  # A single per-row matching key (several match columns crossed) lets the
  # eligibility test be one vector comparison per case. "\r" cannot appear in a
  # printed value, so it is a safe separator that keeps c("1","23") distinct
  # from c("12","3").
  match_key <- NULL
  if (length(match_vars) > 0L) {
    match_key <- do.call(
      paste,
      c(lapply(match_vars, function(v) as.character(df[[v]])), sep = "\r")
    )
  }

  # Pass 1 -- eligibility only, no random draws. Computing every case's eligible
  # pool first means the empty-risk-set abort fires BEFORE any RNG is consumed,
  # so a failed call leaves the random stream untouched.
  elig_list <- vector("list", length(case_rows))
  for (k in seq_along(case_rows)) {
    elig_list[[k]] <- eligible_controls(
      case_row = case_rows[k],
      tvec = tvec,
      entryvec = entryvec,
      match_key = match_key
    )
  }
  empty <- case_rows[lengths(elig_list) == 0L]
  if (length(empty) > 0L) {
    reject_empty_risk_set(empty, tvec, call = call)
  }

  # Pass 2 -- sample min(m, n_eligible) controls per case, accumulating the
  # selected rows and their per-set labels for a single base subset at the end. A
  # late failure time with fewer than m eligible controls keeps all of them (a
  # smaller set), which is correct, not an error.
  all_rows <- integer(0)
  set_id <- integer(0)
  case_ind <- integer(0)
  risk_t <- numeric(0)
  for (k in seq_along(case_rows)) {
    ci <- case_rows[k]
    elig <- elig_list[[k]]
    take <- if (length(elig) > m) elig[sample.int(length(elig), m)] else elig
    members <- c(ci, take)
    all_rows <- c(all_rows, members)
    set_id <- c(set_id, rep.int(k, length(members)))
    case_ind <- c(case_ind, as.integer(members == ci))
    risk_t <- c(risk_t, rep.int(tvec[ci], length(members)))
  }

  out <- df[all_rows, , drop = FALSE]
  out[["set"]] <- set_id
  out[["case"]] <- case_ind
  out[["risk_time"]] <- risk_t
  rownames(out) <- NULL
  data.table::as.data.table(out)
}

#' Compute the eligible control rows for one case's risk set
#'
#' The risk set at the case's failure time `tc` is everyone under observation
#' then -- `entry < tc <= time` -- restricted (when matching) to the case's
#' population stratum, with the case itself removed. `which()` drops the `NA`s a
#' missing time / entry comparison produces, so a subject with unknown timing is
#' simply not eligible.
#'
#' @param case_row Integer row index of the case in the cohort.
#' @param tvec Numeric vector of exit / event times.
#' @param entryvec `NULL` (no delayed entry) or a numeric vector of entry times.
#' @param match_key `NULL` (no matching) or a character vector of per-row
#'   crossed matching keys.
#' @returns An integer vector of eligible control row indices (possibly empty).
#' @family estimators
#' @seealso [sample_ncc()]
#' @noRd
eligible_controls <- function(case_row, tvec, entryvec, match_key) {
  tc <- tvec[case_row]
  at_risk <- tvec >= tc
  if (!is.null(entryvec)) {
    # Delayed entry: a subject is only at risk after entering follow-up.
    at_risk <- at_risk & (entryvec < tc)
  }
  if (!is.null(match_key)) {
    at_risk <- at_risk & (match_key == match_key[case_row])
  }
  at_risk[case_row] <- FALSE
  which(at_risk)
}

#' Coerce a cohort event indicator to a 0/1 integer vector
#'
#' Risk-set sampling needs the cohort's event column as a 0/1 indicator to locate
#' the cases that anchor each risk set. Unlike `resolve_binary_outcome()` (which
#' requires BOTH a case and a control to occur, because the case-control analysis
#' needs a contrast), this resolver only requires at least one event: controls
#' are sampled by being at risk, not by their eventual event status, so an
#' all-events cohort is still samplable (a future case can serve as an earlier
#' case's control).
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

#' Validate the controls-per-case argument for `sample_ncc()`
#'
#' Mirrors `check_ratio()` but is non-optional (the sampler always needs a
#' concrete `m`) and names `m` in the message. The number of controls sampled
#' per case must be a single whole number >= 1.
#'
#' @param m Value to check.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_ratio` on a
#'   `NULL`, non-integer, sub-1, or malformed value.
#' @family validators
#' @noRd
check_sample_m <- function(m, call = rlang::caller_env()) {
  ok <- rlang::is_scalar_double(m) || rlang::is_scalar_integer(m)
  # `!is.finite()` before the modulo: `Inf %% 1` is `NaN` and `NaN != 0` is `NA`,
  # which would raise an unclassed base error from the `if` instead of the
  # classed abort.
  if (is.null(m) || !ok || is.na(m) || !is.finite(m) || m < 1 || m %% 1 != 0) {
    rlang::abort(
      "`m` must be a single whole number >= 1 (controls sampled per case).",
      class = c("matchatr_bad_ratio", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Reject a sample with a case that has no eligible control
#'
#' In the generation path a case whose risk set contains no eligible control
#' cannot form a matched set, which signals a sampling failure rather than a
#' property of user-supplied analysis data (where an empty stratum is merely
#' dropped, with a `matchatr_uninformative_stratum` warning). The abort names how
#' many cases are affected and a few of their failure times so the user can trace
#' the cause.
#'
#' @param empty_rows Integer row indices of the cases with an empty risk set.
#' @param tvec Numeric vector of exit / event times.
#' @param call Caller environment surfaced in the error.
#' @returns Never returns; always aborts with class `matchatr_empty_risk_set`.
#' @family validators
#' @seealso [sample_ncc()]
#' @noRd
reject_empty_risk_set <- function(
  empty_rows,
  tvec,
  call = rlang::caller_env()
) {
  times <- tvec[empty_rows]
  shown <- utils::head(sort(unique(times)), 5L)
  more <- if (length(unique(times)) > length(shown)) ", ..." else ""
  rlang::abort(
    c(
      paste0(
        length(empty_rows),
        " case(s) had no eligible control at their failure time, so a risk set ",
        "could not be formed."
      ),
      i = paste0(
        "Affected failure time(s): ",
        paste(format(shown), collapse = ", "),
        more,
        "."
      ),
      i = paste0(
        "Check the `time` origin/scale, the `entry` column, or whether `match` ",
        "strata are too fine to leave any at-risk control."
      )
    ),
    class = c("matchatr_empty_risk_set", "matchatr_error"),
    call = call
  )
}
