#' Check that a value is a single non-empty string
#'
#' @param x Value to check.
#' @param arg Argument name used in the error message.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` if `x` is
#'   not a length-1 character string with at least one character.
#' @family validators
#' @noRd
check_string <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!rlang::is_string(x) || !nzchar(x)) {
    rlang::abort(
      paste0("`", arg, "` must be a single non-empty character string."),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Check that a value is a non-empty character vector
#'
#' Strata may be specified as one or several column names (e.g. frequency
#' matching on `age` and `sex`), so the design layer accepts a character
#' vector rather than a single string.
#'
#' @param x Value to check.
#' @param arg Argument name used in the error message.
#' @param class Error class to attach (defaults to `matchatr_bad_input`).
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts if `x` is not a character vector whose
#'   entries are all non-empty.
#' @family validators
#' @noRd
check_character <- function(
  x,
  arg = rlang::caller_arg(x),
  class = "matchatr_bad_input",
  call = rlang::caller_env()
) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(!nzchar(x))) {
    rlang::abort(
      paste0(
        "`",
        arg,
        "` must be a non-empty character vector of column names."
      ),
      class = c(class, "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Check that a value is a one-sided formula
#'
#' @param x Value to check.
#' @param arg Argument name used in the error message.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` if `x` is
#'   not a formula.
#' @family validators
#' @noRd
check_formula <- function(
  x,
  arg = rlang::caller_arg(x),
  call = rlang::caller_env()
) {
  if (!inherits(x, "formula")) {
    rlang::abort(
      paste0("`", arg, "` must be a one-sided formula (e.g. `~ age + smoke`)."),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Check that one or more columns exist in `data`
#'
#' Column-existence checks happen in [matcha()] rather than in the design
#' constructors because the constructors never see the data — the same design
#' object is meant to be reusable across data frames.
#'
#' @param data A data.frame or data.table.
#' @param cols Character vector of column names that must be present.
#' @param arg Argument name used in the error message.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_design` listing
#'   every missing column.
#' @family validators
#' @noRd
check_cols_exist <- function(
  data,
  cols,
  arg = rlang::caller_arg(cols),
  call = rlang::caller_env()
) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    rlang::abort(
      paste0(
        "Column(s) for `",
        arg,
        "` not found in `data`: ",
        paste0("`", missing, "`", collapse = ", "),
        "."
      ),
      class = c("matchatr_bad_design", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Validate a prevalence (q0) argument
#'
#' The marginal outcome prevalence q0 anchors case-control weighting (Rose &
#' van der Laan): it is the probability of being a case in the source
#' population, so it must be a single probability strictly inside (0, 1).
#' A value of exactly 0 or 1 implies no cases or no controls in the
#' population and makes the weights degenerate, hence the open interval.
#'
#' @param prevalence `NULL` (no q0 supplied) or a single numeric in (0, 1).
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_prevalence` on an
#'   out-of-range or malformed value.
#' @family validators
#' @noRd
check_prevalence <- function(prevalence, call = rlang::caller_env()) {
  if (is.null(prevalence)) {
    return(invisible(NULL))
  }
  ok <- rlang::is_scalar_double(prevalence) ||
    rlang::is_scalar_integer(prevalence)
  if (!ok || is.na(prevalence) || prevalence <= 0 || prevalence >= 1) {
    rlang::abort(
      "`prevalence` (q0) must be a single number strictly between 0 and 1.",
      class = c("matchatr_bad_prevalence", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Validate a matching-ratio argument
#'
#' The matching ratio m is the number of controls sampled per case (m:1).
#' It must be a single whole number >= 1; a ratio below 1 would mean fewer
#' controls than cases, which the matched / risk-set designs do not represent.
#'
#' @param ratio `NULL` (ratio left unspecified) or a single whole number >= 1.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_ratio` on a
#'   non-integer, sub-1, or malformed value.
#' @family validators
#' @noRd
check_ratio <- function(ratio, call = rlang::caller_env()) {
  if (is.null(ratio)) {
    return(invisible(NULL))
  }
  ok <- rlang::is_scalar_double(ratio) || rlang::is_scalar_integer(ratio)
  # `ratio %% 1 != 0` rejects fractional values like 1.5 (R stores literals
  # such as `2` as doubles, so an integer-class test alone is too strict).
  # The `!is.finite(ratio)` guard must come *before* the modulo: `Inf %% 1` is
  # `NaN` and `NaN != 0` is `NA`, which would make the `if` raise an unclassed
  # base error instead of the classed `matchatr_bad_ratio` (review Issue B1).
  if (!ok || is.na(ratio) || !is.finite(ratio) || ratio < 1 || ratio %% 1 != 0) {
    rlang::abort(
      "`ratio` must be a single whole number >= 1 (controls per case).",
      class = c("matchatr_bad_ratio", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

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

  if (is.logical(y)) {
    return(as.integer(y))
  }

  if (is.factor(y)) {
    if (nlevels(y) != 2L) {
      bad_outcome()
    }
    # First level -> 0, second level -> 1, matching the order the user set.
    return(as.integer(y) - 1L)
  }

  if (is.numeric(y)) {
    vals <- unique(stats::na.omit(y))
    if (!all(vals %in% c(0, 1)) || length(vals) < 2L) {
      # `length(vals) < 2` rejects a degenerate column that is all cases or
      # all controls: there is no contrast to estimate from it.
      bad_outcome()
    }
    return(as.integer(y))
  }

  bad_outcome()
}

#' Warn when a conditional-likelihood stratum is uninformative
#'
#' A matched set (matched CC) or sampled risk set (NCC) contributes nothing to
#' the conditional partial likelihood unless it holds at least one case and at
#' least one control: `survival::clogit` silently drops such strata. Rather
#' than let them vanish unnoticed, the design layer warns — naming the count —
#' so the user can decide whether the sampling was as intended.
#'
#' @param strata_list A list of equal-length vectors, one per stratum-defining
#'   column, identifying the stratum each observation belongs to.
#' @param y01 Integer 0/1 case indicator aligned with the strata vectors.
#' @param call Caller environment surfaced in the warning.
#' @returns `NULL` invisibly; emits a `matchatr_uninformative_stratum` warning
#'   when any stratum lacks a case or a control.
#' @family validators
#' @noRd
warn_uninformative_strata <- function(
  strata_list,
  y01,
  call = rlang::caller_env()
) {
  # Drop rows with a missing stratum value (in any column) or unknown case
  # status; they cannot be assigned to an informative set.
  na_any <- Reduce(`|`, lapply(strata_list, is.na))
  keep <- !(na_any | is.na(y01))
  if (!any(keep)) {
    return(invisible(NULL))
  }
  # `interaction(drop = TRUE)` collapses one or several strata columns into a
  # single factor and drops empty combinations; the 2-D table then has one
  # row per observed stratum and columns for control (0) / case (1).
  kept_strata <- lapply(strata_list, function(v) v[keep])
  strat_f <- interaction(kept_strata, drop = TRUE)
  tab <- table(strat_f, factor(y01[keep], levels = c(0L, 1L)))
  n_bad <- sum(tab[, "0"] == 0L | tab[, "1"] == 0L)
  if (n_bad > 0L) {
    rlang::warn(
      c(
        paste0(
          n_bad,
          " stratum/strata have no cases or no controls and carry no ",
          "information for the conditional likelihood."
        ),
        i = "`survival::clogit` drops these sets; check the matched/risk-set sampling."
      ),
      class = c("matchatr_uninformative_stratum", "matchatr_warning"),
      call = call
    )
  }
  invisible(NULL)
}
