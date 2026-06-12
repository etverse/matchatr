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

#' Reject a data frame with duplicated column names
#'
#' `data.table::as.data.table()` preserves duplicate column names, and `[[` then
#' silently resolves a name to its *first* match — so a duplicated outcome /
#' exposure / strata name would feed a silently-chosen column into the analysis,
#' and `setdiff()`-based existence checks would even report the wrong column as
#' missing. Reject duplicates up front so every role maps to exactly one column.
#'
#' @param data A data.frame or data.table.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` when any
#'   column name occurs more than once.
#' @family validators
#' @noRd
check_unique_colnames <- function(data, call = rlang::caller_env()) {
  nms <- names(data)
  dup <- unique(nms[duplicated(nms)])
  if (length(dup) > 0L) {
    rlang::abort(
      paste0(
        "`data` has duplicated column name(s): ",
        paste0("`", dup, "`", collapse = ", "),
        ". Rename them so each outcome / exposure / design role maps to one column."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
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

#' Validate a prevalence cohort-size (prevalence_n) argument
#'
#' `prevalence_n` is the cohort size q0 was estimated from. When supplied it makes
#' q0 estimated rather than known, so it requires a `prevalence` to attach to and
#' must be a positive whole number.
#'
#' @param prevalence_n `NULL` or a single positive whole number.
#' @param prevalence The `prevalence` (q0) it qualifies; must be non-`NULL` when
#'   `prevalence_n` is supplied.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_prevalence` on an
#'   invalid value or a `prevalence_n` supplied without a `prevalence`.
#' @family validators
#' @noRd
check_prevalence_n <- function(
  prevalence_n,
  prevalence,
  call = rlang::caller_env()
) {
  if (is.null(prevalence_n)) {
    return(invisible(NULL))
  }
  if (is.null(prevalence)) {
    rlang::abort(
      c(
        "`prevalence_n` was supplied without a `prevalence` (q0).",
        i = "`prevalence_n` is the cohort size q0 was estimated from; supply `prevalence` too."
      ),
      class = c("matchatr_bad_prevalence", "matchatr_error"),
      call = call
    )
  }
  ok <- (rlang::is_scalar_double(prevalence_n) ||
    rlang::is_scalar_integer(prevalence_n)) &&
    !is.na(prevalence_n) &&
    prevalence_n > 0 &&
    prevalence_n == round(prevalence_n)
  if (!ok) {
    rlang::abort(
      "`prevalence_n` (the cohort size q0 was estimated from) must be a single positive whole number.",
      class = c("matchatr_bad_prevalence", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Validate a confidence-level argument
#'
#' A confidence level must be a single probability strictly inside (0, 1): a
#' value of 0 or 1 yields a degenerate (zero-width or infinite) interval.
#'
#' @param conf_level Value to check.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` on an
#'   out-of-range or malformed value.
#' @family validators
#' @noRd
check_conf_level <- function(conf_level, call = rlang::caller_env()) {
  ok <- rlang::is_scalar_double(conf_level) ||
    rlang::is_scalar_integer(conf_level)
  if (!ok || is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort(
      "`conf_level` must be a single number strictly between 0 and 1.",
      class = c("matchatr_bad_input", "matchatr_error"),
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
  if (
    !ok || is.na(ratio) || !is.finite(ratio) || ratio < 1 || ratio %% 1 != 0
  ) {
    rlang::abort(
      "`ratio` must be a single whole number >= 1 (controls per case).",
      class = c("matchatr_bad_ratio", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
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
