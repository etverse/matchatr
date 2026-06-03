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

#' Reject a column assigned to two incompatible analysis roles
#'
#' The outcome and the exposure each play a single, fixed role and must be
#' distinct from the adjustment covariates and from the columns that define the
#' sampling design. Reusing the outcome as a confounder regresses / standardises
#' the outcome on itself; reusing it as a matched-set id makes every set
#' single-class; entering the exposure as a confounder double-enters the term of
#' interest; matching or stratifying on the exposure collapses the contrast the
#' design exists to estimate.
#'
#' Confounders and design columns are deliberately allowed to overlap: matching
#' (or frequency-matching) on a variable and additionally adjusting for residual
#' confounding by it is a standard, valid combination.
#'
#' @param outcome Character scalar outcome column name.
#' @param exposure Character scalar exposure column name.
#' @param confounder_vars Character vector of confounder column names (possibly
#'   length 0).
#' @param design_cols Character vector of design-referenced column names
#'   (possibly length 0).
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` when the
#'   outcome or exposure also appears as a confounder or a design column.
#' @family validators
#' @noRd
check_role_collisions <- function(
  outcome,
  exposure,
  confounder_vars,
  design_cols,
  call = rlang::caller_env()
) {
  covariates <- unique(c(confounder_vars, design_cols))
  collide <- function(role, col) {
    rlang::abort(
      c(
        paste0(
          "`",
          role,
          "` column `",
          col,
          "` also appears as a confounder or a design column."
        ),
        i = "The outcome and exposure must each be a distinct column from the covariates and design."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (outcome %in% covariates) {
    collide("outcome", outcome)
  }
  if (exposure %in% covariates) {
    collide("exposure", exposure)
  }
  invisible(NULL)
}

#' Reject an ordered-factor exposure
#'
#' matchatr reports the exposure as a single effect (continuous / trend) or one
#' odds ratio per level (unordered factor). An *ordered* factor is fit by
#' `glm` / `mgcv::gam` with polynomial contrasts (`.L`, `.Q`, ...), whose
#' coefficients are not per-level odds ratios, so it is rejected at the
#' user-facing entry point rather than silently producing polynomial-contrast
#' "ORs".
#'
#' @param data A data.frame or data.table.
#' @param exposure Character scalar exposure column name (already validated to
#'   exist).
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` when the
#'   exposure column is an ordered factor.
#' @family validators
#' @noRd
check_exposure_not_ordered <- function(
  data,
  exposure,
  call = rlang::caller_env()
) {
  col <- data[[exposure]]
  if (is.factor(col) && is.ordered(col)) {
    rlang::abort(
      c(
        paste0(
          "Exposure `",
          exposure,
          "` is an ordered factor, which is fit with polynomial contrasts ",
          "(not per-level odds ratios)."
        ),
        i = "Pass a numeric score for a trend OR, or an unordered factor for per-level ORs."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  invisible(NULL)
}

#' Validate an effect-modifier column for stratum-specific odds ratios
#'
#' Effect modification in a matched case-control analysis is expressed as the
#' `exposure x modifier` interaction in the conditional logistic model
#' (`outcome ~ exposure * modifier + strata(set)`); `contrast()` then reports
#' the exposure odds ratio within each modifier level (the stratum-specific OR).
#' This validator enforces the two structural requirements: the modifier must be
#' categorical (so "per level" is well defined), and the exposure must
#' contribute a single coefficient (binary, continuous, or two-level factor) so
#' each modifier level's OR is the single linear combination
#' `beta_x + beta_{x:level}`.
#'
#' @details
#' A continuous modifier has no discrete levels, so it is rejected with a hint
#' to bin it or wrap it in `factor()`. A factor exposure with three or more
#' levels contributes several coefficients, turning effect modification into a
#' grid of level-by-level interactions rather than one OR per modifier level, so
#' it is rejected (`matchatr_unsupported_combination`). The modifier may coincide
#' with a matching / design column or a confounder — assessing whether the
#' exposure OR differs across the matching variable is the canonical use — but it
#' must be a distinct column from the outcome and the exposure.
#'
#' @param data A data.frame or data.table.
#' @param effect_modifier Character scalar naming the modifier column.
#' @param outcome Character scalar outcome column name.
#' @param exposure Character scalar exposure column name.
#' @param call Caller environment surfaced in the error.
#' @returns `NULL` invisibly; aborts with class `matchatr_bad_input` for a
#'   missing / mis-typed / role-colliding modifier, or
#'   `matchatr_unsupported_combination` for a multi-level factor exposure.
#' @family validators
#' @noRd
check_effect_modifier <- function(
  data,
  effect_modifier,
  outcome,
  exposure,
  call = rlang::caller_env()
) {
  check_string(effect_modifier, arg = "effect_modifier", call = call)
  check_cols_exist(data, effect_modifier, arg = "effect_modifier", call = call)
  # The modifier is a covariate, so it may overlap a matching / design column or
  # a confounder, but it cannot double as the outcome or the exposure.
  if (effect_modifier %in% c(outcome, exposure)) {
    rlang::abort(
      c(
        paste0(
          "`effect_modifier` column `",
          effect_modifier,
          "` must differ from the outcome and the exposure."
        ),
        i = "It enters the model as `exposure * effect_modifier`; the exposure is already the term of interest."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  # A categorical modifier defines the levels the per-level OR is reported over;
  # a continuous modifier has none, so the OR would have to be reported at chosen
  # values (out of scope here).
  modcol <- data[[effect_modifier]]
  if (!(is.factor(modcol) || is.character(modcol) || is.logical(modcol))) {
    rlang::abort(
      c(
        paste0(
          "`effect_modifier` `",
          effect_modifier,
          "` must be categorical (logical, character, or factor)."
        ),
        i = "Bin a continuous modifier or wrap it in `factor()` to report one odds ratio per level."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  # Each modifier level's OR is the single linear combination
  # beta_x + beta_{x:level}, which presumes the exposure contributes ONE
  # coefficient. A 3+-level factor exposure contributes several, making effect
  # modification a level-by-level grid rather than one OR per modifier level.
  xcol <- data[[exposure]]
  x_levels <- if (is.factor(xcol)) {
    nlevels(droplevels(xcol))
  } else if (is.character(xcol)) {
    length(unique(stats::na.omit(xcol)))
  } else {
    NA_integer_
  }
  if (!is.na(x_levels) && x_levels > 2L) {
    rlang::abort(
      c(
        paste0(
          "Effect modification is supported for a single-coefficient exposure ",
          "(binary, continuous, or two-level factor), but `",
          exposure,
          "` has ",
          x_levels,
          " levels."
        ),
        i = "Report per-level odds ratios for a multi-level exposure without an effect modifier, or dichotomise it."
      ),
      class = c("matchatr_unsupported_combination", "matchatr_error"),
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
