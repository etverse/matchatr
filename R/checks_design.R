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

#' Reject a sample with a case that has no eligible control
#'
#' In the [sample_ncc()] generation path a case whose risk set contains no
#' eligible control cannot form a matched set, which signals a sampling failure
#' rather than a property of user-supplied analysis data (where an empty stratum
#' is merely dropped, with a `matchatr_uninformative_stratum` warning). The abort
#' names how many cases are affected and a few of their failure times so the user
#' can trace the cause.
#'
#' @param empty_rows Integer row indices of the cases with an empty risk set.
#' @param tvec Numeric vector of exit / event times.
#' @param call Caller environment surfaced in the error.
#' @returns Never returns; always aborts with class `matchatr_empty_risk_set`.
#' @family validators
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
