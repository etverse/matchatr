#' Fit the 1:1 matched-pair McNemar odds ratio
#'
#' Closed-form conditional odds ratio for an individually matched 1:1
#' case-control sample with a binary exposure (McNemar 1947; Breslow & Day
#' 1980). For one-case-one-control pairs the conditional likelihood reduces to
#' the discordant-pair binomial, so the conditional maximum-likelihood estimate
#' has the closed form OR = n10 / n01 over the discordant pairs — computed here
#' directly, without `survival::clogit`.
#'
#' @details
#' Classify each matched pair by the exposure of its case and its control. Among
#' the pairs discordant on exposure,
#'
#'   n10 = pairs with the case exposed and the control unexposed,
#'   n01 = pairs with the case unexposed and the control exposed.
#'
#' Pairs concordant on exposure (both exposed or both unexposed) cancel from the
#' conditional likelihood and carry no information. The conditional MLE of the
#' odds ratio and the information-matrix variance of its log are
#'
#'   OR = n10 / n01,   Var(log OR) = 1/n10 + 1/n01.
#'
#' This is exactly the 1:1 reduction of the `survival::clogit` conditional
#' likelihood, so it agrees with that engine on the same data while needing no
#' iterative fit. The estimator applies only to genuine 1:1 pairs: a matched set
#' with more than one case or more than one control is M:1 (or richer) matching,
#' which has no two-cell closed form, so it is rejected in favour of the general
#' conditional logistic engine. A pair reduced to a singleton by complete-case
#' dropping, or a set with no case or no control, contributes nothing and is
#' excluded (the `matcha()` entry point already warns about uninformative sets).
#'
#' If every discordant pair points the same way (n10 = 0 or n01 = 0) the
#' boundary MLE is 0 or infinite and the odds ratio is not estimable.
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"mcnemar"`, carrying
#'   the data, the binary `outcome` / `exposure` columns, and the design's
#'   matched-pair `strata`.
#' @returns A list of class `"matchatr_mcnemar"` with the odds ratio (`or`), its
#'   log (`log_or`) and McNemar standard error (`se_log`), the discordant-pair
#'   counts (`n10`, `n01`), the concordant-pair count (`n_concordant`), the
#'   number of complete pairs used (`n_pairs`), the analysis size in individuals
#'   (`n`), and the exposure name.
#' @family estimators
#' @seealso [matcha()], [contrast()], [survival::clogit()]
#' @noRd
fit_mcnemar <- function(fit) {
  y <- resolve_binary_outcome(fit$data, fit$outcome)
  # McNemar is the binary-exposure 1:1 special case; a k>2 / continuous exposure
  # is handled by the general conditional likelihood, so point there (not to the
  # unmatched logistic engine the Mantel-Haenszel resolver suggests).
  x <- resolve_binary_exposure(
    fit$data,
    fit$exposure,
    estimator_label = "McNemar",
    alternative = "estimator = \"clogit\""
  )

  # The matched-set id: one or several crossed columns identify each pair.
  strata_cols <- lapply(fit$design$strata, function(col) fit$data[[col]])

  # Complete-case filter: a pair member missing the outcome, the exposure, or a
  # stratum value cannot be placed in a pair. Drop those rows and report the
  # count, mirroring the conditional-logistic engine.
  na_strata <- Reduce(`|`, lapply(strata_cols, is.na))
  keep <- !(is.na(y) | is.na(x) | na_strata)
  n_dropped <- sum(!keep)
  if (n_dropped > 0L) {
    rlang::warn(
      c(
        paste0(
          n_dropped,
          " row(s) with missing values were dropped from the fit."
        ),
        i = "The odds ratio is estimated on the complete pairs only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  y <- y[keep]
  x <- x[keep]
  stratum <- interaction(
    lapply(strata_cols, function(v) v[keep]),
    drop = TRUE
  )

  # Cases and controls per set. Genuine 1:1 matching has exactly one of each in
  # every set; more than one of either is M:1 (or richer) matching with no
  # two-cell closed form, so reject and point to the conditional logistic engine.
  n_case <- tapply(y == 1L, stratum, sum)
  n_ctrl <- tapply(y == 0L, stratum, sum)
  if (any(n_case > 1L | n_ctrl > 1L, na.rm = TRUE)) {
    n_bad <- sum(n_case > 1L | n_ctrl > 1L, na.rm = TRUE)
    rlang::abort(
      c(
        paste0(
          "The McNemar estimator requires 1:1 matched pairs, but ",
          n_bad,
          " matched set(s) have more than one case or more than one control."
        ),
        i = "Use `estimator = \"clogit\"` for M:1, variable-ratio, or richer matching."
      ),
      class = c("matchatr_not_one_to_one", "matchatr_error"),
      call = fit$call
    )
  }

  # Keep only the complete (one-case, one-control) pairs. Singletons left by
  # complete-case dropping, or sets with no case / no control, are uninformative.
  valid <- names(n_case)[n_case == 1L & n_ctrl == 1L]
  in_valid <- stratum %in% valid
  ys <- y[in_valid]
  xs <- x[in_valid]
  ss <- droplevels(stratum[in_valid])

  # Align the case and control exposure of each pair by ordering both on the set
  # id: each valid set contributes exactly one case row and one control row.
  case_x <- xs[ys == 1L][order(ss[ys == 1L])]
  ctrl_x <- xs[ys == 0L][order(ss[ys == 0L])]

  # Discordant-pair counts and the (inert) concordant count.
  n10 <- sum(case_x == 1L & ctrl_x == 0L)
  n01 <- sum(case_x == 0L & ctrl_x == 1L)
  n_concordant <- sum(case_x == ctrl_x)
  n_pairs <- length(valid)

  # A one-sided (or empty) set of discordant pairs gives a boundary MLE of 0 or
  # infinity: the odds ratio is not identified. Refuse rather than return 0/Inf.
  if (n10 == 0L || n01 == 0L) {
    rlang::abort(
      c(
        "The McNemar odds ratio is not estimable (no discordant pairs in one direction).",
        i = paste0(
          "Discordant pairs: case-exposed/control-unexposed = ",
          n10,
          ", case-unexposed/control-exposed = ",
          n01,
          "."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = fit$call
    )
  }

  or <- n10 / n01
  # Information-matrix variance of the log odds ratio for the discordant-pair
  # binomial (McNemar 1947).
  var_log <- 1 / n10 + 1 / n01

  structure(
    list(
      or = or,
      log_or = log(or),
      se_log = sqrt(var_log),
      n10 = n10,
      n01 = n01,
      n_concordant = n_concordant,
      n_pairs = n_pairs,
      # The analysis size in individuals: the two members of each complete pair.
      n = 2L * n_pairs,
      exposure = fit$exposure
    ),
    class = "matchatr_mcnemar"
  )
}

#' Assemble the McNemar odds-ratio contrast
#'
#' Turns a fitted `matchatr_mcnemar` object into a `matchatr_result` carrying
#' the 1:1 conditional odds ratio with a McNemar Wald interval. The risk
#' difference / risk ratio are rejected (unidentified without a prevalence q0);
#' the variance is the closed-form McNemar estimator, so the robust-sandwich and
#' bootstrap interval methods do not apply.
#'
#' @details
#' The interval is Wald on the log-odds scale, exponentiated, so it is
#' asymmetric on the odds-ratio scale; the OR-scale `se` is the delta-method
#' value OR * SE(log OR), kept for reference, while the reconstructable log-scale
#' estimate and SE live in the result's `estimates`.
#'
#' @param fit A `matchatr_fit` whose `model` is a `matchatr_mcnemar` object.
#' @param type Character contrast scale; only `"or"` is computed.
#' @param ci_method Character interval method; only `"model"` (the McNemar Wald
#'   interval) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` with the single McNemar odds ratio.
#' @family estimators
#' @seealso [contrast()], `fit_mcnemar()`
#' @noRd
contrast_mcnemar <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  # The McNemar odds ratio has one closed-form variance (1/n10 + 1/n01); the
  # model-vs-robust and bootstrap interval choices do not apply.
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the McNemar estimator."
        ),
        i = "It reports the McNemar interval (Var(log OR) = 1/n10 + 1/n01); use `ci_method = \"model\"` (the default)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  mc <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  log_lower <- mc$log_or - z * mc$se_log
  log_upper <- mc$log_or + z * mc$se_log

  estimates <- data.table::data.table(
    term = mc$exposure,
    estimate = mc$log_or, # log OR
    se = mc$se_log,
    ci_lower = log_lower,
    ci_upper = log_upper
  )
  contrasts <- data.table::data.table(
    comparison = mc$exposure,
    estimate = mc$or, # OR
    se = mc$or * mc$se_log, # delta-method SE on the OR scale
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper)
  )
  vcov_mat <- matrix(
    mc$se_log^2,
    nrow = 1,
    dimnames = list(mc$exposure, mc$exposure)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "or",
    estimand = "McNemar OR",
    ci_method = ci_method,
    reference = NULL,
    n = mc$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_mat,
    call = call
  )
}
