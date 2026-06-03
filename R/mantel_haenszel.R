#' Fit the Mantel-Haenszel stratified odds ratio
#'
#' Closed-form Mantel-Haenszel (1959) summary odds ratio for a binary exposure
#' and binary outcome, stratified by the design's `strata` columns (a single
#' crude table when no strata are given). The variance of the log odds ratio is
#' the Robins-Breslow-Greenland (1986) estimator, valid in both the
#' sparse-data (many small strata) and large-strata limits.
#'
#' @details
#' For stratum i with cell counts a_i (exposed cases), b_i (unexposed cases),
#' c_i (exposed controls), d_i (unexposed controls) and total n_i,
#'
#'   OR_MH = sum(a_i d_i / n_i) / sum(b_i c_i / n_i).
#'
#' With R_i = a_i d_i / n_i, S_i = b_i c_i / n_i, P_i = (a_i + d_i) / n_i, and
#' Q_i = (b_i + c_i) / n_i, the Robins-Breslow-Greenland variance of
#' log(OR_MH) is
#'
#'   sum(P_i R_i) / (2 (sum R_i)^2)
#'     + sum(P_i S_i + Q_i R_i) / (2 sum R_i sum S_i)
#'     + sum(Q_i S_i) / (2 (sum S_i)^2).
#'
#' This matches the odds-ratio confidence interval of
#' [stats::mantelhaen.test()] (with `correct = FALSE`).
#'
#' @param fit A `matchatr_fit` whose `engine` resolved to `"mantel_haenszel"`,
#'   carrying the data, the binary `outcome` / `exposure` columns, and the
#'   design's `strata`.
#' @returns A list of class `"matchatr_mh"` with the summary odds ratio (`or`),
#'   its log (`log_or`) and Robins-Breslow-Greenland standard error (`se_log`),
#'   the number of informative strata (`n_strata`), and the analysis size (`n`).
#' @family estimators
#' @seealso [matcha()], [contrast()], [stats::mantelhaen.test()]
#' @noRd
fit_mh <- function(fit) {
  y <- resolve_binary_outcome(fit$data, fit$outcome)
  x <- resolve_binary_exposure(fit$data, fit$exposure)

  # Build the stratifying factor: the crossed design strata, or a single stratum
  # (the crude table) when none are given.
  if (is.null(fit$design$strata)) {
    stratum <- factor(rep(1L, nrow(fit$data)))
  } else {
    strata_cols <- lapply(fit$design$strata, function(col) fit$data[[col]])
    stratum <- interaction(strata_cols, drop = TRUE)
  }

  # Drop rows with a missing outcome, exposure, or stratum.
  keep <- !(is.na(y) | is.na(x) | is.na(stratum))
  y <- y[keep]
  x <- x[keep]
  stratum <- droplevels(stratum[keep])

  cells <- mh_cell_counts(x, y, stratum)
  a <- cells$a
  b <- cells$b
  c <- cells$c
  d <- cells$d
  n <- a + b + c + d

  # Per-stratum products; strata with n_i = 0 are already absent.
  r <- a * d / n
  s <- b * c / n
  sum_r <- sum(r)
  sum_s <- sum(s)
  if (sum_r == 0 || sum_s == 0) {
    rlang::abort(
      c(
        "The Mantel-Haenszel odds ratio is not estimable (a zero exposure-outcome margin).",
        i = "No stratum has both an exposed case and an unexposed control (or vice versa)."
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = fit$call
    )
  }
  or_mh <- sum_r / sum_s

  # Robins-Breslow-Greenland variance of log(OR_MH).
  p <- (a + d) / n
  q <- (b + c) / n
  var_log <- sum(p * r) /
    (2 * sum_r^2) +
    sum(p * s + q * r) / (2 * sum_r * sum_s) +
    sum(q * s) / (2 * sum_s^2)

  structure(
    list(
      or = or_mh,
      log_or = log(or_mh),
      se_log = sqrt(var_log),
      # Informative strata: those contributing to the numerator or denominator
      # (a concordant or all-one-class stratum has r_i = s_i = 0 and is inert).
      n_strata = sum(r > 0 | s > 0),
      n = length(y),
      exposure = fit$exposure
    ),
    class = "matchatr_mh"
  )
}

#' Per-stratum 2x2 cell counts for the Mantel-Haenszel estimator
#'
#' @param x Integer 0/1 exposure (1 = exposed).
#' @param y Integer 0/1 outcome (1 = case).
#' @param stratum A factor giving each observation's stratum.
#' @returns A list of equal-length numeric vectors `a` (exposed cases), `b`
#'   (unexposed cases), `c` (exposed controls), `d` (unexposed controls), one
#'   entry per stratum level.
#' @family estimators
#' @noRd
mh_cell_counts <- function(x, y, stratum) {
  by_stratum <- function(mask) {
    as.numeric(tapply(as.integer(mask), stratum, sum, default = 0L))
  }
  list(
    a = by_stratum(x == 1L & y == 1L),
    b = by_stratum(x == 0L & y == 1L),
    c = by_stratum(x == 1L & y == 0L),
    d = by_stratum(x == 0L & y == 0L)
  )
}

#' Assemble the Mantel-Haenszel odds-ratio contrast
#'
#' Turns a fitted `matchatr_mh` object into a `matchatr_result` carrying
#' the summary odds ratio with a Robins-Breslow-Greenland Wald interval. The
#' risk difference / risk ratio are rejected (unidentified without a prevalence
#' q0); the variance is the closed-form RBG estimator, so the robust-sandwich
#' and bootstrap interval methods do not apply.
#'
#' @param fit A `matchatr_fit` whose `model` is a `matchatr_mh` object.
#' @param type Character contrast scale; only `"or"` is computed.
#' @param ci_method Character interval method; only `"model"` (the RBG Wald
#'   interval) applies.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` with the single Mantel-Haenszel odds ratio.
#' @family estimators
#' @seealso [contrast()], `fit_mh()`
#' @noRd
contrast_mh <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  reject_unidentified_rd_rr(type, call = call)
  # The Mantel-Haenszel odds ratio has one closed-form variance (RBG); the
  # model-vs-robust and bootstrap interval choices do not apply.
  if (!identical(ci_method, "model")) {
    rlang::abort(
      c(
        paste0(
          "`ci_method = \"",
          ci_method,
          "\"` is not available for the Mantel-Haenszel estimator."
        ),
        i = "It reports the Robins-Breslow-Greenland interval; use `ci_method = \"model\"` (the default)."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  mh <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  log_lower <- mh$log_or - z * mh$se_log
  log_upper <- mh$log_or + z * mh$se_log

  estimates <- data.table::data.table(
    term = mh$exposure,
    estimate = mh$log_or, # log OR_MH
    se = mh$se_log,
    ci_lower = log_lower,
    ci_upper = log_upper
  )
  contrasts <- data.table::data.table(
    comparison = mh$exposure,
    estimate = mh$or, # OR_MH
    se = mh$or * mh$se_log, # delta-method SE on the OR scale
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper)
  )
  vcov_mat <- matrix(
    mh$se_log^2,
    nrow = 1,
    dimnames = list(mh$exposure, mh$exposure)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "or",
    estimand = "Mantel-Haenszel OR",
    ci_method = ci_method,
    reference = NULL,
    n = mh$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_mat,
    call = call
  )
}
