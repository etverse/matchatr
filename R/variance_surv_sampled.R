#' Design-preserving bootstrap interval for a design-weighted survival contrast
#'
#' Resamples the sampling design, refits the weighted Cox, re-standardizes the
#' marginal contrast, and returns the percentile interval — the variance source
#' for the design-weighted causal-survival g-computation, whose inclusion weights
#' are not refit (so no information-matrix interval applies).
#'
#' @details
#' For a case-cohort design `matcha()` is called with the full cohort, so the
#' bootstrap resamples the cohort rows with replacement: this redraws which
#' subjects are cases and which fall in the subcohort, propagating the full
#' sampling variability (including the subcohort sampling) into the interval. Each
#' replicate refits the weighted Cox via `fit_surv_gcomp()` and re-standardizes
#' the contrast at the requested times; replicates that fail to fit (e.g. a
#' degenerate resample) are dropped from the percentile.
#'
#' @param fit A `matchatr_fit` whose `engine` is `"surv_gcomp"`.
#' @param type Character contrast scale (`"difference"`, `"ratio"`, `"rmst"`).
#' @param times Numeric vector of evaluation times / horizons.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param n_boot Integer number of bootstrap replicates.
#' @returns A list with `se`, `lower`, and `upper`, each a numeric vector aligned
#'   to `times` (the bootstrap SD and percentile interval bounds).
#' @family causal survival
#' @seealso `contrast_surv_gcomp()`, `fit_surv_gcomp()`
#' @noRd
surv_gcomp_boot_ci <- function(fit, type, times, conf_level, n_boot) {
  design_type <- fit$design$type
  dt <- data.table::as.data.table(fit$data)
  n <- nrow(dt)

  if (!identical(design_type, "case_cohort")) {
    # The nested case-control bootstrap redraws risk sets and recomputes the
    # Samuelsen weights, a distinct resampler from the cohort bootstrap.
    rlang::abort(
      paste0(
        "The design-preserving bootstrap is not wired for design `",
        design_type,
        "`."
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error")
    )
  }

  est <- matrix(NA_real_, n_boot, length(times))
  for (b in seq_len(n_boot)) {
    # Cohort resample: redraws case status and subcohort membership together, so
    # cohort.size (nrow) and the subcohort fraction are preserved in expectation.
    rows <- sample.int(n, n, replace = TRUE)
    boot_fit <- fit
    boot_fit$data <- dt[rows, , drop = FALSE]
    boot_model <- tryCatch(fit_surv_gcomp(boot_fit), error = function(e) NULL)
    if (is.null(boot_model)) {
      next
    }
    boot_fit$model <- boot_model
    std <- tryCatch(surv_gcomp_std_sample(boot_fit), error = function(e) NULL)
    if (is.null(std) || nrow(std$newdata) == 0L) {
      next
    }
    boot_view <- boot_fit
    boot_view$engine <- std$engine
    pt <- tryCatch(
      surv_gcomp_point(boot_view, std, fit$exposure, type, times),
      error = function(e) NULL
    )
    if (!is.null(pt)) {
      est[b, ] <- pt$contrast
    }
  }

  boot_percentile_ci(est, conf_level)
}

#' Percentile interval from a bootstrap-replicate matrix, guarding empty columns
#'
#' Reduces a `replicate x time` matrix of bootstrap contrast estimates to the
#' standard error and percentile bounds, dropping failed (`NA`) replicates. A
#' failed refit leaves an `NA` row; `na.rm = TRUE` would otherwise turn a column
#' with *no* successful replicate into a silent `NA` interval, so this aborts when
#' any column is all-`NA` and warns when some replicates were dropped — the caller
#' (and tests) can then match on the classed condition rather than discover a
#' silently missing interval.
#'
#' @param est A numeric `n_boot` by `length(times)` matrix; failed replicates are
#'   `NA` rows.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in the abort.
#' @returns A list with `se`, `lower`, and `upper`, each aligned to the columns of
#'   `est`. Aborts `matchatr_bootstrap_failed` when a column has no successful
#'   replicate; warns `matchatr_bootstrap_failures` when some replicates failed.
#' @family causal survival
#' @seealso `surv_gcomp_boot_ci()`
#' @noRd
boot_percentile_ci <- function(est, conf_level, call = rlang::caller_env()) {
  n_ok <- colSums(!is.na(est))
  if (any(n_ok == 0L)) {
    rlang::abort(
      c(
        "The design-preserving bootstrap produced no successful replicate.",
        i = "Every weighted-Cox refit failed, so the marginal-survival interval is undefined."
      ),
      class = c("matchatr_bootstrap_failed", "matchatr_error"),
      call = call
    )
  }
  n_failed <- nrow(est) - min(n_ok)
  if (n_failed > 0L) {
    rlang::warn(
      c(
        paste0(
          n_failed,
          " bootstrap replicate(s) failed to refit and were dropped."
        ),
        i = "The percentile interval is computed from the successful replicates only."
      ),
      class = c("matchatr_bootstrap_failures", "matchatr_warning")
    )
  }

  alpha <- 1 - conf_level
  se <- apply(est, 2L, stats::sd, na.rm = TRUE)
  lower <- apply(
    est,
    2L,
    function(col) stats::quantile(col, alpha / 2, na.rm = TRUE, names = FALSE)
  )
  upper <- apply(
    est,
    2L,
    function(col) {
      stats::quantile(col, 1 - alpha / 2, na.rm = TRUE, names = FALSE)
    }
  )
  list(se = se, lower = lower, upper = upper)
}
