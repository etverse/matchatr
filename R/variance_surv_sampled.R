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
