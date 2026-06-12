#' Within-stratum bootstrap confidence interval for a CCW marginal contrast
#'
#' Resamples the case-control sample within the case and control strata
#' separately — drawing `n_case` cases with replacement and `n_control` controls
#' with replacement each replicate — refits the case-control-weighted estimator,
#' and returns the percentile interval of the contrast on the requested scale.
#'
#' @details
#' The stratified resample preserves the design: the case and control counts
#' n1 / n0 are fixed across replicates, so the Rose & van der Laan weights
#' q0 / (n1/n) and (1 − q0) / (n0/n) are constant (the known q0 is treated as
#' fixed; its sampling variability is a separate influence-function term). The
#' percentile interval is taken from the bootstrap distribution of the contrast
#' estimate; the reported `se` is the bootstrap standard deviation. The expected
#' `matchatr_dropped_rows` (and, for TMLE, `matchatr_tmle_convergence`) warnings
#' are muffled per replicate.
#'
#' @param fit A `matchatr_fit` whose `engine` is a `ccw_*` estimator, carrying the
#'   case-control `data` and the resolved analysis roles.
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"or"`.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param n_boot Integer number of bootstrap replicates.
#' @returns A list with `se` (bootstrap standard deviation), `lower`, and `upper`
#'   (percentile interval bounds) of the contrast on the requested scale.
#' @family variance
#' @seealso `contrast_ccw()`, `contrast_ccw_tmle()`
#' @noRd
ccw_bootstrap_ci <- function(fit, type, conf_level, n_boot) {
  # Case / control row indices on the original (un-recoded) data; the design is
  # preserved by resampling each stratum to its own size.
  y <- resolve_binary_outcome(fit$data, fit$outcome)
  case_idx <- which(!is.na(y) & y == 1L)
  ctrl_idx <- which(!is.na(y) & y == 0L)
  n_case <- length(case_idx)
  n_ctrl <- length(ctrl_idx)

  ests <- withCallingHandlers(
    vapply(
      seq_len(n_boot),
      function(b) {
        rows <- c(
          case_idx[sample.int(n_case, n_case, replace = TRUE)],
          ctrl_idx[sample.int(n_ctrl, n_ctrl, replace = TRUE)]
        )
        boot_fit <- fit
        boot_fit$data <- fit$data[rows, , drop = FALSE]
        ccw_boot_point(boot_fit, type)
      },
      numeric(1)
    ),
    # Each replicate re-runs ccw_prepare() (which complete-cases) and possibly the
    # TMLE fluctuation; their expected warnings are noise across 1000s of refits.
    matchatr_dropped_rows = function(w) invokeRestart("muffleWarning"),
    matchatr_tmle_convergence = function(w) invokeRestart("muffleWarning")
  )

  ests <- ests[is.finite(ests)]
  alpha <- 1 - conf_level
  qs <- stats::quantile(ests, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  list(se = stats::sd(ests), lower = qs[[1]], upper = qs[[2]])
}

#' Point estimate of a CCW contrast for one bootstrap replicate
#'
#' Runs the resolved engine on a (resampled) fit and returns the contrast point
#' estimate on the requested scale. Uses `ci_method = "model"` so it does not
#' recurse into the bootstrap.
#'
#' @param fit A `matchatr_fit` (typically a bootstrap resample) whose `model` has
#'   not been fitted; this re-runs the engine on `fit$data`.
#' @param type Character contrast scale.
#' @returns A numeric scalar — the contrast estimate, or `NA_real_` if the
#'   replicate's fit failed.
#' @family variance
#' @seealso `ccw_bootstrap_ci()`
#' @noRd
ccw_boot_point <- function(fit, type) {
  fit$model <- tryCatch(run_engine(fit), error = function(e) NULL)
  if (is.null(fit$model)) {
    return(NA_real_)
  }
  res <- tryCatch(
    if (identical(fit$engine, "ccw_tmle")) {
      contrast_ccw_tmle(fit, type, "model", 0.95)
    } else {
      contrast_ccw(fit, type, "model", 0.95)
    },
    error = function(e) NULL
  )
  if (is.null(res)) NA_real_ else res$contrasts$estimate
}
