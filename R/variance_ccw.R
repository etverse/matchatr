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

  # When q0 is estimated from a cohort of `prevalence_n` members, redraw q0* ~
  # Binomial(N, q0) / N each replicate so the bootstrap captures q̂0's sampling
  # uncertainty on top of the case-control resampling; a known q0 stays fixed.
  q0 <- fit$design$prevalence
  n_cohort <- fit$design$prevalence_n

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
        if (!is.null(n_cohort)) {
          boot_fit$design$prevalence <- stats::rbinom(1L, n_cohort, q0) /
            n_cohort
        }
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
  # Strip `prevalence_n` so the contrast computes the plain point estimate and
  # does NOT re-enter the estimated-q0 / bootstrap variance paths (which both call
  # back here — leaving it set would recurse). The point depends only on q0 (the
  # weights), not on whether q0 is known or estimated.
  fit$design$prevalence_n <- NULL
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

#' Estimated-q0 variance contribution for a CCW marginal contrast
#'
#' When q0 is estimated from a cohort of `prevalence_n` members rather than known,
#' the marginal estimate ψ̂ inherits q̂0's sampling uncertainty through the
#' case-control weights, which the analytic (sandwich / EIF) interval must add.
#'
#' @details
#' By the delta method the extra variance is (∂ψ/∂q0)² Var(q̂0), with
#' Var(q̂0) = q0 (1 − q0) / N for q̂0 = mean(Y) over the N cohort members. The
#' derivative is a central finite difference — ψ̂ refitted at q0 ± h — on the
#' reported scale (the contrast for the risk difference, its log for the risk /
#' odds ratio, matching the scale the analytic SE is formed on).
#'
#' @param fit A `matchatr_fit` whose `design` carries an estimated q0
#'   (`prevalence_n` set).
#' @param type Character contrast scale.
#' @param log_scale Logical; `TRUE` for the ratio / odds-ratio scales, whose SE
#'   lives on the log scale.
#' @returns A numeric scalar — the variance to add to the (log-)scale SE².
#' @family variance
#' @seealso `ccw_bootstrap_ci()`
#' @noRd
ccw_estimated_q0_term <- function(fit, type, log_scale) {
  q0 <- fit$design$prevalence
  n_cohort <- fit$design$prevalence_n
  h <- min(q0, 1 - q0) * 1e-3
  point_at <- function(q) {
    f <- fit
    f$design$prevalence <- q
    p <- ccw_boot_point(f, type)
    if (log_scale) log(p) else p
  }
  d_dq0 <- (point_at(q0 + h) - point_at(q0 - h)) / (2 * h)
  var_q0 <- q0 * (1 - q0) / n_cohort
  d_dq0^2 * var_q0
}

#' Widen an analytic CCW interval for an estimated q0
#'
#' Adds the estimated-q0 variance term (`ccw_estimated_q0_term()`) to the known-q0
#' (sandwich / EIF) standard error and reforms the Wald interval on the reported
#' scale (linear for the risk difference, log for the risk / odds ratio).
#'
#' @param fit A `matchatr_fit` whose `design` carries an estimated q0.
#' @param type Character contrast scale.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param est Numeric point estimate of the contrast.
#' @param lower,upper Numeric known-q0 interval bounds (used to recover the
#'   known-q0 SE on the reported scale).
#' @returns A list with the widened `se`, `lower`, and `upper`.
#' @family variance
#' @seealso `ccw_estimated_q0_term()`
#' @noRd
ccw_apply_estimated_q0 <- function(fit, type, conf_level, est, lower, upper) {
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  log_scale <- type %in% c("ratio", "or")
  term <- ccw_estimated_q0_term(fit, type, log_scale)
  if (log_scale) {
    se_log_known <- (log(upper) - log(lower)) / (2 * z)
    se_log_total <- sqrt(se_log_known^2 + term)
    list(
      se = est * se_log_total,
      lower = exp(log(est) - z * se_log_total),
      upper = exp(log(est) + z * se_log_total)
    )
  } else {
    se_known <- (upper - lower) / (2 * z)
    se_total <- sqrt(se_known^2 + term)
    list(se = se_total, lower = est - z * se_total, upper = est + z * se_total)
  }
}
