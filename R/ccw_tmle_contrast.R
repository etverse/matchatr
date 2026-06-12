#' Assemble the CCW-TMLE marginal contrast
#'
#' Turns a targeted CCW-TMLE fit into a `matchatr_result` reporting the marginal
#' risk difference (`type = "difference"`), risk ratio (`type = "ratio"`), or
#' marginal odds ratio (`type = "or"`) with the efficient-influence-function
#' variance weighted by the case-control weights.
#'
#' @details
#' With the case-control weights treated as fixed (known q0), the variance of a
#' marginal mean is the weighted EIF variance Var(ψ̂) = Σ (wᵢ Dᵢ)² / n², where n
#' is the sample size and Σ wᵢ = n. The risk difference uses D = D₁ − D₀; the
#' risk ratio and odds ratio use the delta-method log-scale influence functions
#' D₁/ψ₁ − D₀/ψ₀ and D₁/(ψ₁(1−ψ₁)) − D₀/(ψ₀(1−ψ₀)), with the interval formed on
#' the log scale and exponentiated. `ci_method = "bootstrap"` instead reports the
#' within-stratum percentile interval (`ccw_bootstrap_ci()`); the point estimate
#' stays the targeted plug-in.
#'
#' @param fit A `matchatr_fit` whose `model` is a `matchatr_ccw_tmle` object.
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"or"`.
#' @param ci_method Character variance source; `"model"` / `"sandwich"` both use
#'   the EIF variance (recorded as `"sandwich"`); `"bootstrap"` reports the
#'   within-stratum percentile interval.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @param n_boot Integer number of bootstrap replicates when
#'   `ci_method = "bootstrap"`.
#' @returns A `matchatr_result` carrying the targeted intervention means and the
#'   marginal contrast with EIF (or within-stratum bootstrap) variance.
#' @family estimators
#' @seealso [contrast()], `fit_ccw_tmle()`, `ccw_bootstrap_ci()`
#' @noRd
contrast_ccw_tmle <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env(),
  n_boot = 1000L
) {
  if (!type %in% c("difference", "ratio", "or")) {
    rlang::abort(
      c(
        paste0(
          "A case-control-weighted estimator reports a marginal effect, not `type = \"",
          type,
          "\"`."
        ),
        i = paste0(
          'Use `type = "difference"` (risk difference), `"ratio"` (risk ',
          'ratio), or `"or"` (marginal odds ratio).'
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  m <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  ey1 <- m$EY1
  ey0 <- m$EY0
  # EIF SE of a weighted marginal mean (or any linear combination of the D's):
  # sqrt(Σ (w D)²) / n, with Σ w = n (the case-control weights treated as fixed).
  se_eif <- function(d) sqrt(sum((m$weights * d)^2)) / m$n

  se_y1 <- se_eif(m$D1)
  se_y0 <- se_eif(m$D0)
  estimates <- data.table::data.table(
    intervention = c("treated", "control"),
    estimate = c(ey1, ey0),
    se = c(se_y1, se_y0),
    ci_lower = c(ey1 - z * se_y1, ey0 - z * se_y0),
    ci_upper = c(ey1 + z * se_y1, ey0 + z * se_y0)
  )

  if (identical(type, "difference")) {
    est <- ey1 - ey0
    se <- se_eif(m$D1 - m$D0)
    lower <- est - z * se
    upper <- est + z * se
    estimand <- "marginal risk difference"
  } else if (identical(type, "ratio")) {
    log_est <- log(ey1) - log(ey0)
    se_log <- se_eif(m$D1 / ey1 - m$D0 / ey0)
    est <- exp(log_est)
    lower <- exp(log_est - z * se_log)
    upper <- exp(log_est + z * se_log)
    # OR/RR-scale `se` is the delta-method value (RR * SE(log RR)); the interval
    # is the log-scale Wald exponentiated, so `se` does not reconstruct it.
    se <- est * se_log
    estimand <- "marginal risk ratio"
  } else {
    log_est <- stats::qlogis(ey1) - stats::qlogis(ey0)
    se_log <- se_eif(m$D1 / (ey1 * (1 - ey1)) - m$D0 / (ey0 * (1 - ey0)))
    est <- exp(log_est)
    lower <- exp(log_est - z * se_log)
    upper <- exp(log_est + z * se_log)
    se <- est * se_log
    estimand <- "marginal odds ratio"
  }

  # The EIF plug-in variance is the influence-function / sandwich variance, the
  # same family the other CCW engines report; `"bootstrap"` overrides the interval
  # with the within-stratum percentile bootstrap (the targeted point is kept).
  recorded_ci <- "sandwich"
  if (identical(ci_method, "bootstrap")) {
    boot <- ccw_bootstrap_ci(fit, type, conf_level, n_boot)
    se <- boot$se
    lower <- boot$lower
    upper <- boot$upper
    recorded_ci <- "bootstrap"
  }

  new_matchatr_result(
    estimates = estimates,
    contrasts = data.table::data.table(
      comparison = "treated vs control",
      estimate = est,
      se = se,
      ci_lower = lower,
      ci_upper = upper
    ),
    type = type,
    estimand = estimand,
    ci_method = recorded_ci,
    reference = "control",
    n = m$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = NULL,
    call = call
  )
}
