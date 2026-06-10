#' Absolute risk from an IPW nested case-control Weibull AFT fit
#'
#' Engine for `absolute_risk.matchatr_fit` when `fit$engine == "ipw_aft"`.
#' Computes the cumulative incidence F_x(t) = 1 − S(t | x) directly from the
#' fitted parametric Weibull accelerated-failure-time survival curve, with a
#' delta-method complementary-log-log confidence interval. Unlike the Cox-type
#' engines (`cch`, `ipw_cox`) there is no Breslow step function: the Weibull AFT
#' gives a closed-form S(t | x), so F is parametric and smooth in `t`.
#'
#' @details
#' For a Weibull AFT, `survival::survreg` parameterises log T = η + σ ε with ε a
#' standard least-extreme-value variate, so the survival function is
#'
#'   S(t | x) = exp(−exp((log t − η) / σ)),
#'
#' where η = x̃ᵀβ̂ is the linear predictor (the model's `(Intercept)` plus
#' covariate effects) and σ̂ = `fit$model$scale`. Writing the cumulative incidence
#' on the complementary log-log scale,
#'
#'   ξ(t | x) = log(−log S(t | x)) = (log t − η) / σ,
#'
#' is linear in the regression coefficients, so the delta method is exact in
#' the gradient. With θ = (β, log σ) the gradient is
#'
#'   ∂ξ/∂β = −x̃ / σ,   ∂ξ/∂(log σ) = −ξ,
#'
#' and Var(ξ) = gᵀ V g, where V is the robust Lin-Wei sandwich
#' `survival::survreg(robust = TRUE)` stores in `vcov()`. The interval ξ ± z·SE(ξ)
#' is inverted to the risk scale by `cloglog_risk_ci()` — the same inversion the
#' Cox-type engines use, sharing the result assembly.
#'
#' Evaluation times at or below the time origin return F = 0 (the survival
#' function is 1 at t = 0; log t is undefined for t ≤ 0).
#'
#' @param fit A `matchatr_fit` with engine `"ipw_aft"` and a non-`NULL` `survreg`
#'   model.
#' @param newdata A data frame of covariate patterns (exposure + confounders).
#' @param times Numeric vector of evaluation times (sorted, de-duplicated here).
#' @param conf_level Numeric confidence level in (0, 1).
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @seealso [absolute_risk()], `fit_ipw_aft()`, `cloglog_risk_ci()`
#' @noRd
absolute_risk_aft <- function(fit, newdata, times, conf_level = 0.95) {
  model <- fit$model
  beta <- stats::coef(model) # regression coefficients, including (Intercept)
  sigma <- model$scale # AFT scale σ (> 0)
  # vcov() is the robust sandwich under robust = TRUE: a (p + 1) × (p + 1) matrix
  # over θ = (β, log σ); the trailing parameter is the log scale.
  vcov_theta <- stats::vcov(model)

  param_names <- colnames(vcov_theta)
  if (is.null(param_names)) {
    # Defensive: align by position when survreg leaves vcov unnamed (β block
    # first, the log-scale parameter(s) last).
    param_names <- c(
      names(beta),
      paste0(".scale", seq_len(ncol(vcov_theta) - length(beta)))
    )
    dimnames(vcov_theta) <- list(param_names, param_names)
  }
  # The scale parameter is whatever sits in vcov beyond the regression block.
  scale_names <- setdiff(param_names, names(beta))

  times <- sort(unique(times))
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Linear predictor η = x̃ᵀβ̂ and the model matrix x̃ (incl. intercept), built
  # from the fitted survreg terms so a data-dependent confounder basis
  # (poly / splines / scale) is reproduced rather than recomputed (see
  # `ar_lp_from_newdata()`).
  lp_info <- ar_lp_from_newdata(fit, newdata, beta)

  rows <- vector("list", nrow(newdata) * length(times))
  k <- 0L
  for (r in seq_len(nrow(newdata))) {
    eta_r <- lp_info$lp[r]
    x_r <- lp_info$mm[r, ] # length p, named by names(beta)

    for (j in seq_along(times)) {
      tt <- times[j]
      if (!is.finite(tt) || tt <= 0) {
        # S(0 | x) = 1 -> F = 0; log t is undefined for t <= 0.
        rows[[k <- k + 1L]] <- list(
          row = r,
          time = tt,
          estimate = 0,
          ci_lower = 0,
          ci_upper = 0
        )
        next
      }

      # ξ = log(-log S) = (log t - η) / σ
      xi <- (log(tt) - eta_r) / sigma

      # Delta-method gradient over θ = (β, log σ), aligned to vcov's columns.
      grad <- numeric(length(param_names))
      names(grad) <- param_names
      grad[names(beta)] <- -x_r / sigma
      grad[scale_names] <- -xi
      var_xi <- as.numeric(t(grad) %*% vcov_theta %*% grad)
      se_xi <- sqrt(max(0, var_xi))

      ci <- cloglog_risk_ci(xi, se_xi, z_crit)
      rows[[k <- k + 1L]] <- list(
        row = r,
        time = tt,
        estimate = ci$estimate,
        ci_lower = ci$ci_lower,
        ci_upper = ci$ci_upper
      )
    }
  }

  new_matchatr_absolute_risk(
    estimates = data.table::rbindlist(rows),
    times = times,
    newdata = newdata,
    conf_level = conf_level,
    ci_method = "delta (log-log)",
    engine = fit$engine,
    estimator = fit$estimator,
    method = "IPW AFT (Weibull)"
  )
}
