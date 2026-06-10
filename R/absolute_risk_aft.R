#' Absolute risk from an IPW nested case-control AFT fit
#'
#' Engine for `absolute_risk.matchatr_fit` when `fit$engine == "ipw_aft"`.
#' Computes the cumulative incidence F_x(t) = 1 − S(t | x) directly from the
#' fitted parametric accelerated-failure-time survival curve, with a delta-method
#' confidence interval. Unlike the Cox-type engines (`cch`, `ipw_cox`) there is no
#' Breslow step function: the parametric AFT gives a closed-form S(t | x), so F is
#' smooth in `t`. All four AFT baselines are supported (`weibull`, `exponential`,
#' `lognormal`, `loglogistic`).
#'
#' @details
#' `survival::survreg` parameterises any of its AFT baselines as a
#' log-location-scale model log T = η + σ ε, with ε a standardised error whose
#' survivor function fixes the shape: extreme-value for `"weibull"` /
#' `"exponential"`, Gaussian for `"lognormal"`, logistic for `"loglogistic"`.
#' Writing the **standardised residual**
#'
#'   z(t | x) = (log t − η) / σ,   η = x̃ᵀβ̂,   σ̂ = `fit$model$scale`,
#'
#' the cumulative incidence is F(t | x) = G(z), where G is the error CDF
#' (G_ev(z) = 1 − exp(−exp(z)) for weibull/exponential, Φ(z) for lognormal,
#' plogis(z) for loglogistic). z is linear in the regression coefficients, so the
#' delta method is exact in the gradient: with θ = (β, log σ),
#'
#'   ∂z/∂β = −x̃ / σ,   ∂z/∂(log σ) = −z,
#'
#' and Var(z) = gᵀ V g, where V is the robust Lin-Wei sandwich
#' `survival::survreg(robust = TRUE)` stores in `vcov()`. The Wald interval
#' z ± k·SE(z) is mapped to the risk scale through the (monotone) G by
#' `aft_risk_ci()`; for weibull/exponential this is exactly the
#' complementary-log-log inversion the Cox-type engines use (`cloglog_risk_ci()`).
#' The `"exponential"` baseline fixes σ = 1, so `vcov()` carries no log-scale
#' parameter and the scale term drops out of the gradient (correctly, since σ is
#' not estimated).
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
#' @seealso [absolute_risk()], `fit_ipw_aft()`, `aft_risk_ci()`
#' @noRd
absolute_risk_aft <- function(fit, newdata, times, conf_level = 0.95) {
  model <- fit$model
  beta <- stats::coef(model) # regression coefficients, including (Intercept)
  sigma <- model$scale # AFT scale σ (> 0; fixed at 1 for the exponential)
  aft_dist <- model$dist # the survreg baseline distribution name
  # vcov() is the robust sandwich under robust = TRUE: over θ = (β, log σ); the
  # trailing parameter is the log scale (absent when σ is fixed, e.g. exponential).
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

      # Standardised residual z = (log t - η) / σ.
      z_std <- (log(tt) - eta_r) / sigma

      # Delta-method gradient of z over θ = (β, log σ), aligned to vcov's columns.
      # `scale_names` is empty for a fixed-scale baseline (exponential), so the
      # scale term is correctly skipped.
      grad <- numeric(length(param_names))
      names(grad) <- param_names
      grad[names(beta)] <- -x_r / sigma
      grad[scale_names] <- -z_std
      var_z <- as.numeric(t(grad) %*% vcov_theta %*% grad)
      se_z <- sqrt(max(0, var_z))

      ci <- aft_risk_ci(z_std, se_z, z_crit, aft_dist)
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
    # The interval is the Wald interval on z mapped through the error CDF; for
    # the extreme-value baselines this is the complementary log-log inversion.
    ci_method = if (aft_dist %in% c("weibull", "exponential")) {
      "delta (log-log)"
    } else {
      "delta (standardised residual)"
    },
    engine = fit$engine,
    estimator = fit$estimator,
    method = paste0("IPW AFT (", aft_dist, ")")
  )
}

#' Risk estimate and CI from an AFT standardised residual
#'
#' Maps a Wald interval on the standardised residual z = (log t − η)/σ to the
#' cumulative incidence F = G(z), where G is the survreg baseline's error CDF.
#' Because G is monotone increasing, z − k·SE gives the lower risk bound. For the
#' extreme-value baselines (`weibull`, `exponential`) G(z) = 1 − exp(−exp(z)), so
#' this reduces to the complementary-log-log inversion the Cox-type engines share
#' (`cloglog_risk_ci()`), reused here to keep those distributions byte-identical.
#'
#' @param z Numeric scalar standardised residual z = (log t − η)/σ.
#' @param se_z Numeric scalar standard error of z (delta method).
#' @param z_crit Numeric critical value (e.g. `qnorm(0.975)`).
#' @param dist Character `survreg` baseline distribution name (one of
#'   [aft_supported_dists()]).
#' @returns A list with `$estimate`, `$ci_lower`, `$ci_upper` on the probability
#'   scale.
#' @family contrasts
#' @noRd
aft_risk_ci <- function(z, se_z, z_crit, dist) {
  if (dist %in% c("weibull", "exponential")) {
    # Extreme-value error: F = 1 − exp(−exp(z)); identical to the cloglog path.
    return(cloglog_risk_ci(z, se_z, z_crit))
  }
  # Gaussian (lognormal) or logistic (loglogistic) error CDF.
  cdf <- switch(
    dist,
    lognormal = stats::pnorm,
    loglogistic = stats::plogis
  )
  list(
    estimate = cdf(z),
    ci_lower = max(0, min(1, cdf(z - z_crit * se_z))),
    ci_upper = max(0, min(1, cdf(z + z_crit * se_z)))
  )
}
