#' Compute IPW Breslow absolute risk for an IPW nested case-control fit
#'
#' Engine for `absolute_risk.matchatr_fit` when `fit$engine == "ipw_cox"`.
#' Computes the inverse-probability-weighted Breslow cumulative baseline hazard
#' over the deduplicated NCC analysis sample (cases + Samuelsen-weighted unique
#' controls) and delegates the F_x(t) evaluation and delta-method CI to the
#' shared `assemble_absolute_risk()`. The coefficient variance is the IPW robust
#' sandwich `survival::coxph(robust = TRUE)` reports.
#'
#' @param fit A `matchatr_fit` with engine `"ipw_cox"` and a non-`NULL` model.
#' @param newdata A data frame of covariate patterns.
#' @param times Numeric vector of evaluation times.
#' @param conf_level Numeric confidence level.
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @noRd
absolute_risk_ncc <- function(fit, newdata, times, conf_level = 0.95) {
  beta <- stats::coef(fit$model)
  vcov_beta <- stats::vcov(fit$model)
  breslow <- ipw_breslow_ncc(fit, beta)
  assemble_absolute_risk(
    fit,
    newdata = newdata,
    times = times,
    beta = beta,
    vcov_beta = vcov_beta,
    breslow = breslow,
    conf_level = conf_level,
    method = "IPW"
  )
}

#' IPW Breslow cumulative baseline hazard for IPW nested case-control data
#'
#' Computes the inverse-probability-weighted Breslow cumulative baseline hazard
#' Λ̂₀(t) = Σ_{k: t_k ≤ t} dΛ̂₀(t_k) and its log-scale variance from the same
#' deduplicated, Samuelsen-weighted NCC analysis sample the `ipw_cox` weighted Cox
#' is fitted on (`ncc_ipw_analysis_data()`): each unique cohort subject appears
#' once, with the case weight 1 and each control's weight 1/π_j.
#'
#' At each unique event time t_k the increment is the Horvitz-Thompson Breslow
#' step
#'   dΛ̂₀(t_k) = (Σ_{i: event at t_k} w_i) / (Σ_{j: at risk at t_k} w_j exp(β̂ᵀ x_j))
#' where the weighted at-risk denominator estimates the full-cohort risk set
#' (the controls, upweighted by 1/π_j, stand in for the unsampled cohort). The
#' time axis is the cohort survival time `fit$design$time` — controls are reused
#' across all event times where they are still at risk, not just their matched
#' set's failure time.
#'
#' The log-scale variance uses the Nelson-Aalen delta-method approximation
#'   Var(log Λ̂₀(t)) ≈ Σ_{k: t_k ≤ t} (dΛ̂₀(t_k))^2 / Λ̂₀(t)^2,
#' the within-sample (Poisson) component only; the additional sampling variance
#' from estimating the inclusion weights is not added (conservative).
#'
#' @param fit A `matchatr_fit` with engine `"ipw_cox"`.
#' @param beta Named numeric vector of fitted coefficients from `coef(fit$model)`.
#' @returns A list with:
#'   - `$times`: sorted numeric vector of unique event times (with a `t = 0`
#'     fence post).
#'   - `$cumhaz`: cumulative baseline hazard at each event time.
#'   - `$var_log_cumhaz`: Var(log Λ̂₀(t)) at each event time.
#' @family contrasts
#' @noRd
ipw_breslow_ncc <- function(fit, beta) {
  # The same deduplicated, case-weight-forced sample the weighted Cox was fit on.
  dt <- ncc_ipw_analysis_data(fit)
  time_col <- fit$design$time
  w <- dt[["ipw_weight"]]
  t_col <- dt[[time_col]]
  is_ev <- as.logical(dt[[fit$outcome]])

  # LP = β̂ᵀ x over the analysis sample; coxph coefficient names follow standard
  # R contrasts, so they index the model.matrix columns directly.
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs <- stats::reformulate(c(fit$exposure, conf_terms))
  mm <- stats::model.matrix(rhs, data = dt)
  lp <- as.vector(mm[, names(beta), drop = FALSE] %*% beta)

  # Weighted Breslow step at each unique event time. The numerator sums the case
  # weights (1) of the events; the denominator is the IPW-weighted at-risk set.
  t_events <- sort(unique(t_col[is_ev]))
  n_times <- length(t_events)
  inc <- numeric(n_times)
  for (i in seq_len(n_times)) {
    tk <- t_events[i]
    at_risk <- t_col >= tk
    denom <- sum(w[at_risk] * exp(lp[at_risk]))
    if (!is.finite(denom) || denom <= 0) {
      next
    }
    num <- sum(w[is_ev & t_col == tk])
    inc[i] <- num / denom
  }

  cumhaz <- cumsum(inc)
  inc_sq_cum <- cumsum(inc^2)
  var_log <- ifelse(cumhaz > 0, inc_sq_cum / cumhaz^2, 0)

  # Fence post at t = 0 so times before the first event map to cumhaz = 0 (F = 0).
  list(
    times = c(0, t_events),
    cumhaz = c(0, cumhaz),
    var_log_cumhaz = c(0, var_log)
  )
}
