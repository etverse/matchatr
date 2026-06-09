#' Compute IPW Breslow absolute risk for a cch fit
#'
#' Engine for `absolute_risk.matchatr_fit` when `fit$engine == "cch"`. Computes
#' the IPW Breslow cumulative baseline hazard from the case-cohort
#' pseudo-likelihood fit and delegates the F_x(t) evaluation and delta-method CI
#' to the shared `assemble_absolute_risk()`.
#'
#' @param fit A `matchatr_fit` with engine `"cch"` and a non-`NULL` model.
#' @param newdata A data frame of covariate patterns.
#' @param times Numeric vector of evaluation times.
#' @param conf_level Numeric confidence level.
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @noRd
absolute_risk_cch <- function(fit, newdata, times, conf_level = 0.95) {
  beta <- stats::coef(fit$model)
  vcov_beta <- stats::vcov(fit$model)
  breslow <- ipw_breslow_cch(fit, beta)
  assemble_absolute_risk(
    fit,
    newdata = newdata,
    times = times,
    beta = beta,
    vcov_beta = vcov_beta,
    breslow = breslow,
    conf_level = conf_level,
    method = fit$design$method %||% "Prentice"
  )
}

#' IPW Breslow cumulative baseline hazard for case-cohort data
#'
#' Computes the inverse-probability-weighted Breslow cumulative baseline hazard
#' Λ̂₀(t) = Σ_{k: t_k ≤ t} dΛ̂₀(t_k) and its log-scale variance for subsequent
#' delta-method CI construction.
#'
#' For simple (unstratified) subcohorts (Prentice / SelfPrentice / LinYing),
#' the IPW denominator at each event time t_k is
#'   Ã(t_k) = (N / n_sub) × Σ_{j ∈ subcohort, t_j ≥ t_k} exp(β̂ᵀ x_j)
#' where N is the full-cohort size and n_sub the subcohort size, giving
#'   dΛ̂₀(t_k) = n_events(t_k) / Ã(t_k).
#' For Borgan I/II (stratified subcohorts), each subject is weighted by the
#' inverse of its stratum-specific subcohort sampling fraction N_s / n_sub_s
#' (Borgan et al. 2000).
#'
#' The log-scale variance uses the Nelson-Aalen delta-method approximation:
#'   Var(log Λ̂₀(t)) ≈ Σ_{k: t_k ≤ t} (dΛ̂₀(t_k))^2 / Λ̂₀(t)^2.
#' This is the within-sample (Poisson) component only; the additional sampling-
#' variance term from the subcohort draw is not included (conservative CI for
#' smaller subcohort fractions).
#'
#' @param fit A `matchatr_fit` with engine `"cch"`.
#' @param beta Named numeric vector of fitted coefficients from `coef(fit$model)`.
#' @returns A list with:
#'   - `$times`: sorted numeric vector of unique event times (with a `t = 0`
#'     fence post).
#'   - `$cumhaz`: cumulative baseline hazard at each event time.
#'   - `$var_log_cumhaz`: Var(log Λ̂₀(t)) at each event time.
#' @family contrasts
#' @noRd
ipw_breslow_cch <- function(fit, beta) {
  dt <- fit$data
  sc_col <- fit$design$subcohort
  time_col <- fit$design$time
  event_col <- fit$outcome
  stratum_cols <- fit$design$stratum
  method <- fit$design$method %||% "Prentice"

  # Subcohort indicator (logical); event indicator (0/1 resolved by matcha())
  is_sc <- {
    v <- dt[[sc_col]]
    if (is.logical(v)) v else v != 0L
  }
  is_ev <- as.logical(dt[[event_col]])
  t_col <- dt[[time_col]]
  N <- nrow(dt)
  n_sub <- sum(is_sc)

  # Build LP = β̂ᵀ x for all subjects in the full cohort. The coefficient names
  # from cch (standard R contrasts) match model.matrix column names directly.
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs <- stats::reformulate(c(fit$exposure, conf_terms))
  mm_full <- stats::model.matrix(rhs, data = dt)
  lp_full <- as.vector(mm_full[, names(beta), drop = FALSE] %*% beta)

  # IPW subject weights: N_s / n_sub_s per stratum (Borgan) or N / n_sub (simple).
  is_borgan <- method %in% c("I.Borgan", "II.Borgan")
  if (is_borgan && !is.null(stratum_cols)) {
    strat_fac <- if (length(stratum_cols) == 1L) {
      factor(dt[[stratum_cols]])
    } else {
      interaction(dt[, stratum_cols, drop = FALSE], sep = ":")
    }
    strat_N <- table(strat_fac)
    strat_n_sub <- table(strat_fac[is_sc])
    ipw_w <- as.numeric(
      strat_N[as.character(strat_fac)] /
        strat_n_sub[as.character(strat_fac)]
    )
  } else {
    ipw_w <- rep(N / n_sub, N)
  }

  # Breslow step function at each unique event time
  t_events <- sort(unique(t_col[is_ev]))
  n_times <- length(t_events)
  inc <- numeric(n_times)

  for (i in seq_len(n_times)) {
    tk <- t_events[i]
    sc_at_risk <- is_sc & t_col >= tk
    if (!any(sc_at_risk)) {
      next
    }
    n_ev_tk <- sum(is_ev & t_col == tk)
    denom <- sum(ipw_w[sc_at_risk] * exp(lp_full[sc_at_risk]))
    inc[i] <- n_ev_tk / denom
  }

  cumhaz <- cumsum(inc)
  # Nelson-Aalen variance for log(Lambda_0): delta-method approximation.
  # Var(Lambda_0(t)) ~ sum_{k <= t} dLambda_0(t_k)^2.
  # Var(log Lambda_0(t)) ~ Var(Lambda_0(t)) / Lambda_0(t)^2 by delta method.
  inc_sq_cum <- cumsum(inc^2)
  var_log <- ifelse(cumhaz > 0, inc_sq_cum / cumhaz^2, 0)

  breslow_step_with_fence(t_events, cumhaz, var_log)
}
