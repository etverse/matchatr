#' Assemble a design-weighted marginal survival contrast by g-computation
#'
#' Turns a fitted design-weighted Cox model (`fit_surv_gcomp()`) into a
#' `matchatr_result` reporting a marginal causal-survival contrast at the
#' requested follow-up times: the risk difference (`type = "difference"`), risk
#' ratio (`type = "ratio"`), or restricted-mean-survival-time difference
#' (`type = "rmst"`). The subject-specific absolute risk F(t | x, W) from the
#' design's weighted Cox + inverse-probability Breslow baseline is standardized
#' (Horvitz-Thompson averaged with the inclusion weights) over the cohort
#' covariate distribution to the treat-all and treat-none risks, then contrasted.
#'
#' @details
#' The marginal risk under the static intervention do(X = a) is
#' F^a(t) = Σ_j w_j F(t | X = a, W_j) / Σ_j w_j, the inclusion-weighted average
#' of the fitted absolute risk over the standardization sample (the subcohort for
#' a case-cohort design, the deduplicated risk-set sample for a nested
#' case-control design). The risk difference and ratio read F^a at each requested
#' time; the RMST difference integrates the marginal survival S^a = 1 − F^a over
#' the interval from 0 to the horizon. Variance is the design-preserving bootstrap (the inclusion
#' weights are not refit, so there is no information-matrix interval); the
#' conditional odds ratio / hazard ratio and the IPW NCC scales are rejected.
#'
#' @param fit A `matchatr_fit` whose `model` is a weighted Cox from
#'   `fit_surv_gcomp()`.
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"rmst"`.
#' @param ci_method Character variance source. `"model"` and `"bootstrap"` both
#'   give the design-preserving bootstrap (the engine's only variance);
#'   `"sandwich"` is not available.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param times Numeric vector of follow-up times (risk difference / ratio) or
#'   horizons (RMST). Required.
#' @param call Caller environment surfaced in any error.
#' @param n_boot Integer number of bootstrap replicates (default 500).
#' @returns A `matchatr_result` carrying the marginal intervention risks and the
#'   per-time marginal contrast with a bootstrap percentile interval.
#' @family causal survival
#' @seealso [contrast()], `fit_surv_gcomp()`, [absolute_risk()]
#' @noRd
contrast_surv_gcomp <- function(
  fit,
  type,
  ci_method,
  conf_level,
  times = NULL,
  call = rlang::caller_env(),
  n_boot = 500L
) {
  if (!type %in% c("difference", "ratio", "rmst")) {
    rlang::abort(
      c(
        paste0(
          "A design-weighted causal-survival fit reports a marginal survival contrast, not `type = \"",
          type,
          "\"`."
        ),
        i = paste0(
          'Use `type = "difference"` (risk difference), `"ratio"` (risk ratio), ',
          'or `"rmst"` (restricted mean survival time difference).'
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (is.null(times) || length(times) == 0L) {
    rlang::abort(
      c(
        "A design-weighted causal-survival contrast requires `times`.",
        i = "Pass the follow-up time(s) (or RMST horizon[s]) at which to report the contrast, e.g. `times = c(2, 5)`."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (!is.numeric(times) || anyNA(times) || any(times <= 0)) {
    rlang::abort(
      "`times` must be a numeric vector of positive, non-missing values.",
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  if (identical(ci_method, "sandwich")) {
    rlang::abort(
      c(
        '`ci_method = "sandwich"` is not available for the design-weighted survival estimator.',
        i = 'Use `ci_method = "bootstrap"` (the design-preserving percentile interval).'
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  # Recode the exposure to 0/1 so the standardization sample and the
  # marginalization match the (0/1-coded) fitted model.
  fit <- surv_gcomp_recode_exposure(fit)
  std <- surv_gcomp_std_sample(fit)
  fit_view <- fit
  fit_view$engine <- std$engine

  point <- surv_gcomp_point(fit_view, std, fit$exposure, type, times)
  boot <- surv_gcomp_boot_ci(fit, type, times, conf_level, n_boot)

  estimand <- switch(
    type,
    difference = "marginal risk difference",
    ratio = "marginal risk ratio",
    rmst = "marginal RMST difference"
  )
  contrasts <- data.table::data.table(
    comparison = "treat-all vs treat-none",
    time = times,
    estimate = point$contrast,
    se = boot$se,
    ci_lower = boot$lower,
    ci_upper = boot$upper
  )
  estimates <- data.table::data.table(
    intervention = rep(c("treat-all", "treat-none"), each = length(times)),
    time = rep(times, 2L),
    risk = c(point$f1, point$f0)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = type,
    estimand = estimand,
    ci_method = "bootstrap",
    reference = NULL,
    n = nrow(std$newdata),
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = NULL,
    call = call
  )
}

#' Marginal-risk point estimates for a design-weighted survival contrast
#'
#' Computes the treat-all / treat-none marginal risks and their contrast on the
#' requested scale by standardizing the fitted absolute risk over the
#' standardization sample.
#'
#' @param fit_view A `matchatr_fit` whose `engine` is the underlying weighted-Cox
#'   engine (`"cch"` / `"ipw_cox"`), used for the `absolute_risk()` dispatch.
#' @param std The list from `surv_gcomp_std_sample()` (`newdata`, `weights`).
#' @param exposure Character scalar exposure column name.
#' @param type Character contrast scale (`"difference"`, `"ratio"`, `"rmst"`).
#' @param times Numeric vector of evaluation times / horizons.
#' @returns A list with `f1`, `f0` (the treat-all / treat-none marginal risk at
#'   each time, or the marginal survival is reported via `f = 1 - risk` for RMST)
#'   and `contrast` (the contrast on the requested scale).
#' @family causal survival
#' @seealso `contrast_surv_gcomp()`
#' @noRd
surv_gcomp_point <- function(fit_view, std, exposure, type, times) {
  if (identical(type, "rmst")) {
    rmst <- vapply(
      times,
      function(h) surv_gcomp_rmst_diff(fit_view, std, exposure, h),
      numeric(3L)
    )
    # rmst rows: RMST^1, RMST^0, difference.
    return(list(f1 = rmst[1L, ], f0 = rmst[2L, ], contrast = rmst[3L, ]))
  }

  mr <- surv_gcomp_marginal_risk(fit_view, std, exposure, times)
  contrast <- if (identical(type, "ratio")) mr$f1 / mr$f0 else mr$f1 - mr$f0
  list(f1 = mr$f1, f0 = mr$f0, contrast = contrast)
}

#' Inclusion-weighted marginal risk under treat-all and treat-none
#'
#' Evaluates the fitted absolute risk for every standardization-sample row with
#' the exposure set to 1 and to 0, then Horvitz-Thompson averages each over the
#' inclusion weights to the marginal risks F^1(t) and F^0(t).
#'
#' @param fit_view A `matchatr_fit` keyed to the underlying weighted-Cox engine.
#' @param std The list from `surv_gcomp_std_sample()`.
#' @param exposure Character scalar exposure column name.
#' @param times Numeric vector of evaluation times.
#' @returns A list with `f1` and `f0`, each a numeric vector aligned to `times`.
#' @family causal survival
#' @seealso `surv_gcomp_point()`
#' @noRd
surv_gcomp_marginal_risk <- function(fit_view, std, exposure, times) {
  beta <- stats::coef(fit_view$model)
  breslow <- surv_gcomp_breslow(fit_view, beta)
  # Cumulative baseline hazard at each requested time (step function: zero before
  # the first failure, last value after the last; right-continuous).
  cumhaz_t <- stats::approx(
    breslow$times,
    breslow$cumhaz,
    xout = times,
    method = "constant",
    f = 0,
    rule = 2
  )$y
  cumhaz_t[is.na(cumhaz_t)] <- 0

  nd1 <- std$newdata
  nd1[[exposure]] <- 1
  nd0 <- std$newdata
  nd0[[exposure]] <- 0
  r1 <- exp(ar_lp_from_newdata(fit_view, nd1, beta)$lp)
  r0 <- exp(ar_lp_from_newdata(fit_view, nd0, beta)$lp)

  w <- std$weights
  w_tot <- sum(w)
  # F^a(t) = 1 - Σ_j w_j exp(-Λ̂₀(t) exp(β̂ᵀ x_j^a)) / Σ_j w_j: one minus the
  # inclusion-weighted mean subject survival under do(X = a).
  marg <- function(r) {
    vapply(
      cumhaz_t,
      function(lambda) 1 - sum(w * exp(-lambda * r)) / w_tot,
      numeric(1)
    )
  }
  list(f1 = marg(r1), f0 = marg(r0))
}

#' Inverse-probability Breslow cumulative baseline hazard for a sampled design
#'
#' Dispatches to the design's IPW Breslow estimator (case-cohort or nested
#' case-control) so the design-weighted g-computation reads the same baseline as
#' [absolute_risk()].
#'
#' @param fit_view A `matchatr_fit` keyed to the underlying weighted-Cox engine.
#' @param beta Named numeric coefficient vector from the weighted Cox.
#' @returns A Breslow step list with `times`, `cumhaz`, and `var_log_cumhaz`.
#' @family causal survival
#' @noRd
surv_gcomp_breslow <- function(fit_view, beta) {
  switch(
    fit_view$engine,
    cch = ipw_breslow_cch(fit_view, beta),
    ipw_cox = ipw_breslow_ncc(fit_view, beta),
    rlang::abort(
      paste0("No IPW Breslow baseline for engine `", fit_view$engine, "`."),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  )
}

#' Marginal RMST difference up to a horizon
#'
#' Integrates the marginal survival S^a = 1 − F^a over the interval from 0 to the horizon under
#' treat-all and treat-none and returns both RMSTs and their difference. The
#' integration grid is the fit's failure times up to the horizon (the marginal
#' survival is a step function on the Breslow jumps), so the trapezoidal sum is
#' the exact area under the step curve.
#'
#' @param fit_view A `matchatr_fit` keyed to the underlying weighted-Cox engine.
#' @param std The list from `surv_gcomp_std_sample()`.
#' @param exposure Character scalar exposure column name.
#' @param horizon Numeric scalar RMST horizon.
#' @returns A length-3 numeric vector: `RMST^1`, `RMST^0`, and their difference.
#' @family causal survival
#' @seealso `surv_gcomp_point()`
#' @noRd
surv_gcomp_rmst_diff <- function(fit_view, std, exposure, horizon) {
  ev <- surv_event_times(fit_view)
  grid <- sort(unique(c(ev[ev > 0 & ev <= horizon], horizon)))
  mr <- surv_gcomp_marginal_risk(fit_view, std, exposure, grid)
  # Prepend t = 0 where S = 1 under both interventions.
  tt <- c(0, grid)
  s1 <- c(1, 1 - mr$f1)
  s0 <- c(1, 1 - mr$f0)
  rmst1 <- sum(diff(tt) * (utils::head(s1, -1L) + utils::tail(s1, -1L)) / 2)
  rmst0 <- sum(diff(tt) * (utils::head(s0, -1L) + utils::tail(s0, -1L)) / 2)
  c(rmst1, rmst0, rmst1 - rmst0)
}

#' Sorted unique failure times of a fit's analysis data
#'
#' @param fit A `matchatr_fit` carrying the outcome / design `time` columns.
#' @returns A sorted numeric vector of distinct event times (`outcome == 1`).
#' @family causal survival
#' @noRd
surv_event_times <- function(fit) {
  tcol <- fit$design$time
  d <- as.integer(fit$data[[fit$outcome]])
  sort(unique(fit$data[[tcol]][!is.na(d) & d == 1L]))
}
