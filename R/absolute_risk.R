#' Absolute risk from a case-cohort fit
#'
#' @description
#' Estimates the cumulative incidence F_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t))
#' from a fitted case-cohort pseudo-likelihood, where Λ̂₀(t) is the inverse-
#' probability-weighted (IPW) Breslow cumulative baseline hazard (Borgan &
#' Liestøl 1990). Pointwise confidence intervals use the delta method on the
#' complementary log-log scale (the "log-log" CI for the survival function,
#' inverted to the risk scale).
#'
#' @param fit A `matchatr_fit` object returned by [matcha()]. Currently only
#'   the case-cohort (`"cch"`) engine is supported.
#' @param ... Unused; present for S3 consistency.
#' @returns A `matchatr_absolute_risk` object. See `absolute_risk.matchatr_fit`
#'   for details on the return structure.
#' @examples
#' \dontrun{
#' fit <- matcha(cohort, outcome = "d", exposure = "x",
#'               design = case_cohort(subcohort = "sc", time = "t"),
#'               confounders = ~z, estimator = "cch")
#' absolute_risk(fit, newdata = data.frame(x = 1, z = 0), times = c(1, 2, 3))
#' }
#' @family contrasts
#' @seealso [matcha()], [contrast()], [case_cohort()]
#' @export
absolute_risk <- function(fit, ...) UseMethod("absolute_risk")

#' @rdname absolute_risk
#'
#' @param newdata A data frame of covariate values at which to evaluate F_x(t).
#'   Must contain columns matching the exposure and confounders used in `fit`.
#'   Each row is one covariate pattern; the result table carries a `row` column
#'   that indexes back to `newdata`.
#' @param times Non-empty numeric vector of evaluation times. Duplicates are
#'   dropped; times are sorted before evaluation. Times before the first event
#'   return F̂ = 0; times after the last event return the last Breslow value
#'   (step-function extrapolation).
#' @param conf_level Numeric in (0, 1). Confidence level for the pointwise
#'   intervals. Default `0.95`.
#'
#' @returns A `matchatr_absolute_risk` list with elements:
#'   - `$estimates`: a `data.table` with columns `row` (newdata row index),
#'     `time`, `estimate` (F̂_x(t)), `ci_lower`, `ci_upper` (delta-method CI
#'     on the probability scale).
#'   - `$times`: the sorted evaluation times.
#'   - `$newdata`: the supplied covariate patterns.
#'   - `$conf_level`, `$ci_method`, `$engine`, `$estimator`, `$method`.
#'
#' @export
absolute_risk.matchatr_fit <- function(
  fit,
  newdata,
  times,
  conf_level = 0.95,
  ...
) {
  if (!identical(fit$engine, "cch")) {
    rlang::abort(
      c(
        paste0(
          "`absolute_risk()` is not yet implemented for the `",
          fit$engine,
          "` engine."
        ),
        i = "Currently only the case-cohort (`cch`) engine is supported."
      ),
      class = c("matchatr_not_implemented", "matchatr_error")
    )
  }
  if (is.null(fit$model)) {
    rlang::abort(
      "The fit has no estimated model. Run `matcha()` before `absolute_risk()`.",
      class = c("matchatr_not_estimated", "matchatr_error")
    )
  }
  if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
    rlang::abort(
      "`newdata` must be a non-empty data frame.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times))) {
    rlang::abort(
      "`times` must be a non-empty numeric vector of finite values.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  absolute_risk_cch(fit, newdata = newdata, times = times, conf_level = conf_level)
}

#' Compute IPW Breslow absolute risk for a cch fit
#'
#' Internal engine for `absolute_risk.matchatr_fit` when `fit$engine == "cch"`.
#' Computes the IPW Breslow cumulative baseline hazard, evaluates F_x(t) for
#' each row of `newdata` at each requested time, and attaches a delta-method
#' log-log CI.
#'
#' @param fit A `matchatr_fit` with engine `"cch"` and a non-`NULL` model.
#' @param newdata A data frame of covariate patterns.
#' @param times Numeric vector of evaluation times (sorted, de-duplicated
#'   internally).
#' @param conf_level Numeric confidence level.
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @noRd
absolute_risk_cch <- function(fit, newdata, times, conf_level = 0.95) {
  times <- sort(unique(times))
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  beta <- stats::coef(fit$model)
  vcov_beta <- stats::vcov(fit$model)

  # IPW Breslow cumulative baseline hazard and its log-scale variance
  breslow <- ipw_breslow_cch(fit, beta)

  # Interpolate Breslow and its variance to the requested evaluation times.
  # rule = 2 (extrapolate with endpoint value) gives the step-function convention:
  # before first event -> cumhaz = 0, after last event -> last cumhaz value.
  cumhaz_t <- stats::approx(
    breslow$times, breslow$cumhaz,
    xout = times, method = "constant", f = 0, rule = 2
  )$y
  var_log_t <- stats::approx(
    breslow$times, breslow$var_log_cumhaz,
    xout = times, method = "constant", f = 0, rule = 2
  )$y
  # Before the first event, the cumhaz is 0 so the variance is also 0.
  cumhaz_t[is.na(cumhaz_t)] <- 0
  var_log_t[is.na(var_log_t)] <- 0

  # Linear predictor for each newdata row; uses the same coefficient names as
  # the cch fit (standard R contrasts, no "aX" prefix that m$x has internally).
  lp_info <- cch_lp_from_newdata(fit, newdata, beta)

  rows <- vector("list", nrow(newdata) * length(times))
  k <- 0L
  for (r in seq_len(nrow(newdata))) {
    lp_r <- lp_info$lp[r]
    # x_r'  V_beta  x_r — the covariate contribution to Var(log Lambda_x(t)).
    # Computed once per row, not per time (times share the same beta variance).
    x_r <- lp_info$mm[r, , drop = FALSE]
    var_beta_r <- as.numeric(x_r %*% vcov_beta %*% t(x_r))

    for (j in seq_along(times)) {
      # Lambda_x(t) = exp(beta'x) * Lambda_0(t)
      lambda_x <- exp(lp_r) * cumhaz_t[j]

      if (!is.finite(lambda_x) || lambda_x <= 0) {
        # cumhaz_t = 0 before the first event -> F = 0 exactly, no CI needed
        rows[[k <- k + 1L]] <- list(
          row = r, time = times[j],
          estimate = 0, ci_lower = 0, ci_upper = 0
        )
        next
      }

      # Delta method on log(-log(1 - F(t))) = log(Lambda_x(t)) = xi
      # Var(xi) = x' V_beta x  +  Var(log Lambda_0(t))
      # Var(log Lambda_0(t)) ~ Greenwood/Nelson-Aalen variance for the log CHF
      se_xi <- sqrt(var_beta_r + var_log_t[j])

      xi <- log(lambda_x)
      f_est <- 1 - exp(-lambda_x)
      # CI on log-log scale inverted to [0,1]: higher xi -> lower S -> higher F
      ci_lo <- 1 - exp(-exp(xi - z_crit * se_xi))
      ci_hi <- 1 - exp(-exp(xi + z_crit * se_xi))
      # Clamp to [0, 1] to defend against numerical overflow at extreme times
      ci_lo <- max(0, min(1, ci_lo))
      ci_hi <- max(0, min(1, ci_hi))

      rows[[k <- k + 1L]] <- list(
        row = r, time = times[j],
        estimate = f_est, ci_lower = ci_lo, ci_upper = ci_hi
      )
    }
  }

  estimates <- data.table::rbindlist(rows)

  structure(
    list(
      estimates = estimates,
      times = times,
      newdata = newdata,
      conf_level = conf_level,
      ci_method = "delta (log-log)",
      engine = fit$engine,
      estimator = fit$estimator,
      method = fit$design$method %||% "Prentice"
    ),
    class = c("matchatr_absolute_risk", "matchatr")
  )
}

#' IPW Breslow cumulative baseline hazard for case-cohort data
#'
#' Computes the inverse-probability-weighted Breslow cumulative baseline hazard
#' Λ̂₀(t) = Σ_{k: t_k ≤ t} dΛ̂₀(t_k) and its log-scale variance for
#' subsequent delta-method CI construction.
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
#'   - `$times`: sorted numeric vector of unique event times.
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
    ipw_w <- as.numeric(strat_N[as.character(strat_fac)] /
                          strat_n_sub[as.character(strat_fac)])
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
    if (!any(sc_at_risk)) next
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

  # Prepend a fence post at t = 0 so that times before the first event map to
  # cumhaz = 0 (and F = 0) rather than extrapolating to the first event's value.
  list(
    times        = c(0, t_events),
    cumhaz       = c(0, cumhaz),
    var_log_cumhaz = c(0, var_log)
  )
}

#' Linear predictor for new covariate patterns from a cch fit
#'
#' Builds the model matrix for `newdata` using the same formula (exposure +
#' confounders) as the fitted cch model, selects the columns corresponding to
#' the fitted coefficients, and returns both the model matrix and the LP vector.
#'
#' @param fit A `matchatr_fit` with engine `"cch"`.
#' @param newdata A data frame of new covariate patterns.
#' @param beta Named numeric vector of fitted coefficients.
#' @returns A list with `$mm` (model matrix, p columns) and `$lp` (numeric
#'   vector of length `nrow(newdata)`).
#' @family contrasts
#' @noRd
cch_lp_from_newdata <- function(fit, newdata, beta) {
  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs <- stats::reformulate(c(fit$exposure, conf_terms))
  mm_new <- tryCatch(
    stats::model.matrix(rhs, data = newdata),
    error = function(e) {
      rlang::abort(
        c(
          "Could not build a model matrix from `newdata`.",
          i = conditionMessage(e)
        ),
        class = c("matchatr_bad_input", "matchatr_error")
      )
    }
  )
  # Guard: all coefficient names must appear in the model matrix
  missing_cols <- setdiff(names(beta), colnames(mm_new))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      c(
        "Columns in `newdata` do not match the fitted model.",
        i = paste0(
          "Missing or mismatched columns: ",
          paste(missing_cols, collapse = ", ")
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  mm_sel <- mm_new[, names(beta), drop = FALSE]
  list(mm = mm_sel, lp = as.vector(mm_sel %*% beta))
}

#' Print a matchatr_absolute_risk object
#'
#' Displays a compact summary of the absolute-risk result: the engine,
#' pseudo-likelihood method, CI method, the evaluation times, and the number
#' of covariate patterns, followed by the estimates table.
#'
#' @param x A `matchatr_absolute_risk` object returned by [absolute_risk()].
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @family contrasts
#' @seealso [absolute_risk()]
#' @export
print.matchatr_absolute_risk <- function(x, ...) {
  cat("<matchatr_absolute_risk>\n")
  cat(" Engine:     ", x$engine, "\n", sep = "")
  cat(" Method:     ", x$method, "\n", sep = "")
  cat(" CI method:  ", x$ci_method, "\n", sep = "")
  cat(
    " Times:      ",
    length(x$times),
    " evaluation time",
    if (length(x$times) == 1L) "" else "s",
    " (",
    paste(round(head(x$times, 5), 2), collapse = ", "),
    if (length(x$times) > 5L) ", ..." else "",
    ")\n",
    sep = ""
  )
  cat(" Patterns:   ", nrow(x$newdata), " covariate pattern",
      if (nrow(x$newdata) == 1L) "" else "s", "\n", sep = "")
  cat("\nAbsolute risk estimates:\n")
  print(x$estimates)
  invisible(x)
}

#' Tidy a matchatr_absolute_risk object
#'
#' Returns the estimates table as a `data.table` with one row per (covariate
#' pattern row, evaluation time), with columns `row`, `time`, `estimate`,
#' `ci_lower`, `ci_upper`.
#'
#' @param x A `matchatr_absolute_risk` object.
#' @param ... Unused.
#' @returns A `data.table`.
#' @family contrasts
#' @seealso [absolute_risk()]
#' @export
tidy.matchatr_absolute_risk <- function(x, ...) {
  x$estimates
}
