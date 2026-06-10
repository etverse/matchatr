#' Time-varying cumulative excess hazard from an additive NCC fit
#'
#' @description
#' Estimates the Aalen cumulative regression functions
#' B_j(t) = ∫₀ᵗ β_j(s) ds — the **time-varying excess cumulative hazard** for each
#' covariate — from a fitted IPW nested case-control additive-hazards
#' (`"ipw_aalen"`) design. Where `contrast(type = "excess")` reports one
#' time-*constant* excess hazard per covariate (the Lin-Ying model), `excess_risk()`
#' relaxes the constant-effect assumption and returns the full cumulative
#' regression function over time — the additive analogue of [absolute_risk()].
#'
#' @details
#' The fully nonparametric additive-hazards model (Aalen 1980) writes
#' λ(t | x) = β₀(t) + Σ_j β_j(t) x_j with time-varying coefficients. Its weighted
#' least-squares estimator gives the cumulative regression functions in increments
#' at each event time t_i,
#'
#'   dB̂(t_i) = (X̃ᵀ W X̃)⁻¹ X̃ᵀ W dN(t_i),
#'
#' where X̃ is the at-risk design matrix (an intercept column plus the covariates),
#' W = diag(w) the Samuelsen inclusion weights, and dN(t_i) the event indicator
#' over the risk set. B̂_j(t) = Σ_{t_i ≤ t} dB̂_j(t_i). Pointwise variances use the
#' Aalen (martingale) estimator
#'
#'   V̂(t) = Σ_{t_i ≤ t} (X̃ᵀ W X̃)⁻¹ {Σ_{i fails} w_i² x̃_i x̃_iᵀ} (X̃ᵀ W X̃)⁻¹,
#'
#' and the pointwise interval is B̂_j(t) ± z·sqrt(V̂_jj(t)) on the linear scale
#' (symmetric, possibly negative — an excess hazard is a rate difference). The
#' weighted estimator reuses controls via 1/π_j and is valid under nested
#' case-control sampling (Borgan & Langholz 1997). `timereg::aalen` (without
#' `const()`) computes the same `cum` and `var.cum`.
#'
#' The estimator is defined only while the weighted design X̃ᵀWX̃ is non-singular;
#' if the risk set shrinks below full column rank at a late event time, B̂(t) is
#' truncated there and a `matchatr_truncated_excess` warning reports the last
#' estimable time. The intercept (baseline cumulative hazard) is not reported —
#' only the covariate excess-hazard functions.
#'
#' @param fit A `matchatr_fit` returned by [matcha()] with
#'   `estimator = "ipw_aalen"`.
#' @param ... Unused; present for S3 consistency.
#' @returns A `matchatr_excess_risk` object. See `excess_risk.matchatr_fit` for
#'   the structure.
#' @examples
#' \dontrun{
#' ncc <- sample_ncc(cohort, time = "t", event = "d", m = 3, incl_prob = TRUE)
#' fit <- matcha(ncc, outcome = "d", exposure = "x",
#'               design = nested_cc(strata = "set", time = "t"),
#'               estimator = "ipw_aalen")
#' excess_risk(fit, times = c(1, 2, 3, 4))
#' }
#' @family contrasts
#' @seealso [matcha()], [contrast()], [absolute_risk()], [sample_ncc()]
#' @export
excess_risk <- function(fit, ...) UseMethod("excess_risk")

#' @rdname excess_risk
#'
#' @param times Non-empty numeric vector of evaluation times. Duplicates are
#'   dropped and times sorted; times before the first event return B̂ = 0, times
#'   after the last estimable event return the last cumulative value
#'   (step-function extrapolation).
#' @param conf_level Numeric in (0, 1). Pointwise confidence level. Default
#'   `0.95`.
#'
#' @returns A `matchatr_excess_risk` list with elements:
#'   - `$estimates`: a `data.table` with columns `term` (covariate),
#'     `time`, `estimate` (B̂_j(t)), `se`, `ci_lower`, `ci_upper` (symmetric Wald
#'     interval on the linear scale).
#'   - `$times`: the sorted evaluation times.
#'   - `$terms`: the reported covariate terms.
#'   - `$conf_level`, `$ci_method`, `$engine`, `$estimator`, `$n`.
#'
#' @export
excess_risk.matchatr_fit <- function(fit, times, conf_level = 0.95, ...) {
  # The cumulative excess hazard is the time-varying counterpart of the constant
  # additive-hazards contrast, so it is defined only for that engine.
  if (!identical(fit$engine, "ipw_aalen")) {
    rlang::abort(
      c(
        paste0(
          "`excess_risk()` is not implemented for the `",
          fit$engine,
          "` engine."
        ),
        i = paste0(
          "It reports the time-varying additive cumulative regression function, ",
          "so it needs an IPW nested case-control additive fit ",
          "(`estimator = \"ipw_aalen\"`)."
        )
      ),
      class = c("matchatr_not_implemented", "matchatr_error")
    )
  }
  if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times))) {
    rlang::abort(
      "`times` must be a non-empty numeric vector of finite values.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }

  time_col <- require_ipw_ncc_columns(fit, "ipw_aalen")
  dt <- ncc_ipw_analysis_data(fit)

  conf_terms <- if (is.null(fit$confounders)) {
    character(0)
  } else {
    attr(stats::terms(fit$confounders), "term.labels")
  }
  rhs <- stats::reformulate(c(fit$exposure, conf_terms))

  # Complete records only: the weighted LS solve needs the full design row.
  needed <- c(time_col, fit$outcome, all.vars(rhs))
  complete <- stats::complete.cases(
    do.call(cbind.data.frame, lapply(needed, function(v) dt[[v]]))
  )
  if (any(!complete)) {
    rlang::warn(
      c(
        paste0(
          sum(!complete),
          " row(s) with missing values were dropped from the excess-risk fit."
        ),
        i = "The cumulative excess hazard is estimated on the complete cases only."
      ),
      class = c("matchatr_dropped_rows", "matchatr_warning")
    )
  }
  dt <- dt[complete, , drop = FALSE]

  # Design with the intercept column (the additive baseline β₀(t)); the
  # covariate columns are the time-varying excess hazards.
  X <- stats::model.matrix(rhs, data = dt)

  ar <- aalen_cumulative(
    time = dt[[time_col]],
    status = as.integer(dt[[fit$outcome]]),
    X = X,
    w = dt[["ipw_weight"]]
  )
  if (ar$truncated) {
    rlang::warn(
      c(
        paste0(
          "The weighted additive design became singular; the cumulative excess ",
          "hazard is truncated at t = ",
          format(ar$last_time, digits = 4),
          "."
        ),
        i = "Later event times have too few at-risk subjects to identify B(t)."
      ),
      class = c("matchatr_truncated_excess", "matchatr_warning")
    )
  }

  assemble_excess_risk(
    fit = fit,
    cumulative = ar,
    times = times,
    conf_level = conf_level,
    n = nrow(X)
  )
}

#' Weighted Aalen cumulative regression functions with martingale variance
#'
#' Hand-rolled weighted least-squares estimator for the Aalen additive-hazards
#' model with time-varying coefficients (Aalen 1980). It backs [excess_risk()],
#' where the weights are the Samuelsen inclusion weights. The point estimate and
#' the martingale (Aalen) pointwise variance reproduce `timereg::aalen`'s `cum`
#' and `var.cum`.
#'
#' @details
#' At each unique event time t_i the increment of the cumulative regression
#' function is dB̂(t_i) = (X̃ᵀWX̃)⁻¹ X̃ᵀW dN(t_i) and the increment of the
#' martingale variance is (X̃ᵀWX̃)⁻¹ {Σ_{i fails} w_i² x̃_i x̃_iᵀ} (X̃ᵀWX̃)⁻¹, both
#' over the risk set {j : T_j ≥ t_i}. The estimator is only defined while X̃ᵀWX̃ is
#' invertible; accumulation stops at the first event time where it is singular.
#'
#' @param time Numeric vector of exit / event times.
#' @param status Integer 0/1 event indicator (1 = event).
#' @param X Numeric design matrix with an intercept column followed by the
#'   covariates (`colnames(X)` name the returned series).
#' @param w Numeric vector of observation (inclusion) weights.
#' @returns A list with `times` (estimable event times), `B` (cumulative
#'   regression matrix, one row per time, one column per design column), `Vdiag`
#'   (the matching pointwise martingale variances), `terms` (`colnames(X)`),
#'   `truncated` (logical), and `last_time` (last estimable time).
#' @family contrasts
#' @seealso [excess_risk()], [timereg::aalen()] (the external oracle)
#' @noRd
aalen_cumulative <- function(time, status, X, w) {
  X <- as.matrix(X)
  p <- ncol(X)
  cn <- colnames(X)
  ev_times <- sort(unique(time[status == 1L]))
  K <- length(ev_times)

  B <- matrix(NA_real_, K, p)
  Vdiag <- matrix(NA_real_, K, p)
  accB <- numeric(p)
  accV <- matrix(0, p, p)
  last_ok <- 0L
  for (i in seq_len(K)) {
    tk <- ev_times[i]
    atrisk <- time >= tk
    Xr <- X[atrisk, , drop = FALSE]
    wr <- w[atrisk]
    xtwx_inv <- tryCatch(solve(crossprod(Xr, wr * Xr)), error = function(e) {
      NULL
    })
    if (is.null(xtwx_inv)) {
      # Risk set below full column rank: B(t) is no longer identified.
      break
    }
    fail <- atrisk & time == tk & status == 1L
    Xf <- X[fail, , drop = FALSE]
    wf <- w[fail]
    # dB = (X̃'WX̃)^-1 Σ_fail w x̃;  dN = 1 only at the failing rows.
    accB <- accB + as.numeric(xtwx_inv %*% crossprod(Xf, wf))
    B[i, ] <- accB
    # Martingale variance increment with the squared weights.
    accV <- accV + xtwx_inv %*% crossprod(Xf, wf^2 * Xf) %*% xtwx_inv
    Vdiag[i, ] <- diag(accV)
    last_ok <- i
  }

  keep <- seq_len(last_ok)
  list(
    times = ev_times[keep],
    B = B[keep, , drop = FALSE],
    Vdiag = Vdiag[keep, , drop = FALSE],
    terms = cn,
    truncated = last_ok < K,
    last_time = if (last_ok > 0L) ev_times[last_ok] else NA_real_
  )
}

#' Assemble a matchatr_excess_risk result from cumulative regression functions
#'
#' Step-interpolates the weighted Aalen cumulative regression functions and their
#' martingale variances to the requested evaluation times, attaches a symmetric
#' pointwise Wald interval, and wraps the per-(covariate, time) table in the
#' `matchatr_excess_risk` S3 structure. The intercept (baseline cumulative hazard)
#' column is dropped — only the covariate excess-hazard functions are reported.
#'
#' @param fit A `matchatr_fit` supplying the `estimator` / `engine` labels.
#' @param cumulative The list returned by `aalen_cumulative()`.
#' @param times Numeric vector of evaluation times (sorted, de-duplicated here).
#' @param conf_level Numeric confidence level.
#' @param n Integer analysis sample size recorded on the result.
#' @returns A `matchatr_excess_risk` object.
#' @family contrasts
#' @noRd
assemble_excess_risk <- function(fit, cumulative, times, conf_level, n) {
  times <- sort(unique(times))
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Covariate terms only (the additive baseline β₀(t) is the intercept column).
  report_cols <- which(cumulative$terms != "(Intercept)")
  term_names <- cumulative$terms[report_cols]

  # A t = 0 fence so times before the first event map to B = 0 (and variance 0),
  # not the first event's value (step-function, rule = 2 extrapolates the tail).
  knot_t <- c(0, cumulative$times)
  rows <- vector("list", length(term_names) * length(times))
  k <- 0L
  for (col in report_cols) {
    b_knot <- c(0, cumulative$B[, col])
    v_knot <- c(0, cumulative$Vdiag[, col])
    b_t <- stats::approx(
      knot_t,
      b_knot,
      xout = times,
      method = "constant",
      f = 0,
      rule = 2
    )$y
    v_t <- stats::approx(
      knot_t,
      v_knot,
      xout = times,
      method = "constant",
      f = 0,
      rule = 2
    )$y
    se_t <- sqrt(pmax(0, v_t))
    for (j in seq_along(times)) {
      rows[[k <- k + 1L]] <- list(
        term = cumulative$terms[col],
        time = times[j],
        estimate = b_t[j],
        se = se_t[j],
        ci_lower = b_t[j] - z_crit * se_t[j],
        ci_upper = b_t[j] + z_crit * se_t[j]
      )
    }
  }

  structure(
    list(
      estimates = data.table::rbindlist(rows),
      times = times,
      terms = term_names,
      conf_level = conf_level,
      ci_method = "Aalen martingale",
      engine = fit$engine,
      estimator = fit$estimator,
      n = n
    ),
    class = c("matchatr_excess_risk", "matchatr")
  )
}

#' Print a matchatr_excess_risk object
#'
#' Displays a compact summary — the engine, variance method, evaluation times,
#' and covariate terms — followed by the cumulative excess-hazard table.
#'
#' @param x A `matchatr_excess_risk` object from [excess_risk()].
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @family contrasts
#' @seealso [excess_risk()]
#' @export
print.matchatr_excess_risk <- function(x, ...) {
  cat("<matchatr_excess_risk>\n")
  cat(" Engine:     ", x$engine, "\n", sep = "")
  cat(" Estimand:   cumulative excess hazard B_j(t)\n")
  cat(" CI method:  ", x$ci_method, "\n", sep = "")
  cat(
    " Times:      ",
    length(x$times),
    " evaluation time",
    if (length(x$times) == 1L) "" else "s",
    "\n",
    sep = ""
  )
  cat(" Terms:      ", paste(x$terms, collapse = ", "), "\n", sep = "")
  cat("\nCumulative excess hazard:\n")
  print(x$estimates)
  invisible(x)
}

#' Tidy a matchatr_excess_risk object
#'
#' Returns the estimates table as a `data.table` with one row per (covariate
#' term, evaluation time): columns `term`, `time`, `estimate`, `se`, `ci_lower`,
#' `ci_upper`.
#'
#' @param x A `matchatr_excess_risk` object.
#' @param ... Unused.
#' @returns A `data.table`.
#' @family contrasts
#' @seealso [excess_risk()]
#' @export
tidy.matchatr_excess_risk <- function(x, ...) {
  x$estimates
}
