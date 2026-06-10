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
#' martingale variance is (X̃ᵀWX̃)⁻¹ (Σ over failures w_i² x̃_i x̃_iᵀ) (X̃ᵀWX̃)⁻¹, both
#' over the risk set (j : T_j ≥ t_i). The estimator is only defined while X̃ᵀWX̃ is
#' invertible; because the risk set is nested-decreasing in t, a singular Gram
#' matrix never recovers, so accumulation stops at the first singular event time.
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
