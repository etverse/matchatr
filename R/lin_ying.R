#' Weighted Lin-Ying constant additive-hazards estimator with robust variance
#'
#' Semiparametric additive-hazards model with time-constant covariate effects
#' (Lin & Ying 1994, Biometrika 81(1)), fit by weighted estimating equations. It
#' backs the IPW nested case-control additive engine (`fit_ipw_aalen()`), where
#' the weights are the Samuelsen inclusion weights.
#'
#' @details
#' The model is λ(t | Z) = λ₀(t) + γᵀ Z with γ time-constant, so γ_j is the
#' excess hazard (additive rate difference) for covariate j. With weights w_i and
#' time-fixed covariates the estimating equation is
#'
#'   U(γ) = Σ_i w_i ∫ {Z_i − Z̄(t)} {dN_i(t) − Y_i(t) γᵀZ_i dt} = 0,
#'
#' where Z̄(t) = Σ_j w_j Y_j(t) Z_j / Σ_j w_j Y_j(t) is the weighted risk-set mean
#' and Y_i, N_i are the at-risk and counting processes. The closed-form solution
#' is γ̂ = A⁻¹ B with
#'
#'   A = ∫ {S²(t) − S¹(t)S¹(t)ᵀ / S⁰(t)} dt,   B = Σ_i w_i δ_i {Z_i − Z̄(T_i)},
#'
#' S⁰/S¹/S² the weighted risk-set moments; the integrand of A is piecewise
#' constant between distinct exit times, so A is a finite sum over them.
#'
#' The variance is the robust sandwich A⁻¹ (Σ_i η̂_i η̂_iᵀ) A⁻¹ from the
#' martingale-residual influence contributions
#'
#'   η̂_i = w_i ( δ_i {Z_i − Z̄(T_i)} − Σ_{t_k ≤ T_i} {Z_i − Z̄(t_k)} a_k − C_i γ̂ ),
#'
#' with a_k = (Σ events at t_k w) / S⁰(t_k) the weighted additive-baseline jump
#' and C_i = Σ_{t_l ≤ T_i} {Z_i − Z̄_l}{Z_i − Z̄_l}ᵀ Δ_l the subject's A-integrand
#' (so Σ_i w_i C_i = A and Σ_i η̂_i = U(γ̂) = 0), the appropriate (slightly
#' conservative) variance under the reused-control IPW weights.
#'
#' Risk-set moments are accumulated by reverse cumulative sums over the
#' exit-time-sorted sample (O(n p²)); ties at a shared exit time share the risk
#' set {j : T_j ≥ t}.
#'
#' @param time Numeric vector of exit / event times.
#' @param status Integer 0/1 event indicator (1 = event).
#' @param Z Numeric matrix of time-fixed covariates (one column per coefficient;
#'   `colnames(Z)` name the returned `gamma` / `robvar`). No intercept column.
#' @param w Numeric vector of observation (inclusion) weights.
#' @returns A list with `gamma` (named numeric vector of excess hazards) and
#'   `robvar` (their robust sandwich variance matrix, dimnamed by `colnames(Z)`).
#' @family estimators
#' @seealso `fit_ipw_aalen()`, [timereg::aalen()] (the external oracle)
#' @noRd
lin_ying_additive <- function(time, status, Z, w) {
  Z <- as.matrix(Z)
  p <- ncol(Z)
  cn <- colnames(Z)

  # Sort by exit time; risk set {j: T_j >= t} is then a suffix of the order.
  o <- order(time)
  time <- time[o]
  status <- status[o]
  Z <- Z[o, , drop = FALSE]
  w <- w[o]
  n <- length(time)

  # Suffix (reverse-cumulative) weighted risk-set moments at every row.
  suffix <- function(v) rev(cumsum(rev(v)))
  rc_S0 <- suffix(w)
  wZ <- w * Z
  rc_S1 <- apply(wZ, 2L, suffix) # n x p
  # S2 suffix sums: one reverse-cumsum per (a, b) entry of w Z Zᵀ.
  rc_S2 <- array(0, dim = c(n, p, p))
  for (a in seq_len(p)) {
    for (b in a:p) {
      s <- suffix(w * Z[, a] * Z[, b])
      rc_S2[, a, b] <- s
      rc_S2[, b, a] <- s
    }
  }

  # Distinct exit times and the sorted-array index where each risk set starts.
  ut <- unique(time)
  K <- length(ut)
  start <- match(ut, time) # first row with time == ut_l (suffix start)
  Delta <- diff(c(0, ut)) # interval lengths, t_(0) = 0

  S0 <- rc_S0[start]
  S1 <- rc_S1[start, , drop = FALSE]
  Zbar <- S1 / S0 # K x p weighted risk-set means

  # A = Σ_l {S²_l − S¹_l S¹_lᵀ / S⁰_l} Δ_l.
  A <- matrix(0, p, p)
  for (l in seq_len(K)) {
    S2l <- matrix(rc_S2[start[l], , ], p, p)
    A <- A + (S2l - tcrossprod(S1[l, ]) / S0[l]) * Delta[l]
  }

  # B = Σ_{events} w_i {Z_i − Z̄(T_i)}.
  tidx <- match(time, ut)
  centered <- Z - Zbar[tidx, , drop = FALSE]
  B <- colSums(w * status * centered)
  # A singular design (a constant covariate, or one collinear with the others)
  # has no solution; surface it as a classed error rather than a raw LAPACK fault.
  Ainv <- tryCatch(
    solve(A),
    error = function(e) {
      rlang::abort(
        c(
          "The additive-hazards model has no estimable solution.",
          i = "A covariate is constant or collinear with the others, so the weighted design is singular."
        ),
        class = c("matchatr_unestimable_exposure", "matchatr_error")
      )
    }
  )
  gamma <- as.numeric(Ainv %*% B)

  # Robust sandwich via the influence contributions η̂_i.
  ev_w <- as.numeric(tapply(
    w * status,
    factor(tidx, levels = seq_len(K)),
    sum
  ))
  ev_w[is.na(ev_w)] <- 0
  a_k <- ev_w / S0 # weighted additive-baseline jumps
  c_l <- as.numeric(Zbar %*% gamma) # Z̄_lᵀ γ̂

  # Cumulative-over-interval quantities evaluated at each subject's exit rank,
  # so the per-subject jump and C_i γ̂ terms are O(n p) not O(n²).
  P <- cumsum(a_k)
  Q <- apply(Zbar * a_k, 2L, cumsum) # K x p
  Ecum <- cumsum(c_l * Delta)
  Fcum <- apply(Zbar * Delta, 2L, cumsum) # K x p
  Hcum <- apply(Zbar * c_l * Delta, 2L, cumsum) # K x p
  Dcum <- cumsum(Delta) # = ut
  b_i <- as.numeric(Z %*% gamma) # Z_iᵀ γ̂

  li <- tidx
  # jump_i = Z_i P(T_i) − Q(T_i);  C_i γ̂ = Z_i b_i D − Z_i E − b_i F + H.
  jump <- Z * P[li] - Q[li, , drop = FALSE]
  Cgamma <- Z *
    (b_i * Dcum[li]) -
    Z * Ecum[li] -
    b_i * Fcum[li, , drop = FALSE] +
    Hcum[li, , drop = FALSE]
  eta <- w * (status * centered - jump - Cgamma)

  robvar <- Ainv %*% crossprod(eta) %*% Ainv
  dimnames(robvar) <- list(cn, cn)
  names(gamma) <- cn
  list(gamma = gamma, robvar = robvar)
}
