# Oracle wrapper for CCW-TMLE: `tmle::tmle()` with the Rose & van der Laan
# case-control weights passed as observation weights. cvQinit = FALSE forces a
# plain glm initial outcome fit (matching matchatr's), gbound = 0.025 matches the
# matchatr default, and Qform / gform pin the same parametric nuisance models, so
# the two implementations target the same Q̄* and should agree closely (the risk
# difference essentially exactly). The single confounder column is named `w`.
# tmle prints a non-fatal internal "object 'w' not found" line to stderr during
# its auxiliary variance step; it does not affect the returned estimates (which
# this wrapper returns) and is not an R condition, so nothing is suppressed.
tmle_ccw_oracle <- function(cc, q0) {
  y01 <- as.integer(cc$case)
  n1 <- sum(y01 == 1L)
  n0 <- sum(y01 == 0L)
  n <- n1 + n0
  wt <- ifelse(y01 == 1L, q0 / (n1 / n), (1 - q0) / (n0 / n))
  fit <- tmle::tmle(
    Y = cc$case,
    A = cc$x,
    W = cc[, "w", drop = FALSE],
    Qform = "Y ~ A + w",
    gform = "A ~ w",
    family = "binomial",
    obsWeights = wt,
    cvQinit = FALSE,
    gbound = 0.025
  )
  fit$estimates
}
