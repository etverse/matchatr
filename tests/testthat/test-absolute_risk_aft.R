# Tests for absolute_risk_aft.R — parametric Weibull AFT cumulative incidence
# F_x(t) from IPW nested case-control (ipw_aft) fits.
#
# Oracles:
#   - predict.survreg(type = "quantile"): the survival package's own inverse CDF.
#     Feeding our F̂ estimates back as probabilities must return the evaluation
#     times (a round-trip through survival's qsurvreg, independent of our forward
#     (log t − η)/σ formula). Machine precision.
#   - numDeriv::grad: the analytic delta-method gradient of ξ = (log t − η)/σ over
#     θ = (β, log σ) reproduced numerically, so the reconstructed estimate + CI
#     match absolute_risk() to machine precision (validates the gradient and the
#     whole CI pipeline).
#   - full-cohort survival::survreg: the NCC subsample F̂_x(t) recovers the
#     full-cohort Weibull AFT curve within sampling tolerance.
# Truth DGP: a Weibull AFT cohort, F_x(t) = 1 − exp(−(t·exp(−η))^(1/σ)) known in
#   closed form.
#
# No Python / delicatessen oracle: F̂_x(t) is a deterministic transform of the
# fitted survreg parameters (a parametric survival curve), not an M-estimator
# delicatessen stacks; the predict.survreg round-trip and the numDeriv gradient
# reconstruction are the independent checks.

# -- Shared DGP ----------------------------------------------------------------

# Weibull AFT cohort with a known linear predictor and shape. Generated as
# T = exp(η) · Weibull(shape, 1), so survreg(dist = "weibull") recovers
# (β₀, β_x, β_z) and σ = 1/shape, and the survival curve is analytic:
#   F(t | x) = 1 − exp(−(t · exp(−η))^shape),  η = β₀ + β_x x + β_z z.
make_aft_ar_cohort <- function(
  n = 4000L,
  beta0 = 0.5,
  beta_x = log(0.6),
  beta_z = 0.3,
  shape = 1.3,
  q_tau = 0.7,
  seed = 505L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    eta <- beta0 + beta_x * x + beta_z * z
    tt <- exp(eta) * stats::rweibull(n, shape = shape, scale = 1)
    tau <- stats::quantile(tt, q_tau)
    cohort <- data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z
    )
    attr(cohort, "truth") <- list(
      beta0 = beta0,
      beta_x = beta_x,
      beta_z = beta_z,
      shape = shape
    )
    cohort
  })
}

# Weibull AFT cohort with a 3-level factor confounder, exercising factor
# contrasts in the absolute-risk linear predictor / gradient.
make_aft_ar_cohort_factor <- function(n = 6000L, seed = 606L) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.45)
    z <- stats::rnorm(n)
    g <- factor(sample(c("a", "b", "c"), n, TRUE, c(0.5, 0.3, 0.2)))
    bg <- c(a = 0, b = 0.4, c = -0.5)
    eta <- 0.4 + log(0.7) * x + 0.25 * z + bg[as.character(g)]
    tt <- exp(eta) * stats::rweibull(n, shape = 1.2, scale = 1)
    tau <- stats::quantile(tt, 0.7)
    data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z,
      g = g
    )
  })
}

# Fit the ipw_aft engine on an NCC sample drawn from such a cohort.
make_aft_ar_fit <- function(
  cohort,
  confounders = ~z,
  m = 3L,
  seed = 9L,
  dist = NULL
) {
  ncc <- withr::with_seed(
    seed,
    sample_ncc(cohort, "t", "d", m = m, incl_prob = TRUE)
  )
  matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = confounders,
    estimator = "ipw_aft",
    dist = dist
  )
}

# The survreg baseline error CDF G: F(t | x) = G((log t − η)/σ).
aft_cdf <- function(dist) {
  switch(
    dist,
    weibull = ,
    exponential = function(z) 1 - exp(-exp(z)),
    lognormal = stats::pnorm,
    loglogistic = stats::plogis
  )
}

# Independent reconstruction of the cumulative incidence and its CI via a
# numerically differentiated z(θ) gradient and the baseline error CDF. The
# exponential fixes σ = 1 (no log-scale parameter), so θ omits it there.
recon_aft_ar <- function(fit, nd_row, times, conf_level = 0.95) {
  beta <- stats::coef(fit$model)
  vmat <- stats::vcov(fit$model)
  sigma <- fit$model$scale
  has_scale <- ncol(vmat) > length(beta)
  G <- aft_cdf(fit$model$dist)
  mt <- stats::delete.response(stats::terms(fit$model))
  mm <- stats::model.matrix(
    mt,
    stats::model.frame(mt, data = nd_row, xlev = fit$model$xlevels)
  )[, names(beta), drop = FALSE]
  theta0 <- if (has_scale) c(beta, log(sigma)) else beta
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)
  do.call(
    rbind,
    lapply(times, function(tt) {
      z <- (log(tt) - as.numeric(mm %*% beta)) / sigma
      zfun <- function(th) {
        if (has_scale) {
          p <- length(th)
          (log(tt) - sum(mm * th[-p])) / exp(th[p])
        } else {
          (log(tt) - sum(mm * th)) / sigma
        }
      }
      g <- numDeriv::grad(zfun, theta0)
      se <- sqrt(as.numeric(t(g) %*% vmat %*% g))
      c(
        estimate = G(z),
        ci_lower = G(z - z_crit * se),
        ci_upper = G(z + z_crit * se)
      )
    })
  )
}

# -- Structural checks ---------------------------------------------------------

test_that("absolute_risk(ipw_aft) returns correct S3 class and structure", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = c(1L, 0L), z = c(0, 0))
  ar <- absolute_risk(fit, newdata = nd, times = c(1, 2, 4))

  expect_s3_class(ar, "matchatr_absolute_risk")
  expect_s3_class(ar, "matchatr")
  expect_true(data.table::is.data.table(ar$estimates))
  expect_named(
    ar$estimates,
    c("row", "time", "estimate", "ci_lower", "ci_upper")
  )
  expect_identical(ar$engine, "ipw_aft")
  expect_identical(ar$method, "IPW AFT (weibull)")
})

test_that("absolute_risk(ipw_aft) returns one row per (newdata row, time)", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = c(1L, 0L), z = c(0.5, -0.5))
  times <- c(1, 2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = times)

  expect_equal(nrow(ar$estimates), nrow(nd) * length(times))
  expect_equal(sort(unique(ar$estimates$row)), 1:2)
})

test_that("absolute_risk(ipw_aft) estimates in [0, 1], CIs ordered and monotone", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(1, 2, 4, 6))
  e <- ar$estimates

  expect_true(all(e$estimate >= 0 & e$estimate <= 1))
  expect_true(all(e$ci_lower >= 0 & e$ci_upper <= 1))
  expect_true(all(e$ci_lower <= e$estimate + 1e-10))
  expect_true(all(e$estimate <= e$ci_upper + 1e-10))
  # Parametric Weibull F is strictly increasing in t.
  expect_true(all(diff(e$estimate) > 0))
})

test_that("tidy/print on an ipw_aft absolute-risk object behave", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(2, 5))
  expect_identical(tidy(ar), ar$estimates)
  expect_invisible(print(ar))
})

# -- Exact point oracle: predict.survreg quantile round-trip -------------------

test_that("AFT F_x(t) round-trips through predict.survreg quantiles", {
  skip_if_not_installed("survival")
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = c(1L, 0L, 1L), z = c(0, 0.5, -1))
  t_eval <- c(1, 2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  # Feed our F̂ back as probabilities; survreg's inverse CDF must return the
  # evaluation times. This validates the forward F̂ against survival's own
  # quantile machinery, not merely a rearrangement of the same formula.
  for (r in seq_len(nrow(nd))) {
    p_hat <- ar$estimates$estimate[ar$estimates$row == r]
    q <- stats::predict(
      fit$model,
      newdata = nd[r, , drop = FALSE],
      type = "quantile",
      p = p_hat
    )
    expect_equal(as.numeric(q), t_eval, tolerance = 1e-7)
  }
})

# -- Exact variance / CI: independent numDeriv reconstruction ------------------

test_that("AFT absolute-risk estimate and CI match a numDeriv reconstruction", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = c(1L, 0L), z = c(0.3, -0.7))
  t_eval <- c(1.5, 3, 5)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  for (r in seq_len(nrow(nd))) {
    recon <- recon_aft_ar(fit, nd[r, , drop = FALSE], t_eval)
    sub <- ar$estimates[ar$estimates$row == r, ]
    expect_equal(sub$estimate, unname(recon[, "estimate"]), tolerance = 1e-7)
    expect_equal(sub$ci_lower, unname(recon[, "ci_lower"]), tolerance = 1e-7)
    expect_equal(sub$ci_upper, unname(recon[, "ci_upper"]), tolerance = 1e-7)
  }
})

# -- Non-Weibull baselines: each distribution's survival curve ----------------

test_that("AFT absolute risk is exact across all four baseline distributions", {
  skip_if_not_installed("survival")
  cohort <- make_aft_ar_cohort()
  nd <- data.frame(x = c(1L, 0L), z = c(0, 0.5))
  t_eval <- c(1, 2, 3, 4)

  for (dd in c("weibull", "exponential", "lognormal", "loglogistic")) {
    fit <- make_aft_ar_fit(cohort, dist = dd)
    ar <- absolute_risk(fit, newdata = nd, times = t_eval)
    expect_identical(ar$method, paste0("IPW AFT (", dd, ")"))

    for (r in seq_len(nrow(nd))) {
      # Independent point oracle: round-trip through survreg's own inverse CDF.
      p_hat <- ar$estimates$estimate[ar$estimates$row == r]
      q <- stats::predict(
        fit$model,
        newdata = nd[r, , drop = FALSE],
        type = "quantile",
        p = p_hat
      )
      expect_equal(as.numeric(q), t_eval, tolerance = 1e-7)
      # Estimate + CI vs the numDeriv reconstruction with the dist's error CDF.
      recon <- recon_aft_ar(fit, nd[r, , drop = FALSE], t_eval)
      sub <- ar$estimates[ar$estimates$row == r, ]
      expect_equal(sub$estimate, unname(recon[, "estimate"]), tolerance = 1e-7)
      expect_equal(sub$ci_lower, unname(recon[, "ci_lower"]), tolerance = 1e-7)
      expect_equal(sub$ci_upper, unname(recon[, "ci_upper"]), tolerance = 1e-7)
    }
  }
})

# -- Complicated design: factor confounder ------------------------------------

test_that("AFT absolute risk round-trips under a 3-level factor confounder", {
  skip_if_not_installed("survival")
  cohort <- make_aft_ar_cohort_factor()
  fit <- make_aft_ar_fit(cohort, confounders = ~ z + g, m = 4L, seed = 11L)

  nd <- data.frame(
    x = c(1L, 0L, 1L, 1L),
    z = c(0, 0, 1, -1),
    g = factor(c("a", "a", "b", "c"), levels = c("a", "b", "c"))
  )
  # Moderate times keep F̂ away from 1, where survreg's inverse-CDF (quantile) is
  # ill-conditioned (the density vanishes in the far tail, so a tiny change in p
  # maps to a large change in t). The forward F̂ is exact at all times — the
  # numDeriv reconstruction below validates it including the extreme tail.
  t_eval <- c(1, 2, 3, 4)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  for (r in seq_len(nrow(nd))) {
    p_hat <- ar$estimates$estimate[ar$estimates$row == r]
    q <- stats::predict(
      fit$model,
      newdata = nd[r, , drop = FALSE],
      type = "quantile",
      p = p_hat
    )
    expect_equal(as.numeric(q), t_eval, tolerance = 1e-7)
    # And the estimate + CI match the numDeriv reconstruction with factor
    # contrasts, here over a wider time grid that includes the extreme tail.
    t_full <- c(1, 2, 4, 6, 9)
    ar_full <- absolute_risk(
      fit,
      newdata = nd[r, , drop = FALSE],
      times = t_full
    )
    recon <- recon_aft_ar(fit, nd[r, , drop = FALSE], t_full)
    expect_equal(
      ar_full$estimates$estimate,
      unname(recon[, "estimate"]),
      tolerance = 1e-7
    )
    expect_equal(
      ar_full$estimates$ci_lower,
      unname(recon[, "ci_lower"]),
      tolerance = 1e-7
    )
    expect_equal(
      ar_full$estimates$ci_upper,
      unname(recon[, "ci_upper"]),
      tolerance = 1e-7
    )
  }
})

# -- Weight independence: GLM working-model weights ---------------------------

test_that("AFT absolute risk is exact under GLM working-model weights", {
  skip_if_not_installed("survival")
  # The absolute-risk machinery reads (β̂, σ̂, V) off the fit; it is agnostic to
  # how the inclusion weights were built. GLM working-model weights change the
  # fit but the parametric F̂ stays exact against survreg's own quantiles.
  cohort <- make_aft_ar_cohort(seed = 808L)
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d", m = 4L, incl_prob = TRUE)
  )
  ncc_glm <- compute_ncc_weights(
    ncc,
    cohort = cohort,
    method = "glm",
    selection_formula = ~ risk_time + z,
    time = "t"
  )
  fit <- matcha(
    ncc_glm,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aft"
  )
  nd <- data.frame(x = c(1L, 0L), z = c(0, 0))
  t_eval <- c(2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  for (r in seq_len(nrow(nd))) {
    p_hat <- ar$estimates$estimate[ar$estimates$row == r]
    q <- stats::predict(
      fit$model,
      newdata = nd[r, , drop = FALSE],
      type = "quantile",
      p = p_hat
    )
    expect_equal(as.numeric(q), t_eval, tolerance = 1e-7)
  }
})

# -- Edge case: an evaluation time at / before the origin ---------------------

test_that("absolute_risk(ipw_aft) is 0 at t <= 0 without warning", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd <- data.frame(x = 1L, z = 0)
  expect_no_warning(
    ar <- absolute_risk(fit, newdata = nd, times = c(0, 1, 3))
  )
  z0 <- ar$estimates[ar$estimates$time == 0, ]
  expect_equal(z0$estimate, 0)
  expect_equal(z0$ci_lower, 0)
  expect_equal(z0$ci_upper, 0)
  expect_true(all(is.finite(ar$estimates$estimate)))
})

# -- Oracle: full-cohort survreg ----------------------------------------------

test_that("absolute_risk(ipw_aft) agrees with full-cohort survreg (sampling tol)", {
  skip_if_not_installed("survival")
  cohort <- make_aft_ar_cohort(n = 4000L, seed = 303L)
  fit <- make_aft_ar_fit(cohort, m = 3L, seed = 13L)

  nd <- data.frame(x = 1L, z = 0)
  t_eval <- c(2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_ipw <- ar$estimates$estimate

  # Full-cohort reference: the parametric Weibull F̂ from a full-data survreg.
  full <- survival::survreg(
    survival::Surv(t, d) ~ x + z,
    data = cohort,
    dist = "weibull"
  )
  eta_full <- as.numeric(stats::predict(full, newdata = nd, type = "lp"))
  f_full <- 1 - exp(-exp((log(t_eval) - eta_full) / full$scale))

  # The NCC subsample (m = 3) reuses controls to estimate the same parametric
  # curve; discrepancies up to 0.05 in F are within sampling variability.
  expect_true(
    all(abs(f_ipw - f_full) < 0.05),
    info = paste("Max discrepancy:", round(max(abs(f_ipw - f_full)), 4))
  )
  expect_true(
    all(
      ar$estimates$ci_lower <= f_full + 0.02 &
        f_full <= ar$estimates$ci_upper + 0.02
    ),
    info = "Full-cohort F_x(t) should lie inside the IPW NCC AFT CI"
  )
})

# -- Truth DGP: Weibull AFT with known F_x(t) ---------------------------------

test_that("absolute_risk(ipw_aft) recovers known F_x(t) in Weibull DGP", {
  cohort <- make_aft_ar_cohort(n = 4000L, seed = 404L)
  truth <- attr(cohort, "truth")
  fit <- make_aft_ar_fit(cohort, m = 3L, seed = 21L)

  # Analytical truth at x = 1, z = 0: η = β₀ + β_x; shape = 1/σ.
  t_eval <- c(2, 4, 6)
  eta1 <- truth$beta0 + truth$beta_x
  f_true <- 1 - exp(-((t_eval * exp(-eta1))^truth$shape))

  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_hat <- ar$estimates$estimate

  expect_true(
    all(abs(f_hat - f_true) < 0.06),
    info = paste("Max discrepancy:", round(max(abs(f_hat - f_true)), 4))
  )
  expect_true(
    all(ar$estimates$ci_lower <= f_true & f_true <= ar$estimates$ci_upper),
    info = "95% CI should cover the analytical truth"
  )
})

# -- Rejections ---------------------------------------------------------------

test_that("absolute_risk(ipw_aft) rejects mismatched newdata columns", {
  fit <- make_aft_ar_fit(make_aft_ar_cohort())
  nd_wrong <- data.frame(x = 1L, foo = 99)
  expect_error(
    absolute_risk(fit, newdata = nd_wrong, times = 4),
    class = "matchatr_bad_input"
  )
})

test_that("absolute_risk() rejects the additive (ipw_aalen) engine", {
  # The constant additive-hazards model reports a scalar excess hazard; it has no
  # survival-curve verb, so absolute_risk() is not implemented for it.
  cohort <- make_aft_ar_cohort(n = 800L, seed = 1L)
  ncc <- withr::with_seed(
    2L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aalen"
  )
  expect_error(
    absolute_risk(fit, newdata = data.frame(x = 1L, z = 0), times = 4),
    class = "matchatr_not_implemented"
  )
})
