# Tests for absolute_risk_ncc.R — IPW Breslow cumulative incidence F_x(t) from
# IPW nested case-control (ipw_cox) fits.
#
# Oracles:
#   - survival::survfit(fit$model): the hand-rolled weighted IPW Breslow F_x(t)
#     must equal survival's own weighted Breslow on the same coxph to machine
#     precision (an exact internal cross-check of the Breslow arithmetic).
#   - full-cohort survival::survfit(coxph): the NCC subsample F_x(t) recovers the
#     full-cohort cumulative incidence within sampling tolerance.
# Truth DGP: exponential survival with known beta and baseline, so the analytical
#   F_x(t) = 1 - exp(-exp(beta_x x + beta_z z) * lambda0 * t) is known in closed form.
#
# No Python oracle: the IPW Breslow estimator is a plug-in step function, not an
# M-estimator, and is not covered by delicatessen.

# -- Shared DGP ----------------------------------------------------------------

# Exponential cohort with a known Cox log-HR and baseline hazard. The observed
# time and event come from a SINGLE latent draw (administrative censoring at tau),
# so censoring is internally consistent.
make_exp_ncc_cohort <- function(
  n = 3000L,
  beta_x = log(2),
  beta_z = 0.3,
  lambda0 = 0.08,
  tau = 8,
  seed = 202L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    rate <- lambda0 * exp(beta_x * x + beta_z * z)
    true_t <- stats::rexp(n, rate)
    d <- as.integer(true_t <= tau)
    t_obs <- pmin(true_t, tau)
    cohort <- data.frame(id = seq_len(n), t = t_obs, d = d, x = x, z = z)
    attr(cohort, "truth") <- c(
      beta_x = beta_x,
      beta_z = beta_z,
      lambda0 = lambda0
    )
    cohort
  })
}

# Fit the ipw_cox engine on an NCC sample drawn from such a cohort.
make_ipw_ar_fit <- function(cohort, m = 3L, seed = 9L) {
  ncc <- withr::with_seed(
    seed,
    sample_ncc(cohort, "t", "d", m = m, incl_prob = TRUE)
  )
  matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_cox"
  )
}

# -- Structural checks ---------------------------------------------------------

test_that("absolute_risk(ipw_cox) returns correct S3 class and structure", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = c(1L, 0L), z = c(0, 0))
  ar <- absolute_risk(fit, newdata = nd, times = c(1, 2, 4))

  expect_s3_class(ar, "matchatr_absolute_risk")
  expect_s3_class(ar, "matchatr")
  expect_true(data.table::is.data.table(ar$estimates))
  expect_named(
    ar$estimates,
    c("row", "time", "estimate", "ci_lower", "ci_upper")
  )
  expect_identical(ar$engine, "ipw_cox")
  expect_identical(ar$method, "IPW")
})

test_that("absolute_risk(ipw_cox) returns one row per (newdata row, time)", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = c(1L, 0L), z = c(0.5, -0.5))
  times <- c(1, 2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = times)

  expect_equal(nrow(ar$estimates), nrow(nd) * length(times))
  expect_equal(sort(unique(ar$estimates$row)), 1:2)
})

test_that("absolute_risk(ipw_cox) estimates are in [0, 1] and CIs are ordered", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(1, 2, 4, 6))
  e <- ar$estimates

  expect_true(all(e$estimate >= 0 & e$estimate <= 1))
  expect_true(all(e$ci_lower >= 0 & e$ci_upper <= 1))
  expect_true(all(e$ci_lower <= e$estimate + 1e-10))
  expect_true(all(e$estimate <= e$ci_upper + 1e-10))
})

test_that("absolute_risk(ipw_cox) is 0 for times before the first event", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(0, 1))

  expect_equal(ar$estimates$estimate[ar$estimates$time == 0], 0)
  expect_equal(ar$estimates$ci_lower[ar$estimates$time == 0], 0)
  expect_equal(ar$estimates$ci_upper[ar$estimates$time == 0], 0)
})

test_that("absolute_risk(ipw_cox) estimate is monotone non-decreasing in time", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(1, 2, 4, 6))
  f <- ar$estimates$estimate
  expect_true(all(diff(f) >= -1e-10))
})

test_that("tidy/print on an ipw_cox absolute-risk object behave", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = c(500 / 100, 10))
  expect_identical(tidy(ar), ar$estimates)
  expect_invisible(print(ar))
})

# -- Exact internal oracle: hand-rolled Breslow == survival::survfit -----------

test_that("IPW Breslow F_x(t) matches survival::survfit on the same coxph", {
  skip_if_not_installed("survival")
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd <- data.frame(x = 1L, z = 0)
  t_eval <- c(1, 2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_mt <- ar$estimates$estimate

  # survfit on the weighted/robust coxph computes the SAME weighted Breslow.
  sf <- survival::survfit(fit$model, newdata = data.frame(x = 1, z = 0))
  f_sf <- 1 - summary(sf, times = t_eval)$surv

  # Both evaluate the identical weighted Breslow step function: agreement is to
  # floating-point precision, not merely sampling tolerance.
  expect_equal(f_mt, f_sf, tolerance = 1e-8)
})

# -- Complicated design: GLM working-model weights + factor confounder --------

# A richer cohort: binary exposure, a continuous confounder, and a 3-level factor
# confounder, so the linear predictor and model matrix exercise factor contrasts.
make_complex_ncc_cohort <- function(n = 6000L, seed = 303L) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.45)
    z <- stats::rnorm(n)
    g <- factor(sample(c("a", "b", "c"), n, TRUE, c(0.5, 0.3, 0.2)))
    bg <- c(a = 0, b = 0.4, c = -0.5)
    lambda0 <- 0.06
    beta_x <- log(2)
    beta_z <- 0.25
    rate <- lambda0 * exp(beta_x * x + beta_z * z + bg[as.character(g)])
    tau <- 9
    tt <- stats::rexp(n, rate)
    cohort <- data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z,
      g = g
    )
    attr(cohort, "truth") <- list(beta_x = beta_x, lambda0 = lambda0, bg = bg)
    cohort
  })
}

test_that("IPW absolute risk matches survfit exactly under GLM weights + factor confounder", {
  skip_if_not_installed("survival")
  cohort <- make_complex_ncc_cohort()
  ncc <- withr::with_seed(
    11L,
    sample_ncc(cohort, "t", "d", m = 4L, incl_prob = TRUE)
  )
  # GLM working-model inclusion weights (the complicated weight path).
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
    confounders = ~ z + g,
    estimator = "ipw_cox"
  )

  # Multiple covariate patterns, each factor level represented.
  nd <- data.frame(
    x = c(1L, 0L, 1L, 1L),
    z = c(0, 0, 1, -1),
    g = factor(c("a", "a", "b", "c"), levels = c("a", "b", "c"))
  )
  t_eval <- c(2, 4, 6, 8)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  # Exact oracle: survfit on the SAME weighted/robust coxph, every covariate
  # pattern. Machine-precision agreement confirms the weighted Breslow and the
  # factor-contrast linear predictor are correct under non-KM weights too.
  dd <- ncc_ipw_analysis_data(fit)
  cx <- survival::coxph(
    survival::Surv(t, d) ~ x + z + g,
    data = dd,
    weights = dd$ipw_weight,
    robust = TRUE
  )
  for (r in seq_len(nrow(nd))) {
    sf <- survival::survfit(cx, newdata = nd[r, , drop = FALSE])
    f_sf <- 1 - summary(sf, times = t_eval)$surv
    f_mt <- ar$estimates$estimate[ar$estimates$row == r]
    expect_equal(f_mt, f_sf, tolerance = 1e-8)
  }
})

test_that("IPW absolute risk recovers truth under a factor confounder (complex DGP)", {
  cohort <- make_complex_ncc_cohort(n = 6000L, seed = 717L)
  truth <- attr(cohort, "truth")
  ncc <- withr::with_seed(
    13L,
    sample_ncc(cohort, "t", "d", m = 4L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~ z + g,
    estimator = "ipw_cox"
  )
  # Analytical F at x = 1, z = 0, g = "a" (the reference level, bg = 0).
  t_eval <- c(2, 4, 6)
  f_true <- 1 - exp(-exp(truth$beta_x * 1) * truth$lambda0 * t_eval)
  nd <- data.frame(x = 1L, z = 0, g = factor("a", levels = c("a", "b", "c")))
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  expect_true(
    all(abs(ar$estimates$estimate - f_true) < 0.06),
    info = paste(
      "Max discrepancy:",
      round(max(abs(ar$estimates$estimate - f_true)), 4)
    )
  )
  expect_true(
    all(ar$estimates$ci_lower <= f_true & f_true <= ar$estimates$ci_upper)
  )
})

# -- Complicated design: tied event times -------------------------------------

test_that("IPW absolute risk matches survfit exactly under tied event times", {
  skip_if_not_installed("survival")
  # Discrete (ceiling) times force heavily tied event times. The coxph fit uses
  # ties = "breslow", so the partial-likelihood coefficients and the hand-rolled
  # Breslow baseline stay mutually consistent; a default Efron fit would leave the
  # plain Breslow baseline inconsistent at ties (off by several points in F).
  cohort <- withr::with_seed(55L, {
    n <- 2500L
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    tt <- stats::rexp(n, 0.15 * exp(log(2) * x + 0.3 * z))
    data.frame(
      id = seq_len(n),
      t = pmin(ceiling(tt), 8),
      d = as.integer(tt <= 8),
      x = x,
      z = z
    )
  })
  expect_true(any(duplicated(cohort$t[cohort$d == 1])))

  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_cox"
  )
  t_eval <- c(2, 4, 6)
  ar <- absolute_risk(fit, newdata = data.frame(x = 1L, z = 0), times = t_eval)
  sf <- survival::survfit(fit$model, newdata = data.frame(x = 1, z = 0))
  f_sf <- 1 - summary(sf, times = t_eval)$surv

  expect_equal(ar$estimates$estimate, f_sf, tolerance = 1e-8)
})

# -- Complicated design: data-dependent confounder basis (poly / splines) -----

test_that("IPW absolute risk matches survfit under a poly() confounder", {
  skip_if_not_installed("survival")
  # poly(z, 2) builds an orthogonal basis that depends on the rows it is computed
  # over. The linear predictor for newdata must reuse the FIT's basis (stored in
  # the terms' predvars), not recompute it from newdata alone -- the latter
  # silently yields a different basis with the same coefficient names and a wrong,
  # often catastrophic, F (off by ~0.4).
  cohort <- make_exp_ncc_cohort(n = 4000L, seed = 88L)
  ncc <- withr::with_seed(
    4L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~ poly(z, 2),
    estimator = "ipw_cox"
  )
  nd <- data.frame(x = c(1L, 1L, 0L), z = c(-1, 1, 0))
  t_eval <- c(2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)

  for (r in seq_len(nrow(nd))) {
    sf <- survival::survfit(fit$model, newdata = nd[r, , drop = FALSE])
    f_sf <- 1 - summary(sf, times = t_eval)$surv
    f_mt <- ar$estimates$estimate[ar$estimates$row == r]
    expect_equal(f_mt, f_sf, tolerance = 1e-8)
  }
})

# -- Edge case: an event exactly at the time origin ---------------------------

test_that("IPW absolute risk handles an event at t = 0 without collapsing the baseline", {
  skip_if_not_installed("survival")
  # Rounded times can place an event exactly at t = 0. The t = 0 fence post must
  # not be prepended on top of a real t = 0 event time, which would duplicate the
  # knot and make approx() collapse it, silently dropping the time-0 increment.
  cohort <- withr::with_seed(71L, {
    n <- 1500L
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    tt <- stats::rexp(n, 0.3 * exp(log(2) * x + 0.3 * z))
    data.frame(
      id = seq_len(n),
      t = pmin(round(tt), 6),
      d = as.integer(tt <= 6),
      x = x,
      z = z
    )
  })
  expect_true(any(cohort$t[cohort$d == 1] == 0))

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
    estimator = "ipw_cox"
  )
  # No "collapsing to unique 'x' values" warning from approx().
  expect_no_warning(
    ar <- absolute_risk(
      fit,
      newdata = data.frame(x = 1L, z = 0),
      times = c(0, 1, 3, 5)
    )
  )
  e <- ar$estimates
  expect_true(all(is.finite(e$estimate)))
  expect_true(all(diff(e$estimate) >= -1e-10))
  # Still exact against survfit on the positive evaluation times.
  sf <- survival::survfit(fit$model, newdata = data.frame(x = 1, z = 0))
  f_sf <- 1 - summary(sf, times = c(1, 3, 5))$surv
  f_mt <- e$estimate[e$time %in% c(1, 3, 5)]
  expect_equal(f_mt, f_sf, tolerance = 1e-8)
})

# -- Oracle: full-cohort survfit ----------------------------------------------

test_that("absolute_risk(ipw_cox) agrees with full-cohort survfit (sampling tol)", {
  skip_if_not_installed("survival")
  cohort <- make_exp_ncc_cohort(n = 4000L, seed = 303L)
  fit <- make_ipw_ar_fit(cohort, m = 3L, seed = 13L)

  nd <- data.frame(x = 1L, z = 0)
  t_eval <- c(2, 4, 6)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_ipw <- ar$estimates$estimate

  # Full-cohort reference: coxph + survfit on the complete cohort.
  cox_full <- survival::coxph(survival::Surv(t, d) ~ x + z, data = cohort)
  sf <- survival::survfit(cox_full, newdata = data.frame(x = 1, z = 0))
  f_full <- 1 - summary(sf, times = t_eval)$surv

  # The NCC subsample (m = 3) reuses controls to estimate the cohort baseline;
  # discrepancies up to 0.05 in F are within expected sampling variability.
  expect_true(
    all(abs(f_ipw - f_full) < 0.05),
    info = paste("Max discrepancy:", round(max(abs(f_ipw - f_full)), 4))
  )
  # The full-cohort F should lie within the (wider) NCC IPW CI.
  expect_true(
    all(
      ar$estimates$ci_lower <= f_full + 0.02 &
        f_full <= ar$estimates$ci_upper + 0.02
    ),
    info = "Full-cohort F_x(t) should lie inside the IPW NCC CI"
  )
})

# -- Truth DGP: exponential survival with known F_x(t) ------------------------

test_that("absolute_risk(ipw_cox) recovers known F_x(t) in exponential DGP", {
  cohort <- make_exp_ncc_cohort(n = 4000L, seed = 404L)
  truth <- attr(cohort, "truth")
  fit <- make_ipw_ar_fit(cohort, m = 3L, seed = 21L)

  # Analytical truth at x = 1, z = 0: F(t) = 1 - exp(-exp(beta_x) * lambda0 * t).
  t_eval <- c(2, 4, 6)
  f_true <- 1 - exp(-exp(truth[["beta_x"]] * 1) * truth[["lambda0"]] * t_eval)

  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_hat <- ar$estimates$estimate

  expect_true(
    all(abs(f_hat - f_true) < 0.06),
    info = paste("Max discrepancy:", round(max(abs(f_hat - f_true)), 4))
  )
  # 95% CI should cover the analytical truth at each evaluation time.
  expect_true(
    all(ar$estimates$ci_lower <= f_true & f_true <= ar$estimates$ci_upper),
    info = "95% CI should cover the analytical truth"
  )
})

# -- Rejection: mismatched newdata --------------------------------------------

test_that("absolute_risk(ipw_cox) rejects mismatched newdata columns", {
  fit <- make_ipw_ar_fit(make_exp_ncc_cohort())
  nd_wrong <- data.frame(x = 1L, foo = 99)
  expect_error(
    absolute_risk(fit, newdata = nd_wrong, times = 4),
    class = "matchatr_bad_input"
  )
})
