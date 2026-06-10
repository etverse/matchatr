# Tests for excess_risk.R — the time-varying Aalen cumulative regression
# functions B_j(t) = ∫β_j(s)ds from an IPW nested case-control additive
# (ipw_aalen) fit.
#
# Oracles:
#   - timereg::aalen (without const()) on the same deduplicated, Samuelsen-weighted
#     analysis sample: B̂_j(t) must equal its `cum` and the pointwise SE the square
#     root of its `var.cum` (the Aalen martingale variance) to machine precision.
#   - Truth DGP: an additive-hazards cohort with a known constant excess hazard
#     β_j, so B_j(t) = β_j · t is known in closed form; the IPW estimator recovers
#     it within a SE-scaled band.
#
# No Python / delicatessen oracle: the weighted Aalen cumulative regression is a
# plug-in step function, not an M-estimator delicatessen stacks.

# -- DGP: additive hazard with a known constant excess hazard ------------------

# λ(t | x, z) = λ0 + β_x x + β_z z is constant in t, so T | x,z ~ Exp(rate) and
# the cumulative regression functions are exactly linear: B_x(t) = β_x t,
# B_z(t) = β_z t.
make_excess_cohort <- function(
  n = 4000L,
  lambda0 = 0.05,
  beta_x = 0.04,
  beta_z = 0.01,
  tau = 8,
  seed = 11L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    rate <- pmax(lambda0 + beta_x * x + beta_z * z, 1e-4)
    tt <- stats::rexp(n, rate)
    cohort <- data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z
    )
    attr(cohort, "truth") <- c(beta_x = beta_x, beta_z = beta_z)
    cohort
  })
}

make_excess_fit <- function(cohort, confounders = ~z, m = 3L, seed = 3L) {
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
    estimator = "ipw_aalen"
  )
}

# -- Structural ----------------------------------------------------------------

test_that("excess_risk returns the matchatr_excess_risk structure", {
  fit <- make_excess_fit(make_excess_cohort())
  er <- excess_risk(fit, times = c(2, 4, 6))

  expect_s3_class(er, "matchatr_excess_risk")
  expect_s3_class(er, "matchatr")
  expect_true(data.table::is.data.table(er$estimates))
  expect_named(
    er$estimates,
    c("term", "time", "estimate", "se", "ci_lower", "ci_upper")
  )
  # Covariate terms only — the additive baseline (intercept) is not reported.
  expect_setequal(unique(er$estimates$term), c("x", "z"))
  expect_false("(Intercept)" %in% er$estimates$term)
  expect_identical(er$engine, "ipw_aalen")
  expect_equal(nrow(er$estimates), 2L * 3L)
})

test_that("excess_risk CIs are ordered and B(t) is 0 before the first event", {
  fit <- make_excess_fit(make_excess_cohort())
  er <- excess_risk(fit, times = c(0, 2, 4))
  e <- er$estimates

  expect_true(all(e$ci_lower <= e$estimate + 1e-12))
  expect_true(all(e$estimate <= e$ci_upper + 1e-12))
  z0 <- e[e$time == 0, ]
  expect_true(all(z0$estimate == 0))
  expect_true(all(z0$se == 0))
})

test_that("tidy/print on a matchatr_excess_risk behave", {
  fit <- make_excess_fit(make_excess_cohort())
  er <- excess_risk(fit, times = c(2, 5))
  expect_identical(tidy(er), er$estimates)
  expect_invisible(print(er))
})

# -- Exact oracle: timereg::aalen (no const()) --------------------------------

test_that("excess_risk matches timereg::aalen cum and var.cum to machine precision", {
  skip_if_not_installed("timereg")
  fit <- make_excess_fit(make_excess_cohort())
  t_eval <- c(2, 4, 6)
  er <- excess_risk(fit, times = t_eval)

  # timereg::aalen on the SAME deduplicated, Samuelsen-weighted analysis sample.
  dt <- ncc_ipw_analysis_data(fit)
  ora <- timereg::aalen(
    survival::Surv(t, d) ~ x + z,
    data = dt,
    weights = dt$ipw_weight,
    robust = 1,
    n.sim = 0
  )
  # Step-evaluate timereg's cumulative coefficients / variance at t_eval.
  step_at <- function(M, tt) {
    t(vapply(
      tt,
      function(tk) M[max(which(M[, 1] <= tk)), ],
      numeric(ncol(M))
    ))
  }
  oc <- step_at(ora$cum, t_eval)
  ov <- step_at(ora$var.cum, t_eval)

  # `tm` (not `term`) avoids colliding with the data.table column name in `i`.
  for (tm in c("x", "z")) {
    sub <- er$estimates[er$estimates$term == tm, ]
    expect_equal(sub$estimate, unname(oc[, tm]), tolerance = 1e-8)
    expect_equal(sub$se, sqrt(unname(ov[, tm])), tolerance = 1e-8)
  }
})

test_that("excess_risk matches timereg for a 3-level factor exposure", {
  skip_if_not_installed("timereg")
  # A factor exposure expands to per-level cumulative excess-hazard functions;
  # each must match timereg's corresponding cum column.
  cohort <- withr::with_seed(21L, {
    n <- 5000L
    g <- factor(sample(c("a", "b", "c"), n, TRUE))
    bg <- c(a = 0, b = 0.03, c = -0.02)
    z <- stats::rnorm(n)
    rate <- pmax(0.06 + bg[as.character(g)] + 0.01 * z, 1e-4)
    tt <- stats::rexp(n, rate)
    data.frame(
      id = seq_len(n),
      t = pmin(tt, 8),
      d = as.integer(tt <= 8),
      g = g,
      z = z
    )
  })
  ncc <- withr::with_seed(
    5L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "g",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aalen"
  )
  t_eval <- c(2, 4, 6)
  er <- excess_risk(fit, times = t_eval)
  expect_setequal(unique(er$estimates$term), c("gb", "gc", "z"))

  dt <- ncc_ipw_analysis_data(fit)
  ora <- timereg::aalen(
    survival::Surv(t, d) ~ g + z,
    data = dt,
    weights = dt$ipw_weight,
    robust = 1,
    n.sim = 0
  )
  step_at <- function(M, tt) {
    t(vapply(tt, function(tk) M[max(which(M[, 1] <= tk)), ], numeric(ncol(M))))
  }
  oc <- step_at(ora$cum, t_eval)
  for (tm in c("gb", "gc")) {
    sub <- er$estimates[er$estimates$term == tm, ]
    expect_equal(sub$estimate, unname(oc[, tm]), tolerance = 1e-8)
  }
})

# -- Truth DGP: known constant excess hazard ----------------------------------

test_that("excess_risk recovers the known cumulative excess hazard B_x(t) = beta_x t", {
  cohort <- make_excess_cohort(n = 6000L, seed = 404L)
  truth <- attr(cohort, "truth")
  fit <- make_excess_fit(cohort, seed = 21L)

  t_eval <- c(2, 4, 6)
  er <- excess_risk(fit, times = t_eval)
  sub <- er$estimates[er$estimates$term == "x", ]
  b_true <- truth[["beta_x"]] * t_eval

  # SE-scaled band (the estimator's sampling SD is the reported martingale SE).
  expect_true(
    all(abs(sub$estimate - b_true) < 3.5 * sub$se),
    info = paste(
      "B_x(t) - beta_x t:",
      paste(round(sub$estimate - b_true, 4), collapse = ", ")
    )
  )
  expect_true(all(sub$ci_lower <= b_true & b_true <= sub$ci_upper))
})

# -- Rejections ----------------------------------------------------------------

test_that("excess_risk rejects a non-additive engine and bad times", {
  cohort <- make_excess_cohort(n = 800L, seed = 1L)
  ncc <- withr::with_seed(
    2L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit_cox <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_cox"
  )
  expect_error(
    excess_risk(fit_cox, times = c(2, 4)),
    class = "matchatr_not_implemented"
  )

  fit_add <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aalen"
  )
  expect_error(
    excess_risk(fit_add, times = numeric(0)),
    class = "matchatr_bad_input"
  )
  expect_error(
    excess_risk(fit_add, times = c(1, Inf)),
    class = "matchatr_bad_input"
  )
})
