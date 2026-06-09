# Tests for absolute_risk.R — IPW Breslow cumulative baseline hazard and
# F_x(t) from case-cohort fits.
#
# Oracle: full-cohort survival::survfit(coxph) on nwtco.
# Truth DGP: exponential survival with known beta; analytical F_x(t) = 1 - exp(-lambda * t).
#
# No Python oracle: the IPW Breslow estimator is a plug-in step function,
# not an M-estimator, and is not covered by delicatessen.

# -- Shared test helpers -------------------------------------------------------

make_nwtco_fit <- function(method = "Prentice") {
  nwtco2 <- survival::nwtco
  nwtco2$subcohort <- as.logical(nwtco2$in.subcohort)
  nwtco2$stage <- factor(nwtco2$stage)
  matcha(
    nwtco2,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "subcohort",
      time = "edrel",
      method = method
    ),
    confounders = ~ stage + age,
    estimator = "cch"
  )
}

# Exponential survival DGP with a known beta; subcohort is 30% of cohort.
# beta_x = 0.5, beta_z = 0.3, baseline hazard lambda0 = 0.05.
# Analytical truth: F_x(t) = 1 - exp(-exp(beta_x * x + beta_z * z) * lambda0 * t).
make_exp_truth_data <- function(n = 2000, seed = 42) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.5)
  z <- rnorm(n)
  beta_x <- 0.5
  beta_z <- 0.3
  lambda0 <- 0.05
  true_t <- rexp(n, rate = lambda0 * exp(beta_x * x + beta_z * z))
  cens_t <- rexp(n, rate = 0.02)
  obs_t <- pmin(true_t, cens_t)
  event <- as.integer(true_t <= cens_t)
  sc <- rep(FALSE, n)
  sc[sort(sample(n, round(n * 0.3)))] <- TRUE
  list(
    data = data.frame(obs_t = obs_t, event = event, x = x, z = z, sc = sc),
    beta_x = beta_x,
    beta_z = beta_z,
    lambda0 = lambda0
  )
}

# -- Structural checks ---------------------------------------------------------

test_that("absolute_risk returns correct S3 class and list structure", {
  fit <- make_nwtco_fit()
  nd <- data.frame(histol = 1L, stage = factor(2, levels = 1:4), age = 3)
  ar <- absolute_risk(fit, newdata = nd, times = c(500, 1000))

  expect_s3_class(ar, "matchatr_absolute_risk")
  expect_s3_class(ar, "matchatr")
  expect_true(data.table::is.data.table(ar$estimates))
  expect_named(
    ar$estimates,
    c("row", "time", "estimate", "ci_lower", "ci_upper")
  )
})

test_that("absolute_risk returns one row per (newdata row, time) combination", {
  fit <- make_nwtco_fit()
  nd <- data.frame(
    histol = c(0L, 1L),
    stage = factor(c(1, 3), levels = 1:4),
    age = c(2, 5)
  )
  times <- c(200, 800, 2000)
  ar <- absolute_risk(fit, newdata = nd, times = times)

  expect_equal(nrow(ar$estimates), nrow(nd) * length(times))
  # row column indexes back to newdata rows
  expect_equal(sort(unique(ar$estimates$row)), 1:2)
})

test_that("tidy.matchatr_absolute_risk returns the estimates data.table", {
  fit <- make_nwtco_fit()
  nd <- data.frame(histol = 0L, stage = factor(1, levels = 1:4), age = 2)
  ar <- absolute_risk(fit, newdata = nd, times = c(500, 1000))

  tidy_out <- tidy(ar)
  expect_s3_class(tidy_out, "data.table")
  expect_identical(tidy_out, ar$estimates)
})

test_that("print.matchatr_absolute_risk runs without error", {
  fit <- make_nwtco_fit()
  nd <- data.frame(histol = 1L, stage = factor(2, levels = 1:4), age = 3)
  ar <- absolute_risk(fit, newdata = nd, times = c(500, 1000))

  expect_invisible(print(ar))
})

# -- Monotonicity and bounds checks -------------------------------------------

test_that("absolute_risk estimates are in [0, 1] and CIs are ordered", {
  fit <- make_nwtco_fit()
  nd <- data.frame(histol = 1L, stage = factor(2, levels = 1:4), age = 3)
  ar <- absolute_risk(fit, newdata = nd, times = c(100, 500, 1000, 2000))
  e <- ar$estimates

  expect_true(all(e$estimate >= 0 & e$estimate <= 1))
  expect_true(all(e$ci_lower >= 0 & e$ci_upper <= 1))
  expect_true(all(e$ci_lower <= e$estimate + 1e-10))
  expect_true(all(e$estimate <= e$ci_upper + 1e-10))
})

test_that("absolute_risk is 0 for times before the first event", {
  fit <- make_nwtco_fit()
  nd <- data.frame(histol = 1L, stage = factor(2, levels = 1:4), age = 3)
  # nwtco event times are > 0; t = 0 is always before the first event
  ar <- absolute_risk(fit, newdata = nd, times = c(0, 1))

  expect_equal(ar$estimates$estimate[ar$estimates$time == 0], 0)
  expect_equal(ar$estimates$ci_lower[ar$estimates$time == 0], 0)
  expect_equal(ar$estimates$ci_upper[ar$estimates$time == 0], 0)
})

# -- Oracle comparison: full-cohort survfit ------------------------------------

test_that("absolute_risk agrees with full-cohort survfit within sampling tolerance", {
  skip_if_not_installed("survival")

  fit_cch <- make_nwtco_fit(method = "Prentice")
  nd <- data.frame(histol = 1L, stage = factor(2, levels = 1:4), age = 3)
  t_eval <- c(500, 1000, 2000)
  ar <- absolute_risk(fit_cch, newdata = nd, times = t_eval)

  # Full-cohort reference: coxph + survfit on the complete nwtco
  nwtco2 <- survival::nwtco
  nwtco2$stage <- factor(nwtco2$stage)
  cox_full <- survival::coxph(
    survival::Surv(edrel, rel) ~ histol + stage + age,
    data = nwtco2
  )
  sf <- survival::survfit(cox_full, newdata = nd)
  f_full <- 1 - summary(sf, times = t_eval)$surv

  # The subcohort is ~16.6% of the cohort; discrepancies up to 0.06 in F are
  # within expected sampling variability. Assert using a loose tolerance —
  # the goal is to confirm the estimator is in the right order of magnitude
  # and that the CI includes the full-cohort reference.
  f_ipw <- ar$estimates$estimate

  expect_true(
    all(abs(f_ipw - f_full) < 0.06),
    info = paste("Max discrepancy:", round(max(abs(f_ipw - f_full)), 4))
  )
  # The full-cohort estimate should lie within the (wider) case-cohort CI
  expect_true(
    all(
      ar$estimates$ci_lower <= f_full + 0.02 &
        f_full <= ar$estimates$ci_upper + 0.02
    ),
    info = "Full-cohort F_x(t) should lie inside the case-cohort CI"
  )
})

# -- Truth DGP: exponential survival with known F_x(t) ------------------------

test_that("absolute_risk recovers known F_x(t) in exponential DGP", {
  skip_if_not_installed("survival")

  dat <- make_exp_truth_data(n = 2000, seed = 42)
  fit <- matcha(
    dat$data,
    outcome = "event",
    exposure = "x",
    design = case_cohort(
      subcohort = "sc",
      time = "obs_t",
      method = "Prentice"
    ),
    confounders = ~z,
    estimator = "cch"
  )

  # Analytical truth at x = 1, z = 0
  t_eval <- c(5, 10, 20)
  truth <- 1 - exp(-exp(dat$beta_x * 1 + dat$beta_z * 0) * dat$lambda0 * t_eval)

  nd <- data.frame(x = 1L, z = 0)
  ar <- absolute_risk(fit, newdata = nd, times = t_eval)
  f_hat <- ar$estimates$estimate

  # With n=2000 and 30% subcohort, expect within ~0.06 of truth
  expect_true(
    all(abs(f_hat - truth) < 0.06),
    info = paste("Max discrepancy:", round(max(abs(f_hat - truth)), 4))
  )

  # 95% CI should cover the analytical truth at each evaluation time
  expect_true(
    all(ar$estimates$ci_lower <= truth & truth <= ar$estimates$ci_upper),
    info = "95% CI should cover the analytical truth"
  )
})

# -- Borgan II (stratified subcohort) -----------------------------------------

test_that("absolute_risk works for Borgan II stratified subcohort", {
  nwtco2 <- survival::nwtco
  nwtco2$subcohort <- as.logical(nwtco2$in.subcohort)
  nwtco2$stage <- factor(nwtco2$stage)

  fit_b2 <- matcha(
    nwtco2,
    outcome = "rel",
    exposure = "histol",
    design = case_cohort(
      subcohort = "subcohort",
      time = "edrel",
      method = "II.Borgan",
      stratum = "stage"
    ),
    confounders = ~age,
    estimator = "cch"
  )
  nd <- data.frame(histol = 1L, age = 3)
  ar <- absolute_risk(fit_b2, newdata = nd, times = c(500, 1000))

  expect_equal(ar$method, "II.Borgan")
  expect_true(all(ar$estimates$estimate > 0 & ar$estimates$estimate < 1))
  expect_true(all(ar$estimates$ci_lower <= ar$estimates$estimate + 1e-10))
})

# -- Rejection paths -----------------------------------------------------------

test_that("absolute_risk rejects non-cch engine with matchatr_not_implemented", {
  fit_fake <- structure(
    list(engine = "clogit", model = list(1)),
    class = "matchatr_fit"
  )
  expect_error(
    absolute_risk(fit_fake, newdata = data.frame(x = 1), times = 1),
    class = "matchatr_not_implemented"
  )
})

test_that("absolute_risk rejects unestimated fit with matchatr_not_estimated", {
  fit_null <- structure(
    list(engine = "cch", model = NULL),
    class = "matchatr_fit"
  )
  expect_error(
    absolute_risk(fit_null, newdata = data.frame(x = 1), times = 1),
    class = "matchatr_not_estimated"
  )
})

test_that("absolute_risk rejects empty newdata with matchatr_bad_input", {
  fit_fake <- structure(
    list(engine = "cch", model = list(1)),
    class = "matchatr_fit"
  )
  expect_error(
    absolute_risk(fit_fake, newdata = data.frame(), times = 1),
    class = "matchatr_bad_input"
  )
})

test_that("absolute_risk rejects non-finite times with matchatr_bad_input", {
  fit_fake <- structure(
    list(engine = "cch", model = list(1)),
    class = "matchatr_fit"
  )
  expect_error(
    absolute_risk(fit_fake, newdata = data.frame(x = 1), times = c(1, Inf)),
    class = "matchatr_bad_input"
  )
  expect_error(
    absolute_risk(fit_fake, newdata = data.frame(x = 1), times = c(1, NA)),
    class = "matchatr_bad_input"
  )
})

test_that("absolute_risk rejects mismatched newdata columns with matchatr_bad_input", {
  fit_cch <- make_nwtco_fit()
  # Missing 'stage' and 'age' columns
  nd_wrong <- data.frame(histol = 1L, foo = 99)
  expect_error(
    absolute_risk(fit_cch, newdata = nd_wrong, times = 500),
    class = "matchatr_bad_input"
  )
})

test_that("absolute_risk snapshots: rejection error messages", {
  fit_fake <- structure(
    list(engine = "clogit", model = list(1)),
    class = "matchatr_fit"
  )
  expect_snapshot(
    absolute_risk(fit_fake, newdata = data.frame(x = 1), times = 1),
    error = TRUE
  )
})
