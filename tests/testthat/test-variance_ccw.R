# The within-stratum bootstrap (ci_method = "bootstrap") resamples cases and
# controls separately and reports the percentile interval. Its oracle is the
# analytic interval each CCW engine already produces — causatr's
# influence-function sandwich (g-formula / IPW / AIPW) or the targeted EIF
# variance (TMLE): a correct stratified bootstrap recovers that standard error,
# and keeps the analytic point estimate unchanged.

test_that("the CCW within-stratum bootstrap recovers the analytic SE", {
  skip_if_not_installed("causatr")
  skip_if_not_installed("withr")

  cc <- make_cohort_ccw(n = 2500L, ratio = 4L, seed = 11L)
  q0 <- attr(cc, "q0")

  for (est in c("ccw_gformula", "ccw_ipw", "ccw_aipw", "ccw_tmle")) {
    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0),
      confounders = ~w,
      estimator = est
    )
    analytic <- contrast(
      fit,
      type = "difference",
      ci_method = "sandwich"
    )$contrasts
    boot <- withr::with_seed(
      7L,
      contrast(fit, type = "difference", ci_method = "bootstrap", n_boot = 400L)
    )

    # The bootstrap keeps the analytic point estimate (only the interval changes).
    expect_equal(boot$contrasts$estimate, analytic$estimate)
    # The stratified-bootstrap SD recovers the analytic risk-difference SE.
    expect_equal(boot$contrasts$se, analytic$se, tolerance = 0.25)
    expect_identical(boot$ci_method, "bootstrap")
  }
})

test_that("the CCW bootstrap accepts a custom replicate count", {
  skip_if_not_installed("causatr")
  skip_if_not_installed("withr")

  cc <- make_cohort_ccw(n = 1500L, ratio = 3L, seed = 4L)
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_aipw"
  )
  # Two independent bootstraps with the same small n_boot agree closely; the
  # interval is finite and ordered. (n_boot is threaded through contrast()'s `...`.)
  b1 <- withr::with_seed(
    1L,
    contrast(fit, ci_method = "bootstrap", n_boot = 200L)
  )
  b2 <- withr::with_seed(
    2L,
    contrast(fit, ci_method = "bootstrap", n_boot = 200L)
  )
  expect_lt(b1$contrasts$ci_lower, b1$contrasts$ci_upper)
  expect_equal(b1$contrasts$se, b2$contrasts$se, tolerance = 0.3)
})
