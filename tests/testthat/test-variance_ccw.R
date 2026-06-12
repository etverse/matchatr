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

# When q0 is estimated from a cohort of `prevalence_n` members (rather than
# known), the marginal estimate inherits q-hat0's sampling uncertainty through the
# weights. The analytic (delta-method) interval and the bootstrap (which redraws
# q0* per replicate) implement that variance two different ways, so they are each
# other's oracle; and as the cohort grows the extra term must vanish.

test_that("an estimated q0 widens the interval, matching the bootstrap", {
  skip_if_not_installed("causatr")
  skip_if_not_installed("withr")

  cc <- make_cohort_ccw(n = 3000L, ratio = 3L, seed = 7L)
  q0 <- attr(cc, "q0")
  se_of <- function(r) {
    (r$contrasts$ci_upper - r$contrasts$ci_lower) / (2 * stats::qnorm(0.975))
  }

  for (est in c("ccw_gformula", "ccw_tmle")) {
    # q0 estimated from a small cohort -> a visible extra term.
    fit_est <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0, prevalence_n = 1500L),
      confounders = ~w,
      estimator = est
    )
    expect_false(fit_est$details$prevalence_known)

    se_analytic <- se_of(contrast(
      fit_est,
      "difference",
      ci_method = "sandwich"
    ))
    se_boot <- se_of(withr::with_seed(
      3L,
      contrast(fit_est, "difference", ci_method = "bootstrap", n_boot = 400L)
    ))
    # The analytic delta-method term and the q0*-redraw bootstrap agree.
    expect_equal(se_analytic, se_boot, tolerance = 0.15)

    # As the cohort grows the q-hat0 term vanishes: the estimated-q0 interval
    # collapses onto the known-q0 one.
    fit_known <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0),
      confounders = ~w,
      estimator = est
    )
    fit_huge <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0, prevalence_n = 100000000L),
      confounders = ~w,
      estimator = est
    )
    se_known <- se_of(contrast(fit_known, "difference"))
    se_huge <- se_of(contrast(fit_huge, "difference"))
    expect_equal(se_huge, se_known, tolerance = 0.005)
  }
})

test_that("prevalence_n is validated", {
  expect_error(
    unmatched_cc(prevalence_n = 1000L),
    class = "matchatr_bad_prevalence" # no prevalence to attach to
  )
  expect_error(
    unmatched_cc(prevalence = 0.1, prevalence_n = 1.5),
    class = "matchatr_bad_prevalence" # not a whole number
  )
  expect_error(
    unmatched_cc(prevalence = 0.1, prevalence_n = -10L),
    class = "matchatr_bad_prevalence" # not positive
  )
})
