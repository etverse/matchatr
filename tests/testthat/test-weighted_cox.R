# Counter-matched nested case-control analysis. sample_ncc_counter_matched()
# draws m controls per case from the *opposite* surrogate stratum and appends
# Langholz-Borgan (1995) log-weights; matcha(design = counter_matched(...),
# estimator = "weighted_cox") fits survival::coxph with those weights as a Cox
# offset and reports the exposure's hazard ratio. Oracles:
#   (1) structural invariants: one case per set, all controls from the opposite
#       surrogate stratum, log-weight formula verified on a controlled micro-cohort
#       where n_same and n_opp are known exactly;
#   (2) full-cohort survival::coxph for the exposure log-HR (the counter-matched
#       subsample targets the same conditional MLE under the Cox model);
#   (3) make_ncc_cohort() DGP with known Cox log-HR for coverage.

# --- structural validity: sample_ncc_counter_matched() -------------------

test_that("a counter-matched sample has one case and m opposite-stratum controls per set", {
  co <- make_ncc_cohort(n = 1200L)
  co$z_bin <- co$x
  m <- 2L
  ncc <- withr::with_seed(11L, sample_ncc_counter_matched(
    co, time = "t", event = "d", surrogate = "z_bin", m = m
  ))

  expect_s3_class(ncc, "data.table")
  # One set per cohort event, exactly one case per set.
  expect_identical(length(unique(ncc$set)), sum(co$d))
  expect_true(all(tapply(ncc$case, ncc$set, sum) == 1L))
  # Set sizes between 2 and m+1 (a smaller set only when the opposite stratum
  # is thinner than m).
  set_sizes <- tabulate(ncc$set)
  expect_true(all(set_sizes >= 2L & set_sizes <= m + 1L))
  # Every control must be from the opposite surrogate stratum to the case.
  by_set <- split(as.data.frame(ncc), ncc$set)
  all_opposite <- all(vapply(by_set, function(s) {
    z_case <- s$z_bin[s$case == 1L]
    all(s$z_bin[s$case == 0L] != z_case)
  }, logical(1)))
  expect_true(all_opposite)
  # The output carries a finite log_w column.
  expect_true("log_w" %in% names(ncc))
  expect_true(all(is.finite(ncc$log_w)))
})

test_that("log_w values match the Langholz-Borgan formula exactly", {
  # Controlled micro-cohort: two events at t=1 and t=3 with fully determined
  # risk-set composition, so expected log-weights can be computed by hand.
  #
  # At t=1 (case: row 1, z_bin=0):
  #   at-risk ex. case: rows 2-7
  #   same stratum (z_bin=0): rows 3,4     => n_same=2
  #   opp  stratum (z_bin=1): rows 2,5,6,7 => n_opp=4, m_take=1
  #   log_w_case = log(2+1) = log(3)
  #   log_w_ctrl = log(4/1) = log(4)
  #
  # At t=3 (case: row 2, z_bin=1):
  #   at-risk ex. case: rows 3,4,5,6,7
  #   same stratum (z_bin=1): rows 5,6,7   => n_same=3
  #   opp  stratum (z_bin=0): rows 3,4     => n_opp=2, m_take=1
  #   log_w_case = log(3+1) = log(4)
  #   log_w_ctrl = log(2/1) = log(2)
  micro <- data.frame(
    id    = 1:7,
    t     = c(1, 3, 5, 5, 5, 5, 5),
    d     = c(1, 1, 0, 0, 0, 0, 0),
    x     = c(0, 1, 0, 0, 1, 1, 1),
    z_bin = c(0, 1, 0, 0, 1, 1, 1)
  )
  ncc <- withr::with_seed(1L, sample_ncc_counter_matched(
    micro, time = "t", event = "d", surrogate = "z_bin", m = 1L
  ))
  s1 <- as.data.frame(ncc[ncc$set == 1L, ])
  expect_equal(s1$log_w[s1$case == 1L], log(3), tolerance = 1e-12)
  expect_equal(s1$log_w[s1$case == 0L], log(4), tolerance = 1e-12)

  s2 <- as.data.frame(ncc[ncc$set == 2L, ])
  expect_equal(s2$log_w[s2$case == 1L], log(4), tolerance = 1e-12)
  expect_equal(s2$log_w[s2$case == 0L], log(2), tolerance = 1e-12)
})

test_that("sample_ncc_counter_matched does not mutate the input cohort", {
  co <- make_ncc_cohort(n = 300L)
  co$z_bin <- co$x
  before <- data.table::copy(data.table::as.data.table(co))
  withr::with_seed(1L, sample_ncc_counter_matched(
    co, time = "t", event = "d", surrogate = "z_bin", m = 1L
  ))
  expect_identical(data.table::as.data.table(co), before)
})

test_that("a control that is a later case may appear in the counter-matched sample", {
  # The defining NCC structure: a subject sampled as a control can itself fail
  # later in the cohort -- it serves as a control before its own event. With a
  # large-enough cohort this must happen (any control with t > risk_time and
  # d == 1 is such a subject).
  co <- make_ncc_cohort(n = 1000L)
  co$z_bin <- as.integer(co$z > 0)  # surrogate independent of x
  ncc <- withr::with_seed(3L, sample_ncc_counter_matched(
    co, time = "t", event = "d", surrogate = "z_bin", m = 1L
  ))
  ctrl <- as.data.frame(ncc[ncc$case == 0L, ])
  expect_true(any(ctrl$d == 1L))
})

test_that("a smaller set is returned when fewer than m opposite-stratum controls exist", {
  # A late failure time where only one opposite-stratum subject is at risk but
  # m=2 is requested: the sampler returns the single available control, not an
  # error. This mirrors the sample_ncc() contract for thin risk sets.
  coh <- data.frame(
    id    = 1:5,
    t     = c(4, 5, 6, 7, 8),
    d     = c(0, 1, 0, 0, 0),
    x     = c(0, 1, 1, 1, 0),
    z_bin = c(0, 1, 1, 1, 0)
  )
  # Case at t=5 (row 2, z_bin=1): at-risk = rows 3,4,5; opp (z_bin=0): row 5 only
  ncc <- withr::with_seed(1L, sample_ncc_counter_matched(
    coh, time = "t", event = "d", surrogate = "z_bin", m = 2L
  ))
  expect_equal(nrow(ncc), 2L)   # case + 1 control (not 2)
  expect_equal(sum(ncc$case), 1L)
})

# --- rejection: sample_ncc_counter_matched() ---------------------------

test_that("a continuous surrogate is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$z_cont <- co$z  # continuous: not binary
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "z_cont"),
    class = "matchatr_bad_input"
  )
})

test_that("a numeric surrogate with more than two distinct values is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$z_tri <- sample(0:2, nrow(co), replace = TRUE)
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "z_tri"),
    class = "matchatr_bad_input"
  )
})

test_that("a three-level factor surrogate is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$z_fac <- factor(sample(c("a", "b", "c"), nrow(co), replace = TRUE))
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "z_fac"),
    class = "matchatr_bad_input"
  )
})

test_that("NA in the surrogate column is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$z_bin <- co$x
  co$z_bin[1L] <- NA_integer_
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "z_bin"),
    class = "matchatr_bad_input"
  )
})

test_that("no eligible opposite-stratum control aborts with matchatr_empty_risk_set", {
  # Subject 1 (z_bin=0) is censored at t=1 and already gone when the case
  # fails at t=3; subjects 3 and 4 (z_bin=1) are at risk at t=3 but share the
  # case's surrogate value. The opposite-stratum pool is therefore empty, so the
  # sampler must abort rather than return a singleton set.
  bad <- data.frame(
    id    = 1:4,
    t     = c(1, 3, 8, 9),
    d     = c(0, 1, 0, 0),
    x     = c(0, 1, 1, 1),
    z_bin = c(0L, 1L, 1L, 1L)
  )
  expect_error(
    sample_ncc_counter_matched(bad, time = "t", event = "d", surrogate = "z_bin"),
    class = "matchatr_empty_risk_set"
  )
})

test_that("a log_w column clash is rejected", {
  co <- make_ncc_cohort(n = 200L)
  co$z_bin <- co$x
  co$log_w <- 0
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "z_bin"),
    class = "matchatr_bad_input"
  )
})

test_that("a missing surrogate column is rejected", {
  co <- make_ncc_cohort(n = 200L)
  expect_error(
    sample_ncc_counter_matched(co, time = "t", event = "d", surrogate = "nope"),
    class = "matchatr_bad_design"
  )
})

# --- weighted Cox: truth-based simulation --------------------------------

test_that("the counter-matched weighted HR recovers the cohort Cox log-HR (truth-based)", {
  # make_ncc_cohort() generates a Cox PH cohort with known beta_x = log(2.2).
  # Surrogate = x (the true exposure) gives maximum counter-matching efficiency.
  # The weighted partial likelihood (Langholz & Borgan 1995) targets the same
  # Cox log-HR as the full cohort under any risk-set design.
  co <- make_ncc_cohort(beta_x = log(2.2), beta_z = log(1.5))
  truth <- attr(co, "truth")
  co$z_bin <- co$x
  ncc <- sample_ncc_counter_matched_fixture(cohort = co, m = 1L)

  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    confounders = ~z,
    estimator = "weighted_cox"
  )
  res <- contrast(fit)

  expect_equal(res$type, "hr")
  # SELF-SCALING band of 3.5 reported SEs (the truth-DGP convention used across
  # all NCC tests): the correct weighted partial likelihood recovers beta_x.
  expect_lt(
    abs(log(res$contrasts$estimate) - unname(truth["beta_x"])),
    3.5 * res$estimates$se
  )
})

# --- weighted Cox: full-cohort coxph oracle ------------------------------

test_that("the counter-matched weighted HR agrees with the full-cohort coxph", {
  # The NCC subsample draws a design-based subset of the cohort's partial
  # likelihood; under the Cox model the full-cohort and subsample estimates
  # both target beta_x, so they must lie within a combined-SE band. The band
  # is conservative because the two estimates are positively correlated
  # (they share the case events).
  co <- make_ncc_cohort(n = 3000L, beta_x = log(2.2), beta_z = log(1.5))
  cox <- survival::coxph(survival::Surv(t, d) ~ x + z, data = co)
  b_cohort <- unname(stats::coef(cox)["x"])
  se_cohort <- sqrt(stats::vcov(cox)["x", "x"])

  co$z_bin <- co$x
  ncc <- sample_ncc_counter_matched_fixture(cohort = co, m = 1L)
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    confounders = ~z,
    estimator = "weighted_cox"
  )
  res <- contrast(fit)
  b_cm <- log(res$contrasts$estimate)
  se_cm <- res$estimates$se

  expect_lt(abs(b_cm - b_cohort), 3.5 * sqrt(se_cm^2 + se_cohort^2))
})

# --- default contrast type and result structure --------------------------

test_that("default contrast type is 'hr' for counter_matched", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  res <- contrast(fit)  # no type = argument
  expect_equal(res$type, "hr")
})

test_that("contrast(type = 'hr') returns finite bounded estimates", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  res <- contrast(fit, type = "hr")
  expect_true(is.finite(res$contrasts$estimate))
  expect_true(is.finite(res$contrasts$ci_lower))
  expect_true(is.finite(res$contrasts$ci_upper))
  expect_gt(res$contrasts$ci_lower, 0)
  expect_gt(res$contrasts$ci_upper, res$contrasts$ci_lower)
})

# --- rejection: contrast() scale and ci_method -------------------------

test_that("contrast(type = 'or') is rejected as unidentified", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  expect_error(contrast(fit, type = "or"), class = "matchatr_unidentified_estimand")
})

test_that("contrast(type = 'difference') is rejected as unidentified", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("ci_method = 'sandwich' is rejected for counter-matched", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  expect_error(
    contrast(fit, type = "hr", ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("ci_method = 'bootstrap' is rejected for counter-matched", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  expect_error(
    contrast(fit, type = "hr", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

# --- rejection: missing weights column in design -------------------------

test_that("a counter-matched design with no weights column aborts at fit time", {
  # counter_matched(weights = NULL) stores no weights slot; fit_weighted_cox()
  # detects the missing column and aborts rather than silently fitting an
  # unweighted model.
  ncc <- sample_ncc_counter_matched_fixture()
  expect_error(
    matcha(
      ncc, "case", "x",
      counter_matched(strata = "set", time = "risk_time"),  # no weights =
      estimator = "weighted_cox"
    ),
    class = "matchatr_bad_design"
  )
})

# --- snapshot tests for key error messages -------------------------------

test_that("the empty-opposite-stratum error reads clearly", {
  bad <- data.frame(
    id = 1:4, t = c(1, 3, 8, 9), d = c(0, 1, 0, 0),
    x = c(0, 1, 1, 1), z_bin = c(0L, 1L, 1L, 1L)
  )
  expect_snapshot(
    sample_ncc_counter_matched(bad, time = "t", event = "d", surrogate = "z_bin"),
    error = TRUE
  )
})

test_that("the missing weights column error reads clearly", {
  ncc <- sample_ncc_counter_matched_fixture()
  expect_snapshot(
    matcha(
      ncc, "case", "x",
      counter_matched(strata = "set", time = "risk_time"),
      estimator = "weighted_cox"
    ),
    error = TRUE
  )
})

test_that("the OR-from-counter-matched error reads clearly", {
  ncc <- sample_ncc_counter_matched_fixture()
  fit <- matcha(
    ncc, "case", "x",
    counter_matched(strata = "set", time = "risk_time", weights = "log_w"),
    estimator = "weighted_cox"
  )
  expect_snapshot(contrast(fit, type = "or"), error = TRUE)
})
