# Design-weighted marginal causal survival (PHASE_10). The surv_gcomp engine
# g-computes a marginal survival contrast from a sampled cohort: it fits the
# design's own weighted Cox (cch for case-cohort, ipw_cox for nested
# case-control), reads the inverse-probability Breslow baseline, and Horvitz-
# Thompson standardizes the subject-specific absolute risk F(t | x, W) over the
# cohort covariate distribution. contrast() reports the marginal risk difference
# / risk ratio / RMST difference at the requested times with a design-preserving
# bootstrap interval. Oracles:
#   (1) full-cohort g-computation (everyone in the subcohort) is the marginal
#       truth; a case-cohort sample retains every case, so the failure-time grid
#       is identical and the design-weighted estimate must recover it (Monte-Carlo
#       unbiasedness over subcohort draws);
#   (2) the marginalized closed form agrees with absolute_risk() per subject (the
#       standardization is just a weighted average of the validated F_x(t)).
# No Python (delicatessen) oracle: a design-weighted marginal survival contrast
# is not among its M-estimator templates.

# Closed-form marginal point estimate (no bootstrap) for the recovery oracles.
gcomp_point <- function(fit, type, times) {
  std <- surv_gcomp_std_sample(fit)
  fv <- fit
  fv$engine <- std$engine
  surv_gcomp_point(fv, std, fit$exposure, type, times)$contrast
}

# Redraw the subcohort in a fresh cohort copy (cases always retained).
redraw_subcohort <- function(co, frac, seed) {
  out <- co
  out$sub <- 0L
  idx <- withr::with_seed(seed, sample.int(nrow(co), round(frac * nrow(co))))
  out$sub[idx] <- 1L
  out
}

# --- oracle 1: Monte-Carlo unbiasedness vs the full-cohort marginal truth -----

test_that("case-cohort surv_gcomp is unbiased for the full-cohort marginal RD(t)", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 2500L, sub_frac = 0.4, seed = 1L)
  co$all <- 1L
  truth_fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "all", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  eval_t <- as.numeric(stats::quantile(
    co$t[co$d == 1],
    c(0.4, 0.7),
    names = FALSE
  ))
  truth <- gcomp_point(truth_fit, "difference", eval_t)

  reps <- 50L
  rd <- matrix(NA_real_, reps, length(eval_t))
  for (r in seq_len(reps)) {
    cr <- redraw_subcohort(co, 0.4, seed = 100L + r)
    f <- matcha(
      cr,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "sub", time = "t"),
      confounders = ~z,
      estimator = "surv_gcomp"
    )
    rd[r, ] <- gcomp_point(f, "difference", eval_t)
  }
  # The mean over draws recovers the truth tightly (bias << the per-draw SD);
  # treatment is genuinely protective.
  expect_lt(max(abs(colMeans(rd) - truth)), 0.01)
  expect_true(all(truth < 0))
})

test_that("case-cohort surv_gcomp is unbiased for the full-cohort marginal RR(t)", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 2500L, sub_frac = 0.4, seed = 2L)
  co$all <- 1L
  truth_fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "all", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  eval_t <- as.numeric(stats::quantile(
    co$t[co$d == 1],
    c(0.4, 0.7),
    names = FALSE
  ))
  truth <- gcomp_point(truth_fit, "ratio", eval_t)

  reps <- 50L
  rr <- matrix(NA_real_, reps, length(eval_t))
  for (r in seq_len(reps)) {
    cr <- redraw_subcohort(co, 0.4, seed = 200L + r)
    f <- matcha(
      cr,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "sub", time = "t"),
      confounders = ~z,
      estimator = "surv_gcomp"
    )
    rr[r, ] <- gcomp_point(f, "ratio", eval_t)
  }
  expect_lt(max(abs(colMeans(rr) - truth)), 0.02)
  expect_true(all(truth < 1))
})

test_that("case-cohort surv_gcomp is unbiased for the full-cohort marginal RMST diff", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 2500L, sub_frac = 0.4, seed = 3L)
  co$all <- 1L
  truth_fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "all", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  horizon <- as.numeric(stats::quantile(co$t[co$d == 1], 0.7, names = FALSE))
  truth <- gcomp_point(truth_fit, "rmst", horizon)

  reps <- 50L
  rm_ <- numeric(reps)
  for (r in seq_len(reps)) {
    cr <- redraw_subcohort(co, 0.4, seed = 300L + r)
    f <- matcha(
      cr,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "sub", time = "t"),
      confounders = ~z,
      estimator = "surv_gcomp"
    )
    rm_[r] <- gcomp_point(f, "rmst", horizon)
  }
  # Treatment extends restricted mean survival; the mean recovers the truth.
  expect_lt(abs(mean(rm_) - truth), 0.03)
  expect_true(truth > 0)
})

test_that("stratified case-cohort weights recover the full-cohort marginal RD(t)", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 3000L, sub_frac = 0.4, seed = 4L)
  co$all <- 1L
  truth_fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "all", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  eval_t <- as.numeric(stats::quantile(
    co$t[co$d == 1],
    c(0.4, 0.7),
    names = FALSE
  ))
  truth <- gcomp_point(truth_fit, "difference", eval_t)

  reps <- 40L
  rd <- matrix(NA_real_, reps, length(eval_t))
  for (r in seq_len(reps)) {
    cr <- redraw_subcohort(co, 0.4, seed = 400L + r)
    f <- matcha(
      cr,
      outcome = "d",
      exposure = "x",
      design = case_cohort(subcohort = "sub", time = "t", stratum = "stratum"),
      confounders = ~z,
      estimator = "surv_gcomp"
    )
    rd[r, ] <- gcomp_point(f, "difference", eval_t)
  }
  expect_lt(max(abs(colMeans(rd) - truth)), 0.015)
})

# --- oracle 2: the marginalized closed form agrees with absolute_risk() --------

test_that("the marginalized risk equals the absolute_risk() average per subject", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 2000L, sub_frac = 0.5, seed = 5L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  std <- surv_gcomp_std_sample(fit)
  fv <- fit
  fv$engine <- std$engine
  times <- as.numeric(stats::quantile(
    co$t[co$d == 1],
    c(0.3, 0.6),
    names = FALSE
  ))

  mr <- surv_gcomp_marginal_risk(fv, std, "x", times)

  # Independent reconstruction via the public absolute_risk() verb: F^1(t) is the
  # inclusion-weighted mean of F_x(t) over the standardization sample with x = 1.
  nd1 <- std$newdata
  nd1$x <- 1
  ar1 <- absolute_risk(fv, newdata = nd1, times = times)$estimates
  w <- std$weights
  f1_ref <- vapply(
    times,
    function(tt) {
      sub <- ar1[ar1$time == tt, ]
      sum(w * sub$estimate[order(sub$row)]) / sum(w)
    },
    numeric(1)
  )
  expect_equal(mr$f1, f1_ref, tolerance = 1e-8)
})

test_that("the RMST integral equals a dense-grid integral of marginal survival", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, sub_frac = 0.5, seed = 14L)
  fit <- matcha(co, outcome = "d", exposure = "x",
                design = case_cohort(subcohort = "sub", time = "t"),
                confounders = ~ z, estimator = "surv_gcomp")
  std <- surv_gcomp_std_sample(fit)
  fview <- fit
  fview$engine <- std$engine
  horizon <- as.numeric(stats::quantile(co$t[co$d == 1], 0.6, names = FALSE))

  eng <- surv_gcomp_rmst_diff(fview, std, "x", horizon)

  # Independent arbiter: a near-continuous left-Riemann sum of S^a over [0, H],
  # blind to the engine's event-time grid. The marginal survival is a step, so
  # the engine's left-Riemann on the failure-time grid is the exact area; a
  # trapezoidal rule (the prior bug) would disagree by ~1%.
  ug <- seq(0, horizon, length.out = 5000L)
  mr <- surv_gcomp_marginal_risk(fview, std, "x", ug[-1])
  s1 <- c(1, 1 - mr$f1)
  s0 <- c(1, 1 - mr$f0)
  fine1 <- sum(diff(ug) * utils::head(s1, -1L))
  fine0 <- sum(diff(ug) * utils::head(s0, -1L))
  expect_equal(eng[1], fine1, tolerance = 5e-3)
  expect_equal(eng[2], fine0, tolerance = 5e-3)
  expect_equal(eng[3], fine1 - fine0, tolerance = 5e-3)
})

# --- result structure ---------------------------------------------------------

test_that("the result carries per-time contrasts, intervention risks, and tidies", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1800L, seed = 6L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  eval_t <- as.numeric(stats::quantile(
    co$t[co$d == 1],
    c(0.3, 0.6),
    names = FALSE
  ))
  res <- contrast(fit, type = "difference", times = eval_t, n_boot = 40L)

  expect_s3_class(res, "matchatr_result")
  expect_identical(nrow(res$contrasts), 2L)
  expect_identical(res$contrasts$time, eval_t)
  expect_identical(res$contrasts$comparison, rep("treat-all vs treat-none", 2L))
  expect_true(all(is.finite(res$contrasts$ci_lower)))
  expect_true(all(res$contrasts$ci_lower < res$contrasts$estimate))
  expect_true(all(res$contrasts$estimate < res$contrasts$ci_upper))
  # Intervention risks: treat-all and treat-none at each time.
  expect_identical(nrow(res$estimates), 4L)
  expect_identical(
    sort(unique(res$estimates$intervention)),
    c("treat-all", "treat-none")
  )
  expect_identical(res$ci_method, "bootstrap")

  td <- tidy(res)
  expect_true(all(c("time", "type") %in% names(td)))
  expect_identical(td$time, eval_t)
  expect_identical(unique(td$type), "difference")
})

test_that("a two-level factor exposure is accepted", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 7L)
  co$xf <- factor(
    ifelse(co$x == 1L, "treated", "control"),
    levels = c("control", "treated")
  )
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "xf",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  res <- contrast(
    fit,
    type = "difference",
    times = stats::median(co$t[co$d == 1]),
    n_boot = 30L
  )
  expect_true(is.finite(res$contrasts$estimate))
})

test_that("the default contrast scale for surv_gcomp is the risk difference", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 8L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  res <- contrast(fit, times = stats::median(co$t[co$d == 1]), n_boot = 20L)
  expect_identical(res$type, "difference")
})

# --- rejections ---------------------------------------------------------------

test_that("a non-survival scale is rejected as unidentified", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 9L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  tt <- stats::median(co$t[co$d == 1])
  for (bad in c("or", "hr", "af", "excess")) {
    expect_error(
      contrast(fit, type = bad, times = tt),
      class = "matchatr_unidentified_estimand"
    )
  }
})

test_that("times are required and must be positive numerics", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 10L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  expect_error(contrast(fit, type = "difference"), class = "matchatr_bad_input")
  expect_error(
    contrast(fit, type = "difference", times = c(1, -2)),
    class = "matchatr_bad_input"
  )
  expect_error(
    contrast(fit, type = "difference", times = c(1, NA)),
    class = "matchatr_bad_input"
  )
})

test_that("sandwich variance is not available for this engine", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 11L)
  fit <- matcha(
    co,
    outcome = "d",
    exposure = "x",
    design = case_cohort(subcohort = "sub", time = "t"),
    confounders = ~z,
    estimator = "surv_gcomp"
  )
  expect_error(
    contrast(
      fit,
      type = "difference",
      times = stats::median(co$t[co$d == 1]),
      ci_method = "sandwich"
    ),
    class = "matchatr_unsupported_variance"
  )
})

test_that("a non-binary exposure is rejected", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1500L, seed = 12L)
  co$xc <- stats::rnorm(nrow(co))
  expect_error(
    matcha(
      co,
      outcome = "d",
      exposure = "xc",
      design = case_cohort(subcohort = "sub", time = "t"),
      confounders = ~z,
      estimator = "surv_gcomp"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("surv_gcomp is not available for an unmatched or nested design", {
  skip_if_not_installed("survival")
  co <- make_surv_cohort(n = 1000L, seed = 13L)
  # nested_cc has no surv_gcomp row (its marginal survival ships in a later chunk).
  expect_error(
    matcha(
      co,
      outcome = "d",
      exposure = "x",
      design = nested_cc(strata = "stratum", time = "t"),
      confounders = ~z,
      estimator = "surv_gcomp"
    ),
    class = "matchatr_bad_estimator"
  )
  cc <- make_cohort_cc(n = 400L)
  expect_error(
    matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(),
      confounders = ~age,
      estimator = "surv_gcomp"
    ),
    class = "matchatr_bad_estimator"
  )
})
