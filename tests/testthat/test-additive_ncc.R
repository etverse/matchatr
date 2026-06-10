# IPW additive-hazards (Lin-Ying) for nested case-control (PHASE_7 final chunk).
# matcha(estimator = "ipw_aalen") fits the weighted constant additive-hazards
# estimator (lin_ying_additive()) on the deduplicated Samuelsen-weighted NCC
# sample and contrast(type = "excess") reports the excess hazard gamma (additive
# rate difference) with the robust sandwich variance.
#
# Oracles:
#   - Truth DGP: a constant-additive-hazard cohort (lambda = a0 + ax*x + az*z,
#     z binary, so the additive structure is exact and lambda stays positive);
#     the IPW additive estimator recovers the known excess hazard ax within a
#     3.5-SE band, and a known per-level structure for a factor exposure.
#   - timereg::aalen on the same weighted analysis sample is the external oracle
#     for the estimator: the excess-hazard point estimate agrees to machine
#     precision, and the robust SE to finite-sample order (the two robust
#     variances are asymptotically equivalent). No multipleNCC additive oracle
#     exists (Cox-only).

# ---- DGP: constant additive hazard with a known excess hazard ---------------

# lambda(t | x, z) = a0 + ax * x + az * z, with z binary so lambda > 0 and the
# additive (const) model is exactly specified. Constant hazard -> exponential
# event times. The exposure's excess hazard is ax.
make_additive_cohort <- function(
  n = 4000L,
  a0 = 0.10,
  ax = 0.08,
  az = 0.05,
  tau = 5,
  seed = 50L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rbinom(n, 1L, 0.5)
    lambda <- a0 + ax * x + az * z
    tt <- stats::rexp(n, lambda)
    cohort <- data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z
    )
    attr(cohort, "truth") <- c(ax = ax, az = az)
    cohort
  })
}

aalen_excess <- function(ncc, exposure = "x") {
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = exposure,
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aalen"
  )
  res <- contrast(fit)
  list(
    gamma = res$contrasts$estimate,
    se = res$contrasts$se,
    type = res$type,
    reference = res$reference,
    terms = res$contrasts$comparison
  )
}

# ---- structural -------------------------------------------------------------

test_that("ipw_aalen wires to a matchatr_aalen fit and defaults to the excess scale", {
  cohort <- make_additive_cohort(n = 800L, seed = 1L)
  ncc <- withr::with_seed(
    2L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aalen"
  )
  expect_identical(fit$engine, "ipw_aalen")
  expect_s3_class(fit$model, "matchatr_aalen")
  res <- contrast(fit)
  expect_identical(res$type, "excess")
  # The excess hazard is reported on the bare exposure term.
  expect_identical(res$contrasts$comparison, "x")
})

# ---- truth recovery ---------------------------------------------------------

test_that("ipw_aalen recovers the known additive excess hazard", {
  cohort <- make_additive_cohort(n = 4000L, ax = 0.08, seed = 51L)
  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  est <- aalen_excess(ncc)

  # The IPW weighted additive estimator is consistent for the excess hazard ax;
  # 3.5-SE band (the truth-DGP convention).
  expect_lt(abs(est$gamma[1] - 0.08), 3.5 * est$se[1])
})

test_that("ipw_aalen recovers per-level excess hazards for a factor exposure", {
  # A three-level exposure with a known per-level additive structure.
  cohort <- withr::with_seed(60L, {
    n <- 4000L
    g <- factor(
      sample(c("lo", "mid", "hi"), n, TRUE),
      levels = c("lo", "mid", "hi")
    )
    z <- stats::rbinom(n, 1L, 0.5)
    code <- as.integer(g) - 1L # lo=0, mid=1, hi=2
    lambda <- 0.08 + 0.04 * code + 0.05 * z # excess: mid-lo = 0.04, hi-lo = 0.08
    tt <- stats::rexp(n, lambda)
    data.frame(
      id = seq_len(n),
      t = pmin(tt, 5),
      d = as.integer(tt <= 5),
      x = g,
      z = z
    )
  })
  ncc <- withr::with_seed(
    7L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aalen"
  )
  est <- contrast(fit)

  expect_identical(est$contrasts$comparison, c("xmid", "xhi"))
  expect_identical(est$reference, "lo")
  expect_lt(abs(est$contrasts$estimate[1] - 0.04), 3.5 * est$contrasts$se[1])
  expect_lt(abs(est$contrasts$estimate[2] - 0.08), 3.5 * est$contrasts$se[2])

  # Both per-level excess hazards match timereg::aalen exactly.
  skip_if_not_installed("timereg")
  const <- timereg::const
  dt <- ncc[!duplicated(ncc$.cohort_row), ]
  dt$ipw_weight[dt$d == 1L] <- 1
  ora <- timereg::aalen(
    survival::Surv(t, d) ~ const(x) + const(z),
    data = as.data.frame(dt),
    weights = dt$ipw_weight,
    n.sim = 0
  )
  ora_x <- ora$gamma[grepl("const(x)", rownames(ora$gamma), fixed = TRUE), 1L]
  expect_equal(unname(est$contrasts$estimate), unname(ora_x), tolerance = 1e-6)
})

# ---- external oracle: timereg::aalen on the same weighted sample ------------

test_that("ipw_aalen matches timereg::aalen on the same weighted analysis sample", {
  skip_if_not_installed("timereg")
  # Make the const() formula special resolvable without attaching timereg.
  const <- timereg::const
  cohort <- make_additive_cohort(n = 2500L, seed = 51L)
  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  est <- aalen_excess(ncc)

  # The deduplicated, case-weight-forced sample the engine fits on.
  dt <- ncc[!duplicated(ncc$.cohort_row), ]
  dt$ipw_weight[dt$d == 1L] <- 1
  ora <- timereg::aalen(
    survival::Surv(t, d) ~ const(x) + const(z),
    data = as.data.frame(dt),
    weights = dt$ipw_weight,
    n.sim = 0
  )
  ora_g <- unname(ora$gamma[rownames(ora$gamma) == "const(x)", 1L])
  ora_se <- sqrt(ora$robvar.gamma[
    rownames(ora$gamma) == "const(x)",
    rownames(ora$gamma) == "const(x)"
  ])

  # The point estimate solves the same estimating equation, so it agrees to
  # machine precision; the two robust sandwich variances are asymptotically
  # equivalent and agree to finite-sample order (within 5%).
  expect_equal(unname(est$gamma[1]), unname(ora_g), tolerance = 1e-6)
  expect_equal(unname(est$se[1]), unname(ora_se), tolerance = 0.05)
})

test_that("ipw_aalen matches timereg::aalen for a complex covariate set", {
  skip_if_not_installed("timereg")
  const <- timereg::const
  # Continuous exposure, a continuous and a three-level-factor confounder, heavy
  # ties at the censoring cap: a setting with no closed-form truth where the
  # external oracle still gives the exact fit.
  cohort <- make_complex_ncc_cohort(seed = 101L)
  ncc <- withr::with_seed(
    9L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "xc",
    nested_cc(strata = "set", time = "t"),
    confounders = ~ z1 + z2,
    estimator = "ipw_aalen"
  )
  res <- contrast(fit)

  dt <- ncc[!duplicated(ncc$.cohort_row), ]
  dt$ipw_weight[dt$d == 1L] <- 1
  ora <- timereg::aalen(
    survival::Surv(t, d) ~ const(xc) + const(z1) + const(z2),
    data = as.data.frame(dt),
    weights = dt$ipw_weight,
    n.sim = 0
  )
  ora_xc <- unname(ora$gamma[rownames(ora$gamma) == "const(xc)", 1L])
  ora_xc_se <- sqrt(ora$robvar.gamma[
    rownames(ora$gamma) == "const(xc)",
    rownames(ora$gamma) == "const(xc)"
  ])

  # The reported (continuous) exposure excess hazard matches timereg exactly, and
  # so does the full coefficient vector (exposure + every confounder column).
  expect_equal(unname(res$contrasts$estimate[1]), ora_xc, tolerance = 1e-6)
  expect_equal(unname(res$contrasts$se[1]), unname(ora_xc_se), tolerance = 0.05)
  expect_equal(
    sort(unname(fit$model$gamma)),
    sort(unname(ora$gamma[, 1L])),
    tolerance = 1e-6
  )
})

# ---- rejections -------------------------------------------------------------

test_that("ipw_aalen rejects off-scale contrast types", {
  cohort <- make_additive_cohort(n = 500L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aalen"
  )
  for (ty in c("hr", "or", "af")) {
    expect_error(
      contrast(fit, type = ty),
      class = "matchatr_unidentified_estimand"
    )
  }
  expect_error(
    contrast(fit, type = "ratio"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("ipw_aalen rejects sandwich and bootstrap variance", {
  cohort <- make_additive_cohort(n = 500L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aalen"
  )
  expect_error(
    contrast(fit, ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
  expect_error(
    contrast(fit, ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("ipw_aalen requires incl_prob data and a nested design", {
  cohort <- make_additive_cohort(n = 500L, seed = 1L)
  ncc_no <- withr::with_seed(1L, sample_ncc(cohort, "t", "d", m = 2L))
  expect_error(
    matcha(
      ncc_no,
      "d",
      "x",
      nested_cc(strata = "set", time = "t"),
      estimator = "ipw_aalen"
    ),
    class = "matchatr_missing_ipw_weights"
  )
  expect_error(
    matcha(
      cohort,
      "d",
      "x",
      matched_cc(strata = "id"),
      estimator = "ipw_aalen"
    ),
    class = "matchatr_bad_estimator"
  )
})
