# IPW accelerated failure time for nested case-control (PHASE_7 final chunk).
# matcha(estimator = "ipw_aft") fits survival::survreg(weights = ipw_weight,
# dist = "weibull", robust = TRUE) on the deduplicated Samuelsen-weighted NCC
# sample and contrast(type = "af") reports the time ratio exp(beta) (acceleration
# factor) with the Lin-Wei robust sandwich variance.
#
# Oracles:
#   - Truth DGP: a Weibull AFT cohort with a known time-ratio coefficient; the
#     IPW AFT recovers the full-cohort survreg coefficient within a 3.5-SE band.
#   - Independent reconstruction: KMprob (Samuelsen inclusion probabilities) +
#     survival::survreg on the hand-assembled weighted sample matches matchatr to
#     machine precision, validating the weight assembly and coefficient/variance
#     extraction. No multipleNCC AFT oracle exists (multipleNCC is Cox-only).

# ---- DGP: Weibull AFT with a known time-ratio coefficient -------------------

# log T = 0.5 + beta_x x + beta_z z + sigma * extreme-value error, generated as
# T = exp(eta) * Weibull(shape, 1). survreg(dist = "weibull") recovers beta_x;
# exp(beta_x) is the time ratio (acceleration factor) for the exposure.
make_aft_cohort <- function(
  n = 4000L,
  beta_x = log(0.6),
  beta_z = 0.3,
  shape = 1.3,
  seed = 40L
) {
  withr::with_seed(seed, {
    x <- stats::rbinom(n, 1L, 0.4)
    z <- stats::rnorm(n)
    eta <- 0.5 + beta_x * x + beta_z * z
    tt <- exp(eta) * stats::rweibull(n, shape = shape, scale = 1)
    tau <- stats::quantile(tt, 0.7)
    cohort <- data.frame(
      id = seq_len(n),
      t = pmin(tt, tau),
      d = as.integer(tt <= tau),
      x = x,
      z = z
    )
    attr(cohort, "truth") <- c(beta_x = beta_x)
    cohort
  })
}

aft_estimate <- function(ncc) {
  fit <- matcha(
    ncc,
    outcome = "d",
    exposure = "x",
    design = nested_cc(strata = "set", time = "t"),
    confounders = ~z,
    estimator = "ipw_aft"
  )
  res <- contrast(fit)
  # contrast() reports exp(beta) (time ratio); estimates carries the log-scale
  # coefficient and SE, which is what the survreg oracle is compared on.
  list(
    log_tr = res$estimates$estimate[1],
    log_se = res$estimates$se[1],
    tr = res$contrasts$estimate[1],
    type = res$type
  )
}

# ---- structural -------------------------------------------------------------

test_that("ipw_aft wires to a survreg fit and defaults to the time-ratio scale", {
  cohort <- make_aft_cohort(n = 600L, seed = 1L)
  ncc <- withr::with_seed(
    2L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aft"
  )
  expect_identical(fit$engine, "ipw_aft")
  expect_s3_class(fit$model, "survreg")
  expect_identical(contrast(fit)$type, "af")
})

# ---- truth recovery ---------------------------------------------------------

test_that("ipw_aft recovers the full-cohort AFT time-ratio coefficient", {
  cohort <- make_aft_cohort(n = 4000L, seed = 41L)
  full <- survival::survreg(
    survival::Surv(t, d) ~ x + z,
    data = cohort,
    dist = "weibull"
  )
  target <- unname(stats::coef(full)["x"])

  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  est <- aft_estimate(ncc)

  # The IPW weighted AFT is unbiased for the full-cohort AFT coefficient; 3.5-SE
  # band (the truth-DGP convention).
  expect_lt(abs(est$log_tr - target), 3.5 * est$log_se)
})

# ---- exact independent reconstruction --------------------------------------

test_that("ipw_aft matches an independent KMprob + survreg reconstruction", {
  skip_if_not_installed("multipleNCC")
  cohort <- make_aft_cohort(n = 2500L, seed = 41L)
  ncc <- withr::with_seed(
    3L,
    sample_ncc(cohort, "t", "d", m = 3L, incl_prob = TRUE)
  )
  est <- aft_estimate(ncc)

  n <- nrow(cohort)
  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  ss <- rep(0L, n)
  ss[ctrl_rows] <- 1L
  ss[cohort$d == 1L] <- 2L
  pi_km <- multipleNCC::KMprob(survtime = cohort$t, samplestat = ss, m = 3)
  keep <- sort(unique(c(which(cohort$d == 1L), ctrl_rows)))
  w <- ifelse(cohort$d[keep] == 1L, 1, 1 / pi_km[keep])
  ora <- survival::survreg(
    survival::Surv(t, d) ~ x + z,
    data = cohort[keep, ],
    weights = w,
    dist = "weibull",
    robust = TRUE
  )

  expect_equal(est$log_tr, unname(stats::coef(ora)["x"]), tolerance = 1e-6)
  expect_equal(est$log_se, sqrt(stats::vcov(ora)["x", "x"]), tolerance = 1e-6)
})

test_that("ipw_aft matches a KMprob + survreg reconstruction for a complex covariate set", {
  skip_if_not_installed("multipleNCC")
  # Continuous exposure, a continuous and a three-level-factor confounder, heavy
  # ties at the censoring cap.
  cohort <- make_complex_ncc_cohort(seed = 102L)
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
    estimator = "ipw_aft"
  )
  res <- contrast(fit)

  n <- nrow(cohort)
  ctrl_rows <- unique(ncc$.cohort_row[ncc$case == 0L])
  ss <- rep(0L, n)
  ss[ctrl_rows] <- 1L
  ss[cohort$d == 1L] <- 2L
  pi_km <- multipleNCC::KMprob(survtime = cohort$t, samplestat = ss, m = 3)
  keep <- sort(unique(c(which(cohort$d == 1L), ctrl_rows)))
  w <- ifelse(cohort$d[keep] == 1L, 1, 1 / pi_km[keep])
  ora <- survival::survreg(
    survival::Surv(t, d) ~ xc + z1 + z2,
    data = cohort[keep, ],
    weights = w,
    dist = "weibull",
    robust = TRUE
  )
  pos <- match("xc", names(stats::coef(ora)))

  # The reported exposure coefficient/SE and the full coefficient vector match
  # the independent reconstruction across all confounder columns.
  expect_equal(
    unname(res$estimates$estimate[1]),
    unname(stats::coef(ora)["xc"]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(res$estimates$se[1]),
    sqrt(stats::vcov(ora)[pos, pos]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(stats::coef(fit$model)),
    unname(stats::coef(ora)),
    tolerance = 1e-6
  )
})

# ---- rejections -------------------------------------------------------------

test_that("ipw_aft rejects off-scale contrast types", {
  cohort <- make_aft_cohort(n = 500L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aft"
  )
  for (ty in c("hr", "or", "excess")) {
    expect_error(
      contrast(fit, type = ty),
      class = "matchatr_unidentified_estimand"
    )
  }
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("ipw_aft rejects sandwich and bootstrap variance", {
  cohort <- make_aft_cohort(n = 500L, seed = 1L)
  ncc <- withr::with_seed(
    1L,
    sample_ncc(cohort, "t", "d", m = 2L, incl_prob = TRUE)
  )
  fit <- matcha(
    ncc,
    "d",
    "x",
    nested_cc(strata = "set", time = "t"),
    estimator = "ipw_aft"
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

test_that("ipw_aft requires incl_prob data and a nested design", {
  cohort <- make_aft_cohort(n = 500L, seed = 1L)
  ncc_no <- withr::with_seed(1L, sample_ncc(cohort, "t", "d", m = 2L))
  expect_error(
    matcha(
      ncc_no,
      "d",
      "x",
      nested_cc(strata = "set", time = "t"),
      estimator = "ipw_aft"
    ),
    class = "matchatr_missing_ipw_weights"
  )
  expect_error(
    matcha(cohort, "d", "x", matched_cc(strata = "id"), estimator = "ipw_aft"),
    class = "matchatr_bad_estimator"
  )
})
