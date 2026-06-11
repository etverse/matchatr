# Tests for the case-control-weighted g-formula engine (estimator =
# "ccw_gformula"). Two oracles:
#
#   1. PSEUDO-COHORT (exact): matchatr only builds the case-control weights and
#      forwards to causatr, so building the q0-weighted causat() fit + contrast()
#      by hand and comparing must agree to machine precision. The weights are
#      computed from the raw Rose & van der Laan formula here (independently of
#      cc_weights(), which test-weights_cc.R validates separately).
#   2. TRUTH-BASED: a cohort with an analytically known marginal RD / RR / mOR,
#      sampled into an unmatched case-control study; case-control weighting must
#      recover the MARGINAL truth, which differs from the conditional odds ratio a
#      logistic fit reports.

test_that("ccw_gformula forwards exactly to a hand-weighted causatr g-formula", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 4000L, ratio = 4L, seed = 5L)
  q0 <- attr(cc, "q0")

  # Hand-built oracle: the raw case-control weights, then causatr directly. The
  # `quasibinomial` family matches fit_ccw() (the right family for fractional
  # weights; identical mean model to binomial, but silent on non-integer
  # successes).
  y01 <- as.integer(cc$case)
  n1 <- sum(y01 == 1L)
  n0 <- sum(y01 == 0L)
  n <- n1 + n0
  wt <- ifelse(y01 == 1L, q0 / (n1 / n), (1 - q0) / (n0 / n))
  oracle_fit <- causatr::causat(
    cc,
    outcome = "case",
    treatment = "x",
    confounders = ~w,
    estimator = "gcomp",
    family = "quasibinomial",
    weights = as.numeric(wt),
    model_fn = stats::glm
  )

  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = q0),
    confounders = ~w,
    estimator = "ccw_gformula"
  )

  ints <- list(
    treated = causatr::static(1),
    control = causatr::static(0)
  )
  for (ty in c("difference", "ratio", "or")) {
    oracle <- causatr::contrast(
      oracle_fit,
      interventions = ints,
      type = ty,
      reference = "control",
      ci_method = "sandwich"
    )
    res <- contrast(fit, type = ty)
    expect_equal(
      res$contrasts$estimate,
      oracle$contrasts$estimate,
      tolerance = 1e-8
    )
    expect_equal(res$contrasts$se, oracle$contrasts$se, tolerance = 1e-8)
    expect_equal(
      res$estimates$estimate,
      oracle$estimates$estimate,
      tolerance = 1e-8
    )
  }
})

test_that("ccw_gformula recovers the marginal truth, distinct from the conditional OR", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw()
  truth <- attr(cc, "truth")
  q0 <- attr(cc, "q0")

  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = q0),
    confounders = ~w,
    estimator = "ccw_gformula"
  )
  rd <- contrast(fit, type = "difference")$contrasts$estimate
  rr <- contrast(fit, type = "ratio")$contrasts$estimate
  mor <- contrast(fit, type = "or")$contrasts$estimate

  # The case-control-weighted g-formula recovers the analytical marginal truth.
  # Tolerances are relative (testthat semantics): ~2% on the small-magnitude RD,
  # sub-tolerance on the RR / mOR.
  expect_equal(rd, unname(truth["rd"]), tolerance = 0.02)
  expect_equal(rr, unname(truth["rr"]), tolerance = 0.04)
  expect_equal(mor, unname(truth["mor"]), tolerance = 0.05)

  # The conditional odds ratio (a plain logistic fit) recovers exp(beta_x), which
  # is a DIFFERENT number from the marginal OR — the non-collapsibility pin that
  # motivates case-control weighting.
  cond_fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~w,
    estimator = "logistic"
  )
  cond_or <- contrast(cond_fit, type = "or")$contrasts$estimate
  expect_equal(cond_or, unname(truth["cond_or"]), tolerance = 0.05)
  # The marginal OR (truth ~2.30) and the conditional OR (truth 2.50) are
  # genuinely different targets: the two two-sided recoveries land on disjoint
  # tolerance bands (|2.30 - 2.50| = 0.20 > 0.05 + 0.05), so case-control
  # weighting and a plain logistic fit provably estimate different estimands.
  expect_false(
    isTRUE(all.equal(
      unname(truth["mor"]),
      unname(truth["cond_or"]),
      tolerance = 0.05
    ))
  )
})

test_that("ccw_gformula returns a well-formed marginal result", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 4000L, ratio = 4L, seed = 9L)
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_gformula"
  )
  res <- contrast(fit, type = "difference")

  expect_s3_class(res, "matchatr_result")
  expect_identical(res$estimand, "marginal risk difference")
  expect_identical(res$estimator, "ccw_gformula")
  # Two intervention means (treat-all / treat-none) and one contrast row.
  expect_equal(nrow(res$estimates), 2L)
  expect_equal(nrow(res$contrasts), 1L)
  expect_true(all(is.finite(res$contrasts$estimate)))
  expect_true(all(is.finite(res$contrasts$se)))
  # The default contrast scale for a ccw fit is the risk difference.
  expect_identical(contrast(fit)$type, "difference")
})

test_that("the S3 surface works on a ccw fit (no causatr_fit vcov crash)", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 3000L, ratio = 3L, seed = 7L)
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_gformula"
  )

  # The fitted model is a causatr_fit with no coef()/vcov(); the fit-level tidy /
  # summary must surface the marginal contrast instead of calling vcov().
  expect_no_error(print(fit))
  td <- tidy(fit)
  expect_s3_class(td, "data.table")
  # tidy(fit) reports the marginal risk-difference contrast, equal to the value
  # contrast() returns directly.
  rd <- contrast(fit, type = "difference")$contrasts$estimate
  expect_equal(td$estimate, rd)
  expect_identical(td$type, "difference")

  sm <- summary(fit)
  expect_s3_class(sm, "matchatr_result")
  expect_identical(sm$estimand, "marginal risk difference")
})

test_that("ccw_gformula records the variance it actually computed", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 3000L, ratio = 3L, seed = 7L)
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_gformula"
  )
  # A marginal g-formula contrast has only causatr's influence-function /
  # sandwich variance, so both ci_method inputs map to it and the result records
  # "sandwich" (not the requested "model"), with identical SEs.
  m <- contrast(fit, type = "difference", ci_method = "model")
  s <- contrast(fit, type = "difference", ci_method = "sandwich")
  expect_identical(m$ci_method, "sandwich")
  expect_identical(s$ci_method, "sandwich")
  expect_equal(m$contrasts$se, s$contrasts$se)
})

test_that("ccw_gformula rejects a non-binary exposure", {
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  cc$xc <- rnorm(nrow(cc))
  # Confounders are supplied so this isolates the exposure check (the g-formula
  # also requires confounders, tested separately).
  expect_snapshot(
    error = TRUE,
    matcha(
      cc,
      outcome = "case",
      exposure = "xc",
      design = unmatched_cc(prevalence = attr(cc, "q0")),
      confounders = ~w,
      estimator = "ccw_gformula"
    )
  )
})

test_that("ccw_gformula requires confounders to standardize over", {
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  expect_snapshot(
    error = TRUE,
    matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = attr(cc, "q0")),
      estimator = "ccw_gformula"
    )
  )
})

test_that("ccw_gformula rejects a missing prevalence q0", {
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  expect_error(
    matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "ccw_gformula"
    ),
    class = "matchatr_missing_prevalence"
  )
})

test_that("ccw_gformula rejects bootstrap variance and off-scale contrasts", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_gformula"
  )
  expect_snapshot(
    error = TRUE,
    contrast(fit, type = "difference", ci_method = "bootstrap")
  )
  expect_snapshot(error = TRUE, contrast(fit, type = "hr"))
})
