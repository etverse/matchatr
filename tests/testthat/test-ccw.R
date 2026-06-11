# Tests for the case-control-weighted causal engines: ccw_gformula (g-computation),
# ccw_ipw (inverse-probability weighting), and ccw_aipw (doubly-robust augmented
# IPW). Oracles:
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
#   3. DOUBLE ROBUSTNESS: ccw_aipw recovers the marginal truth when either (but
#      not both) of the outcome / propensity working models is misspecified.

# matchatr CCW estimator -> causatr causat() estimator, for the hand-built oracle.
.ccw_causat_map <- c(
  ccw_gformula = "gcomp",
  ccw_ipw = "ipw",
  ccw_aipw = "aipw"
)

test_that("the CCW engines forward exactly to a hand-weighted causatr fit", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 4000L, ratio = 4L, seed = 5L)
  q0 <- attr(cc, "q0")

  # Raw Rose & van der Laan case-control weights (independent of cc_weights()).
  y01 <- as.integer(cc$case)
  n1 <- sum(y01 == 1L)
  n0 <- sum(y01 == 0L)
  n <- n1 + n0
  wt <- ifelse(y01 == 1L, q0 / (n1 / n), (1 - q0) / (n0 / n))
  ints <- list(treated = causatr::static(1), control = causatr::static(0))

  for (est in names(.ccw_causat_map)) {
    # Hand-built oracle: causatr directly, with the same quasibinomial family and
    # (for ipw / aipw) the named propensity fitter that fit_ccw() uses.
    oracle_args <- list(
      cc,
      outcome = "case",
      treatment = "x",
      confounders = ~w,
      estimator = .ccw_causat_map[[est]],
      family = "quasibinomial",
      weights = as.numeric(wt),
      model_fn = stats::glm
    )
    if (.ccw_causat_map[[est]] %in% c("ipw", "aipw")) {
      oracle_args$propensity_model_fn <- stats::glm
    }
    oracle_fit <- do.call(causatr::causat, oracle_args)

    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0),
      confounders = ~w,
      estimator = est
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

test_that("ccw_ipw and ccw_aipw recover the marginal truth", {
  skip_if_not_installed("causatr")

  # Both working models are correctly specified here (linear in w), so all three
  # CCW estimators are consistent for the marginal RD / RR / mOR.
  cc <- make_cohort_ccw()
  truth <- attr(cc, "truth")
  q0 <- attr(cc, "q0")

  for (est in c("ccw_ipw", "ccw_aipw")) {
    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = q0),
      confounders = ~w,
      estimator = est
    )
    rd <- contrast(fit, type = "difference")$contrasts$estimate
    rr <- contrast(fit, type = "ratio")$contrasts$estimate
    mor <- contrast(fit, type = "or")$contrasts$estimate
    expect_equal(rd, unname(truth["rd"]), tolerance = 0.02)
    expect_equal(rr, unname(truth["rr"]), tolerance = 0.05)
    expect_equal(mor, unname(truth["mor"]), tolerance = 0.06)
  }
})

test_that("ccw_aipw is doubly robust (consistent if either model is correct)", {
  skip_if_not_installed("causatr")

  rd_of <- function(cc, est) {
    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = attr(cc, "q0")),
      confounders = ~w,
      estimator = est
    )
    contrast(fit, type = "difference")$contrasts$estimate
  }

  # The marginal risk difference here is ~0.04-0.08, below the magnitude at which
  # all.equal()'s relative tolerance applies, so recovery is asserted as an
  # absolute-error band |estimate - truth| < BAND around the analytical truth (a
  # two-sided check, the same idiom causatr's DR tests use). The DGP gives a clean
  # 10x gap: the consistent estimators land within ~0.003, the misspecified one
  # outside ~0.03, so BAND = 0.01 separates them robustly.
  BAND <- 0.01

  # Scenario A: the OUTCOME model is misspecified (`~ w` misses a w^2 term) and the
  # propensity is correct. CCW-AIPW (and CCW-IPW, which uses the correct
  # propensity) recover the marginal truth; CCW-g-formula (wrong outcome) does not.
  ccA <- make_dr_cohort_ccw("out_wrong")
  truthA <- attr(ccA, "truth")
  expect_lt(abs(rd_of(ccA, "ccw_aipw") - truthA), BAND)
  expect_lt(abs(rd_of(ccA, "ccw_ipw") - truthA), BAND)
  expect_gt(abs(rd_of(ccA, "ccw_gformula") - truthA), BAND)

  # Scenario B: the PROPENSITY model is misspecified and the outcome is correct.
  # CCW-AIPW (and CCW-g-formula, correct outcome) recover the truth; CCW-IPW
  # (wrong propensity) does not.
  ccB <- make_dr_cohort_ccw("prop_wrong")
  truthB <- attr(ccB, "truth")
  expect_lt(abs(rd_of(ccB, "ccw_aipw") - truthB), BAND)
  expect_lt(abs(rd_of(ccB, "ccw_gformula") - truthB), BAND)
  expect_gt(abs(rd_of(ccB, "ccw_ipw") - truthB), BAND)
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

test_that("the CCW rejections fire for ipw and aipw too", {
  skip_if_not_installed("causatr")

  # The fit_ccw() / contrast_ccw() guards are shared across the CCW family, so the
  # same classed errors fire for ccw_ipw and ccw_aipw as for ccw_gformula.
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  cc$xc <- rnorm(nrow(cc))
  q0 <- attr(cc, "q0")

  for (est in c("ccw_ipw", "ccw_aipw")) {
    expect_error(
      matcha(cc, "case", "x", unmatched_cc(prevalence = q0), estimator = est),
      class = "matchatr_bad_input" # no confounders
    )
    expect_error(
      matcha(
        cc,
        "case",
        "xc",
        unmatched_cc(prevalence = q0),
        confounders = ~w,
        estimator = est
      ),
      class = "matchatr_bad_input" # non-binary exposure
    )
    expect_error(
      matcha(cc, "case", "x", unmatched_cc(), estimator = est),
      class = "matchatr_missing_prevalence"
    )
    fit <- matcha(
      cc,
      "case",
      "x",
      unmatched_cc(prevalence = q0),
      confounders = ~w,
      estimator = est
    )
    expect_error(
      contrast(fit, ci_method = "bootstrap"),
      class = "matchatr_unsupported_variance"
    )
    expect_error(
      contrast(fit, type = "hr"),
      class = "matchatr_unidentified_estimand"
    )
  }
})
