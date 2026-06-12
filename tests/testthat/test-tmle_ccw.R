# Tests for the case-control-weighted targeted maximum likelihood estimator
# (estimator = "ccw_tmle"), the one hand-rolled engine in the CCW family.
# Oracles:
#   1. tmle::tmle(obsWeights = ) on the same case-control weights, with a plain
#      glm initial fit (cvQinit = FALSE) and matching gbound — the targeted risk
#      difference must agree essentially exactly, the risk / odds ratios closely.
#   2. TRUTH-BASED recovery of the analytical marginal RD / RR / mOR.
#   3. DOUBLE ROBUSTNESS: like CCW-AIPW, CCW-TMLE recovers the marginal truth when
#      either the outcome or the propensity working model is misspecified.

test_that("ccw_tmle matches tmle::tmle() with case-control obsWeights", {
  skip_if_not_installed("causatr")
  skip_if_not_installed("tmle")

  cc <- make_cohort_ccw(n = 4000L, ratio = 4L, seed = 5L)
  q0 <- attr(cc, "q0")
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = q0),
    confounders = ~w,
    estimator = "ccw_tmle"
  )
  oracle <- tmle_ccw_oracle(cc, q0)

  z <- stats::qnorm(0.975)
  rd <- contrast(fit, type = "difference")$contrasts
  rr <- contrast(fit, type = "ratio")$contrasts
  or <- contrast(fit, type = "or")$contrasts

  # The risk difference is the targeted parameter; matchatr's single-fluctuation
  # TMLE and tmle's agree on it (and its SE) essentially exactly.
  rd_se <- (rd$ci_upper - rd$ci_lower) / (2 * z)
  expect_equal(rd$estimate, oracle$ATE$psi, tolerance = 1e-3)
  expect_equal(rd_se, sqrt(oracle$ATE$var.psi), tolerance = 1e-3)

  # The risk / odds ratios are read off the same targeted Q̄*; they agree with
  # tmle within ~1% (tmle separately fine-tunes the two means), SEs within ~5%.
  rr_se_log <- (log(rr$ci_upper) - log(rr$ci_lower)) / (2 * z)
  or_se_log <- (log(or$ci_upper) - log(or$ci_lower)) / (2 * z)
  expect_equal(rr$estimate, oracle$RR$psi, tolerance = 0.01)
  expect_equal(or$estimate, oracle$OR$psi, tolerance = 0.01)
  expect_equal(rr_se_log, sqrt(oracle$RR$var.log.psi), tolerance = 0.05)
  expect_equal(or_se_log, sqrt(oracle$OR$var.log.psi), tolerance = 0.05)
})

test_that("ccw_tmle recovers the marginal truth", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw()
  truth <- attr(cc, "truth")
  fit <- matcha(
    cc,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(prevalence = attr(cc, "q0")),
    confounders = ~w,
    estimator = "ccw_tmle"
  )
  expect_equal(
    contrast(fit, type = "difference")$contrasts$estimate,
    unname(truth["rd"]),
    tolerance = 0.02
  )
  expect_equal(
    contrast(fit, type = "ratio")$contrasts$estimate,
    unname(truth["rr"]),
    tolerance = 0.05
  )
  expect_equal(
    contrast(fit, type = "or")$contrasts$estimate,
    unname(truth["mor"]),
    tolerance = 0.06
  )
})

test_that("ccw_tmle is doubly robust (consistent if either model is correct)", {
  skip_if_not_installed("causatr")

  rd_of <- function(cc) {
    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = attr(cc, "q0")),
      confounders = ~w,
      estimator = "ccw_tmle"
    )
    contrast(fit, type = "difference")$contrasts$estimate
  }

  # As in the CCW-AIPW double-robustness test: the small-magnitude marginal RD is
  # asserted as an absolute-error band around the analytical truth (all.equal()'s
  # relative tolerance is unreliable at this scale). TMLE recovers the truth
  # whether the outcome (out_wrong) or the propensity (prop_wrong) working model
  # is the misspecified one.
  BAND <- 0.01
  ccA <- make_dr_cohort_ccw("out_wrong")
  expect_lt(abs(rd_of(ccA) - attr(ccA, "truth")), BAND)
  ccB <- make_dr_cohort_ccw("prop_wrong")
  expect_lt(abs(rd_of(ccB) - attr(ccB, "truth")), BAND)
})
