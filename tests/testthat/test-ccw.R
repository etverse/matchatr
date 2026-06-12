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

test_that("the S3 surface works on a ccw fit (no vcov crash)", {
  skip_if_not_installed("causatr")

  cc <- make_cohort_ccw(n = 3000L, ratio = 3L, seed = 7L)
  # ccw_gformula's model is a causatr_fit; ccw_tmle's is a matchatr_ccw_tmle —
  # neither has coef()/vcov(), so the fit-level tidy / summary must surface the
  # marginal contrast instead of calling vcov(). Both route through the same
  # engine-based branch.
  for (est in c("ccw_gformula", "ccw_tmle")) {
    fit <- matcha(
      cc,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(prevalence = attr(cc, "q0")),
      confounders = ~w,
      estimator = est
    )

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
  }
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

test_that("ccw_gformula rejects an off-scale contrast", {
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
  expect_snapshot(error = TRUE, contrast(fit, type = "hr"))
})

test_that("the CCW rejections fire for ipw, aipw, and tmle too", {
  skip_if_not_installed("causatr")

  # The shared ccw_prepare() guard and the contrast_ccw() / contrast_ccw_tmle()
  # guards fire the same classed errors across the CCW family as for ccw_gformula.
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  cc$xc <- rnorm(nrow(cc))
  q0 <- attr(cc, "q0")

  for (est in c("ccw_ipw", "ccw_aipw", "ccw_tmle")) {
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
      contrast(fit, type = "hr"),
      class = "matchatr_unidentified_estimand"
    )
  }
})

test_that("the CCW family complete-cases rows with missing data", {
  skip_if_not_installed("causatr")

  # The CCW family is hand-rolled (TMLE) or delegated (causatr), neither of which
  # tolerates NA the same way, so ccw_prepare() complete-cases the whole family up
  # front. A row with a missing confounder is dropped with a matchatr_dropped_rows
  # warning, and the estimate equals the fit on the same data with that row
  # pre-dropped -- the weights are recomputed on the analysed (complete-case)
  # sample, so the answer is exactly listwise deletion.
  cc <- make_cohort_ccw(n = 3000L, ratio = 3L, seed = 7L)
  q0 <- attr(cc, "q0")
  cc_na <- cc
  cc_na$w[c(3L, 17L, 42L, 88L, 150L)] <- NA
  cc_pre <- cc[!is.na(cc_na$w), ]
  rownames(cc_pre) <- NULL

  # Muffle ONLY the expected matchatr_dropped_rows warning (asserted separately),
  # by class -- not a blanket suppression.
  rd <- function(d, est) {
    withCallingHandlers(
      contrast(
        matcha(
          d,
          outcome = "case",
          exposure = "x",
          design = unmatched_cc(prevalence = q0),
          confounders = ~w,
          estimator = est
        ),
        type = "difference"
      )$contrasts$estimate,
      matchatr_dropped_rows = function(w) invokeRestart("muffleWarning")
    )
  }

  for (est in c("ccw_gformula", "ccw_ipw", "ccw_aipw", "ccw_tmle")) {
    expect_warning(
      matcha(
        cc_na,
        outcome = "case",
        exposure = "x",
        design = unmatched_cc(prevalence = q0),
        confounders = ~w,
        estimator = est
      ),
      class = "matchatr_dropped_rows"
    )
    expect_equal(rd(cc_na, est), rd(cc_pre, est))
  }
})

test_that("matched-CC CCW recovers the marginal truth (matching variable adjusted)", {
  skip_if_not_installed("causatr")

  # Frequency-matched data: the controls are sampled within M-strata, so they are
  # not a representative population control sample. The marginal CCW estimators
  # recover the marginal risk difference when the matching variable M is in the
  # confounders (adjusting for it, not conditioning on the matched sets; Rose &
  # van der Laan 2009). q0 is the cohort prevalence.
  mcc <- make_matched_cohort_ccw()
  truth <- attr(mcc, "truth")
  q0 <- attr(mcc, "q0")

  for (est in c("ccw_gformula", "ccw_aipw", "ccw_tmle")) {
    rd <- contrast(
      matcha(
        mcc,
        outcome = "case",
        exposure = "x",
        design = matched_cc(strata = "set", prevalence = q0),
        confounders = ~M,
        estimator = est
      ),
      type = "difference"
    )$contrasts$estimate
    expect_equal(rd, truth, tolerance = 0.02)
  }
})

test_that("CCW is rejected on a nested case-control (risk-set) design", {
  skip_if_not_installed("causatr")

  # A nested CC is risk-set (incidence-density) sampled, so its controls are not a
  # case-control sample and the q0 binary reweighting does not identify a marginal
  # effect; the inclusion-weighted ipw_cox hazard ratio is the right tool.
  cc <- make_cohort_ccw(n = 2000L, ratio = 3L, seed = 3L)
  cc$set <- rep(seq_len(nrow(cc) / 2L), length.out = nrow(cc))
  for (est in c("ccw_gformula", "ccw_ipw", "ccw_aipw", "ccw_tmle")) {
    expect_error(
      matcha(
        cc,
        outcome = "case",
        exposure = "x",
        design = nested_cc(strata = "set", time = "set", ratio = 1L),
        confounders = ~w,
        estimator = est
      ),
      class = "matchatr_bad_estimator"
    )
  }
})
