# The nested case-control (NCC) risk-set analysis. matcha(design = nested_cc(),
# estimator = "clogit") fits the conditional partial likelihood with each sampled
# risk set as a stratum -- the same engine as matched CC -- and contrast()
# reports the result as a HAZARD RATIO (type = "hr"), because OR = HR exactly
# under risk-set (incidence-density) sampling (Prentice & Breslow 1978; no
# rare-disease assumption). Oracles: survival::clogit (exact pass-through on the
# NCC data), the full-cohort survival::coxph beta (the design-faithful HR target
# the NCC subsample recovers), and a cohort DGP with a known Cox log-HR sampled
# by sample_ncc_riskset(). The relative-efficiency m/(m+1) at the null is the
# classical NCC information result (Goldstein & Langholz 1992).

# --- exact pass-through oracle on a sampled NCC ---------------------------

test_that("the NCC clogit reproduces survival::clogit exactly", {
  co <- make_ncc_cohort()
  ncc <- sample_ncc_riskset(co, m = 3L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ x + z + strata(set), data = ncc)
  # The wrapper must build the same conditional likelihood as a hand-fit clogit,
  # including the variance (the partial-likelihood information matrix).
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  expect_equal(stats::vcov(fit$model), stats::vcov(oracle))

  res <- contrast(fit)
  hx <- res$contrasts[res$contrasts$comparison == "x", ]
  b <- stats::coef(oracle)["x"]
  se <- sqrt(stats::vcov(oracle)["x", "x"])
  z <- stats::qnorm(0.975)
  expect_equal(hx$estimate, unname(exp(b)), tolerance = 1e-8)
  expect_equal(
    c(hx$ci_lower, hx$ci_upper),
    unname(exp(b + c(-1, 1) * z * se)),
    tolerance = 1e-8
  )
})

# --- truth-based: the CMLE recovers the cohort Cox log-HR -----------------

test_that("the NCC hazard ratio recovers the cohort Cox log-HR (truth-based)", {
  co <- make_ncc_cohort(beta_x = log(2.2))
  truth <- attr(co, "truth")
  ncc <- sample_ncc_riskset(co, m = 3L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  res <- contrast(fit)
  # The risk-set conditional likelihood is unbiased for the cohort Cox beta. Use
  # a SELF-SCALING band of 3.5 reported SEs (the truth-DGP convention): a fixed
  # absolute tolerance below one SD would pass only by luck of the seed.
  expect_lt(
    abs(log(res$contrasts$estimate) - unname(truth["beta_x"])),
    3.5 * res$estimates$se
  )
})

test_that("an adjusted NCC recovers both the exposure and confounder log-HRs", {
  co <- make_ncc_cohort(beta_x = log(2), beta_z = log(1.6))
  truth <- attr(co, "truth")
  ncc <- sample_ncc_riskset(co, m = 4L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  # Pass-through against the hand-fit clogit (same model), then truth recovery.
  oracle <- survival::clogit(case ~ x + z + strata(set), data = ncc)
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  td <- tidy(fit) # log scale, per-term std.error
  for (term in c("x", "z")) {
    est <- td$estimate[td$term == term]
    se <- td$std.error[td$term == term]
    truth_term <- if (term == "x") truth["beta_x"] else truth["beta_z"]
    expect_lt(abs(est - unname(truth_term)), 3.5 * se)
  }
})

# --- OR = HR: the NCC subsample targets the full-cohort Cox beta ----------

test_that("the NCC log-HR agrees with the full-cohort coxph beta (OR = HR)", {
  co <- make_ncc_cohort(beta_x = log(2.2), beta_z = log(1.5))
  # The full-cohort Cox fit is the design-faithful HR the NCC subsample must
  # recover -- this is the OR = HR equivalence, not a rare-disease approximation.
  cox <- survival::coxph(survival::Surv(t, d) ~ x + z, data = co)
  b_cohort <- unname(stats::coef(cox)["x"])
  se_cohort <- sqrt(stats::vcov(cox)["x", "x"])

  ncc <- sample_ncc_riskset(co, m = 4L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  res <- contrast(fit)
  b_ncc <- log(res$contrasts$estimate)
  se_ncc <- res$estimates$se
  # The NCC estimate is a subsample of the cohort, so the two are correlated; a
  # combined-SE band is conservative. A failure here would mean the conditional
  # likelihood is NOT targeting the cohort Cox beta (a broken OR = HR claim).
  expect_lt(abs(b_ncc - b_cohort), 3.5 * sqrt(se_ncc^2 + se_cohort^2))
})

# --- relative efficiency m/(m+1) at the null -----------------------------

test_that("NCC efficiency at the null is m/(m+1) of the full cohort", {
  # At beta = 0 the NCC partial-likelihood information per case is m/(m+1) of the
  # full-cohort information (Goldstein & Langholz 1992), so the estimator
  # variance ratio Var_ncc / Var_cohort is (m+1)/m. This is a Monte-Carlo pin
  # (single seed, large cohort); the band is generous to absorb the sampling of
  # which controls land in each risk set.
  m <- 2L
  co0 <- make_ncc_cohort(n = 8000L, beta_x = 0, beta_z = 0, seed = 99L)
  cox0 <- survival::coxph(survival::Surv(t, d) ~ x, data = co0)
  se_cohort <- sqrt(stats::vcov(cox0)["x", "x"])

  ncc0 <- sample_ncc_riskset(co0, m = m, seed = 5L)
  fit <- matcha(
    ncc0,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  se_ncc <- contrast(fit)$estimates$se
  var_ratio <- (se_ncc / se_cohort)^2
  # (m+1)/m = 1.5 for m = 2; allow a wide band around the asymptotic target.
  expect_gt(var_ratio, 1.30)
  expect_lt(var_ratio, 1.75)
})

# --- factor exposure (per-level HR vs the reference) ---------------------

test_that("a factor exposure reports the per-level hazard ratio and reference", {
  co <- make_ncc_cohort(beta_x = log(2.5))
  co$xf <- factor(
    ifelse(co$x == 1L, "exposed", "unexposed"),
    levels = c("unexposed", "exposed")
  )
  ncc <- sample_ncc_riskset(co, m = 3L)
  fit <- matcha(
    ncc,
    "case",
    "xf",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ xf + strata(set), data = ncc)
  res <- contrast(fit)
  expect_identical(res$type, "hr")
  expect_identical(res$contrasts$comparison, "xfexposed")
  expect_identical(res$reference, "unexposed")
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)["xfexposed"])),
    tolerance = 1e-8
  )
})

# --- effect modification on NCC -> stratum-specific hazard ratios --------

test_that("an effect modifier on NCC reports stratum-specific hazard ratios", {
  co <- make_ncc_cohort(beta_x = log(2))
  co$grp <- factor(ifelse(co$z > 0, "hi", "lo"), levels = c("lo", "hi"))
  ncc <- sample_ncc_riskset(co, m = 4L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    effect_modifier = "grp",
    estimator = "clogit"
  )
  res <- contrast(fit)
  # The new scale label must thread through the stratum-specific assembly: the
  # nested design reports HRs per level, not ORs.
  expect_identical(res$type, "hr")
  expect_identical(res$estimand, "stratum-specific hazard ratio")
  expect_identical(res$reference, "lo")
  # Hand-built linear combos from the interaction clogit: beta_x at the "lo"
  # reference, beta_x + beta_{x:hi} at "hi".
  oracle <- survival::clogit(case ~ x * grp + strata(set), data = ncc)
  b <- stats::coef(oracle)
  expect_equal(
    res$contrasts$estimate,
    unname(exp(c(b["x"], b["x"] + b["x:grphi"]))),
    tolerance = 1e-8
  )
})

# --- result structure / labels ------------------------------------------

test_that("the NCC result is labelled and sized as a hazard ratio", {
  co <- make_ncc_cohort()
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  res <- contrast(fit)
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "hr")
  expect_identical(res$estimand, "hazard ratio")
  expect_identical(res$estimator, "clogit")
  expect_identical(res$engine, "clogit")
  # Analysis n is the rows clogit used (coxph nobs() counts events, not rows).
  expect_identical(res$n, fit$model$n)
  # contrast() with no type defaults to the HR for a nested case-control design.
  expect_identical(contrast(fit)$type, "hr")
})

test_that("tidy reports the NCC hazard ratios matching the oracle", {
  co <- make_ncc_cohort()
  ncc <- sample_ncc_riskset(co, m = 3L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ x + z + strata(set), data = ncc)
  td <- tidy(fit, exponentiate = TRUE)
  # No intercept row (the conditional likelihood has none), one row per term.
  expect_setequal(td$term, c("x", "z"))
  expect_equal(
    td$estimate,
    unname(exp(stats::coef(oracle))[td$term]),
    tolerance = 1e-8
  )
})

test_that("summary labels the NCC table as hazard ratios", {
  co <- make_ncc_cohort(n = 600L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  # The design-aware header reads "hazard ratios", not "odds ratios".
  expect_output(summary(fit), "Conditional hazard ratios")
})

# --- rejections ----------------------------------------------------------

test_that("a non-binary outcome is rejected for the NCC design", {
  co <- make_ncc_cohort(n = 400L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  # A continuous "outcome" is not a case indicator.
  ncc_cont <- ncc
  ncc_cont$case <- ncc_cont$z
  expect_error(
    matcha(
      ncc_cont,
      "case",
      "x",
      nested_cc(strata = "set", time = "risk_time"),
      estimator = "clogit"
    ),
    class = "matchatr_bad_outcome"
  )
  # A 3-level outcome is a polytomous problem, not a risk-set case indicator.
  ncc_multi <- ncc
  ncc_multi$case <- factor(c("a", "b", "c")[(seq_len(nrow(ncc)) %% 3L) + 1L])
  expect_error(
    matcha(
      ncc_multi,
      "case",
      "x",
      nested_cc(strata = "set", time = "risk_time"),
      estimator = "clogit"
    ),
    class = "matchatr_bad_outcome"
  )
})

test_that("requesting an odds ratio from an NCC design is rejected", {
  co <- make_ncc_cohort(n = 600L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  # The risk-set design identifies the hazard ratio, not an odds ratio, even
  # though exp(beta) is numerically the same value.
  expect_error(
    contrast(fit, type = "or"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("requesting a hazard ratio from a matched design is rejected", {
  df <- make_matched_cc(n_sets = 80L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  # A matched case-control design has no risk-set / time structure, so a hazard
  # ratio is a genuinely different (un-targeted) estimand there.
  expect_error(
    contrast(fit, type = "hr"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("RD / RR are rejected as unidentified for the NCC design", {
  co <- make_ncc_cohort(n = 600L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
  expect_error(
    contrast(fit, type = "ratio"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("sandwich / bootstrap CIs are not available for the NCC design", {
  co <- make_ncc_cohort(n = 600L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  expect_error(
    contrast(fit, type = "hr", ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
  expect_error(
    contrast(fit, type = "hr", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("a risk set with no control triggers the uninformative-stratum warning", {
  # Set 2 holds only a case (no eligible control survived to its failure time),
  # so the conditional likelihood drops it. Sets 1 and 3 are informative
  # discordant pairs pointing in opposite directions, so the CMLE on the rest is
  # finite (beta -> 0) and clogit converges cleanly.
  bad <- data.frame(
    case = c(1L, 0L, 1L, 1L, 0L),
    x = c(1L, 0L, 1L, 0L, 1L),
    set = c(1L, 1L, 2L, 3L, 3L),
    t = c(1, 1, 2, 3, 3)
  )
  expect_warning(
    matcha(
      bad,
      "case",
      "x",
      nested_cc(strata = "set", time = "t"),
      estimator = "clogit"
    ),
    class = "matchatr_uninformative_stratum"
  )
})

test_that("an exposure with no within-risk-set variation is not estimable", {
  co <- make_ncc_cohort(n = 300L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  # Make x constant within every set: it has no conditional contribution, so
  # clogit aliases its coefficient to NA and no hazard ratio is identified.
  ncc$x <- as.integer(ncc$set %% 2L)
  fit <- suppressWarnings(matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  ))
  expect_error(
    contrast(fit),
    class = "matchatr_unestimable_exposure"
  )
})

# --- message / print snapshots ------------------------------------------

test_that("NCC contrast rejection and print messages read clearly", {
  co <- make_ncc_cohort(n = 600L)
  ncc <- sample_ncc_riskset(co, m = 2L)
  fit <- matcha(
    ncc,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    estimator = "clogit"
  )
  expect_snapshot(contrast(fit, type = "or"), error = TRUE)
  matched <- matcha(
    make_matched_cc(n_sets = 60L),
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  expect_snapshot(contrast(matched, type = "hr"), error = TRUE)
  expect_snapshot(print(contrast(fit)))
})
