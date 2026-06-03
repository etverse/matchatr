# The 1:1 matched-pair McNemar closed-form odds-ratio engine. matcha() computes
# OR = n10/n01 with Var(log OR) = 1/n10 + 1/n01 WITHOUT survival::clogit, and
# contrast(type = "or") reports it. Oracles, in increasing independence:
#   1. survival::clogit on the same 1:1 binary data -- the conditional likelihood
#      reduces exactly to McNemar's, so the OR, SE, and CI must match it.
#   2. The closed form n10/n01 and 1/n10 + 1/n01 computed by hand from the
#      discordant-pair counts -- the variance oracle that does NOT touch clogit.
#   3. A matched-pair DGP with a known conditional log-OR (Breslow & Day 1980).
# The OR^2-bias test pins the hard-rules invariant: the unconditional 1:1 MLE
# (a parameter per pair) is exactly twice the conditional estimate.

# --- McNemar == survival::clogit on the same 1:1 binary data ---------------

test_that("McNemar reproduces survival::clogit on 1:1 binary data", {
  df <- make_matched_cc(n_sets = 300L, ratio = 1L, beta_x = log(2.5))
  res <- contrast(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set", ratio = 1),
      estimator = "mcnemar"
    ),
    type = "or"
  )
  # survival::clogit is an independent engine; its 1:1 conditional likelihood is
  # McNemar's, so the OR, the log-OR SE, and the CI must agree to numerical
  # precision -- a genuine cross-engine check, not a forwarding tautology.
  oracle <- survival::clogit(case ~ x + strata(set), data = df)
  b <- unname(stats::coef(oracle)["x"])
  se <- sqrt(stats::vcov(oracle)["x", "x"])
  z <- stats::qnorm(0.975)
  expect_equal(res$contrasts$estimate, exp(b), tolerance = 1e-7)
  expect_equal(res$estimates$se, se, tolerance = 1e-7)
  expect_equal(
    c(res$contrasts$ci_lower, res$contrasts$ci_upper),
    exp(b + c(-1, 1) * z * se),
    tolerance = 1e-7
  )
})

# --- closed-form OR AND variance, independent of clogit --------------------

test_that("the McNemar OR and variance match the discordant-pair closed form", {
  df <- make_matched_cc(n_sets = 350L, ratio = 1L, beta_x = log(3))
  # Hand-count the discordant pairs, WITHOUT any model: n10 = case exposed &
  # control unexposed, n01 = case unexposed & control exposed.
  sets <- split(df, df$set)
  n10 <- sum(vapply(
    sets,
    function(s) s$x[s$case == 1L] == 1L && s$x[s$case == 0L] == 0L,
    logical(1)
  ))
  n01 <- sum(vapply(
    sets,
    function(s) s$x[s$case == 1L] == 0L && s$x[s$case == 0L] == 1L,
    logical(1)
  ))

  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "mcnemar"
  )
  res <- contrast(fit, type = "or")
  # Point OR = n10 / n01 exactly.
  expect_equal(res$contrasts$estimate, n10 / n01, tolerance = 1e-10)
  # Var(log OR) = 1/n10 + 1/n01 exactly -- the variance oracle no clogit
  # comparison provides (it would only re-check forwarding).
  expect_equal(res$estimates$se^2, 1 / n10 + 1 / n01, tolerance = 1e-10)
  # The fitted object exposes the raw counts and pair total.
  expect_identical(fit$model$n10, n10)
  expect_identical(fit$model$n01, n01)
  expect_identical(fit$model$n_pairs, 350L)
  expect_identical(fit$model$n, 700L)
})

# --- truth-based recovery --------------------------------------------------

test_that("the McNemar OR recovers the matched-pair log-OR (truth-based)", {
  df <- make_matched_cc(n_sets = 500L, ratio = 1L, beta_x = log(2.5))
  truth <- attr(df, "truth")
  res <- contrast(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )
  # The set-level exposure prevalence is conditioned away, so the CMLE recovers
  # exp(beta_x) up to sampling error. Self-scaling band (a few reported SEs)
  # rather than a fixed tolerance, per the matched-CC variance hard rule.
  expect_lt(
    abs(log(res$contrasts$estimate) - unname(truth["beta_x"])),
    3.5 * res$estimates$se
  )
})

# --- the OR^2-bias invariant (unconditional MLE = 2x conditional) ----------

test_that("the unconditional 1:1 MLE doubles the conditional log-OR (OR^2 bias)", {
  # Hard-rules invariant: for 1:1 matching, unconditional logistic regression
  # with a parameter per pair estimates 2*beta (so the OR is squared); the
  # conditional / McNemar estimate is beta. Concordant pairs are fit perfectly
  # by their pair intercept and drop out, leaving the discordant pairs, on which
  # the unconditional log-OR is algebraically exactly twice the conditional one.
  df <- make_matched_cc(n_sets = 400L, ratio = 1L, beta_x = log(2.5))
  beta_cond <- contrast(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )$estimates$estimate

  unconditional <- suppressWarnings(stats::glm(
    case ~ x + factor(set),
    data = df,
    family = stats::binomial()
  ))
  beta_uncond <- unname(stats::coef(unconditional)["x"])

  # Exactly double (to glm's convergence tolerance) -- the OR is squared.
  expect_equal(beta_uncond, 2 * beta_cond, tolerance = 1e-3)
  # And the bias is large, not a rounding artefact: the unconditional estimate
  # is far from the true beta (a sign the matched-set-indicator fit is wrong).
  expect_gt(abs(beta_uncond - beta_cond), 0.5 * abs(beta_cond))
})

# --- exposure coding equivalence -------------------------------------------

test_that("logical and two-level-factor exposures give the same McNemar OR", {
  df <- make_matched_cc(n_sets = 300L, ratio = 1L)
  or_num <- contrast(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )$contrasts$estimate

  df$xl <- df$x == 1L # logical
  df$xf <- factor(
    ifelse(df$x == 1L, "exposed", "unexposed"), # 2-level factor
    levels = c("unexposed", "exposed")
  )
  or_log <- contrast(
    matcha(df, "case", "xl", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )$contrasts$estimate
  or_fac <- contrast(
    matcha(df, "case", "xf", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )$contrasts$estimate
  expect_equal(or_log, or_num, tolerance = 1e-12)
  expect_equal(or_fac, or_num, tolerance = 1e-12)
})

# --- concordant pairs are inert --------------------------------------------

test_that("a pair concordant on exposure leaves the OR and SE unchanged", {
  df <- make_matched_cc(n_sets = 200L, ratio = 1L)
  base <- contrast(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    type = "or"
  )
  # A 1:1 pair whose case and control share the same exposure has no within-pair
  # exposure contrast, so it is neither n10 nor n01 and cannot move the estimate.
  concordant <- data.frame(
    case = c(1L, 0L),
    x = c(1L, 1L),
    z = 0L,
    set = max(df$set) + 1L
  )
  aug <- contrast(
    matcha(
      rbind(df, concordant),
      "case",
      "x",
      matched_cc(strata = "set"),
      estimator = "mcnemar"
    ),
    type = "or"
  )
  expect_equal(
    aug$contrasts$estimate,
    base$contrasts$estimate,
    tolerance = 1e-12
  )
  expect_equal(aug$estimates$se, base$estimates$se, tolerance = 1e-12)
  # The concordant pair IS counted among the pairs analysed, though.
  expect_identical(aug$n, base$n + 2L)
})

# --- result structure / labels ---------------------------------------------

test_that("the McNemar result is labelled and sized correctly", {
  df <- make_matched_cc(n_sets = 250L, ratio = 1L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "mcnemar"
  )
  res <- contrast(fit, type = "or")
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "McNemar OR")
  expect_identical(res$estimator, "mcnemar")
  expect_identical(res$engine, "mcnemar")
  expect_identical(res$contrasts$comparison, "x")
  # n is the individuals in the complete 1:1 pairs (2 per pair).
  expect_identical(res$n, fit$model$n)
  # contrast() with no type defaults to the OR for the McNemar engine.
  expect_identical(contrast(fit)$type, "or")
  # The log-scale estimates reconstruct the CI; the OR-scale se is delta-method.
  z <- stats::qnorm(0.975)
  expect_equal(
    res$contrasts$ci_lower,
    exp(res$estimates$estimate - z * res$estimates$se),
    tolerance = 1e-12
  )
})

# --- rejections ------------------------------------------------------------

test_that("M:1 (or richer) matching is rejected for McNemar", {
  df <- make_matched_cc(n_sets = 80L, ratio = 3L)
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    class = "matchatr_not_one_to_one"
  )
  # A sample mixing 1:1 and 1:2 sets is still not pure 1:1 -> rejected.
  mixed <- rbind(
    make_matched_cc(n_sets = 50L, ratio = 1L),
    transform(make_matched_cc(n_sets = 50L, ratio = 2L), set = set + 1000L)
  )
  expect_error(
    matcha(
      mixed,
      "case",
      "x",
      matched_cc(strata = "set"),
      estimator = "mcnemar"
    ),
    class = "matchatr_not_one_to_one"
  )
})

test_that("a non-binary exposure is rejected by McNemar (points to clogit)", {
  df <- make_matched_cc(n_sets = 100L, ratio = 1L)
  df$x <- factor(sample(c("lo", "mid", "hi"), nrow(df), replace = TRUE))
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    class = "matchatr_bad_input"
  )
})

test_that("one-sided (or absent) discordant pairs are not estimable", {
  df <- make_matched_cc(n_sets = 150L, ratio = 1L)
  # Force the exposure to equal the case indicator: every discordant pair is
  # case-exposed/control-unexposed (n01 = 0), so OR = Inf is not identified.
  df$x <- df$case
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("RD / RR are rejected as unidentified for McNemar", {
  df <- make_matched_cc(n_sets = 200L, ratio = 1L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "mcnemar"
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

test_that("sandwich / bootstrap CIs are not available for McNemar", {
  df <- make_matched_cc(n_sets = 200L, ratio = 1L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "mcnemar"
  )
  expect_error(
    contrast(fit, type = "or", ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
  expect_error(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("a missing pair member triggers the dropped-rows warning with a count", {
  df <- make_matched_cc(n_sets = 200L, ratio = 1L)
  df$x[1:3] <- NA
  expect_warning(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    class = "matchatr_dropped_rows",
    regexp = "3 row"
  )
})

test_that("McNemar rejection messages read clearly", {
  df3 <- make_matched_cc(n_sets = 80L, ratio = 3L)
  expect_snapshot(
    matcha(df3, "case", "x", matched_cc(strata = "set"), estimator = "mcnemar"),
    error = TRUE
  )
  df <- make_matched_cc(n_sets = 200L, ratio = 1L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "mcnemar"
  )
  expect_snapshot(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    error = TRUE
  )
})
