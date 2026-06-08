# Homogeneity of an exposure's odds ratios across disease subtypes, from a
# polytomous (multinomial) case-control fit.
#
# Oracles, strongest first:
#  1. Exact linear algebra. The reported chi-squared and pooled OR equal the
#     hand-built C V C' Wald form and the GLS combination computed from the SAME
#     multinom variance the contrast layer exposes -- pins the implementation's
#     algebra to its definition to machine precision.
#  2. Closed-form 2x2 Woolf oracle, INDEPENDENT of multinom's vcov. For a binary
#     exposure and a saturated 3-group outcome the two subtype log-ORs and their
#     covariance have a closed form (Woolf variances on the diagonal, the shared
#     reference cells 1/n_ref1 + 1/n_ref0 off-diagonal). The reference cells
#     cancel in the difference variance, so df = 1 and the homogeneity
#     chi-squared is (b_A - b_B)^2 / (1/n_A1 + 1/n_A0 + 1/n_B1 + 1/n_B0).
#  3. riskclustr::eh_test_subtype (mlogit engine, independent codebase) -- the
#     canonical applied implementation of this exact test.
#  4. Operating characteristics: a truth DGP with EQUAL subtype ORs gives a test
#     of correct size; with UNEQUAL ORs, high power; and under a true common OR
#     the pooled SE is smaller than each subtype SE (Begg & Gray, 1984).

# ---- Exact linear-algebra oracle --------------------------------------------

test_that("chi-squared and pooled OR equal the hand-built C V C' / GLS forms", {
  d <- make_3group_table()
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  h <- test_homogeneity(fit)

  # The contrast layer exposes the same stacked subtype log-ORs and covariance
  # the homogeneity test consumes; rebuild both statistics by hand from them.
  res <- contrast(fit, type = "or")
  b <- res$estimates$estimate # caseA: x, caseB: x log-ORs
  V <- res$vcov
  cmat <- matrix(c(1, -1), nrow = 1)
  cb <- cmat %*% b
  chisq_hand <- as.numeric(t(cb) %*% solve(cmat %*% V %*% t(cmat)) %*% cb)
  vinv <- solve(V)
  ones <- c(1, 1)
  bc_hand <- as.numeric(t(ones) %*% vinv %*% b) /
    as.numeric(t(ones) %*% vinv %*% ones)
  se_hand <- sqrt(1 / as.numeric(t(ones) %*% vinv %*% ones))

  expect_equal(h$homogeneity$statistic, chisq_hand, tolerance = 1e-10)
  expect_equal(h$homogeneity$df, 1L)
  expect_equal(h$homogeneity$common_or, exp(bc_hand), tolerance = 1e-10)
  z <- stats::qnorm(0.975)
  expect_equal(
    h$homogeneity$ci_lower,
    exp(bc_hand - z * se_hand),
    tolerance = 1e-10
  )
  expect_equal(
    h$homogeneity$ci_upper,
    exp(bc_hand + z * se_hand),
    tolerance = 1e-10
  )
})

# ---- Closed-form Woolf oracle (independent of multinom's vcov) ---------------

test_that("homogeneity chi-squared and pooled OR match the closed-form 2x2 Woolf values", {
  d <- make_3group_table() # control 80/120, caseA 60/40, caseB 30/70
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  h <- test_homogeneity(fit)

  # Closed-form subtype log-ORs and their Woolf variance / covariance.
  b_a <- log((60 * 120) / (40 * 80))
  b_b <- log((30 * 120) / (70 * 80))
  v_a <- 1 / 60 + 1 / 40 + 1 / 80 + 1 / 120
  v_b <- 1 / 30 + 1 / 70 + 1 / 80 + 1 / 120
  cov_ab <- 1 / 80 + 1 / 120 # shared reference cells
  # The reference cells cancel in the difference variance: df = 1.
  var_diff <- v_a + v_b - 2 * cov_ab
  chisq_cf <- (b_a - b_b)^2 / var_diff
  p_cf <- stats::pchisq(chisq_cf, df = 1, lower.tail = FALSE)

  V <- matrix(c(v_a, cov_ab, cov_ab, v_b), 2, 2)
  vinv <- solve(V)
  ones <- c(1, 1)
  denom <- as.numeric(t(ones) %*% vinv %*% ones)
  bc_cf <- as.numeric(t(ones) %*% vinv %*% c(b_a, b_b)) / denom

  # The relative gap is bounded by multinom's optimiser tolerance (the saturated
  # multinomial vcov equals the Woolf form only up to convergence).
  expect_equal(h$homogeneity$statistic, chisq_cf, tolerance = 1e-2)
  expect_equal(h$homogeneity$p.value, p_cf, tolerance = 1e-2)
  expect_equal(h$homogeneity$common_or, exp(bc_cf), tolerance = 1e-3)

  # And the engine's exposure vcov sub-block IS the closed-form V.
  ex <- multinom_exposure_or(fit$model, "x")
  expect_equal(unname(ex$vcov), V, tolerance = 1e-4)
})

# ---- Closed-form oracle with M = 3 subtypes (df = 2) -------------------------

test_that("the df = 2 homogeneity test matches the closed-form 4-group Woolf values", {
  # Exercises the general (M-1) x M contrast and the >=2x2 C V C' inversion that
  # the binary df = 1 cases never reach. caseA / caseB share OR = 2.25, caseC has
  # OR = 0.5, so the homogeneity test rejects.
  d <- make_4group_table()
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  h <- test_homogeneity(fit)
  expect_equal(nrow(h$homogeneity), 1L) # one exposure column
  expect_equal(h$homogeneity$df, 2L)

  # Closed-form subtype log-ORs and the 3x3 Woolf covariance: the diagonal is the
  # per-subtype Woolf sum, every off-diagonal the shared reference cells.
  b <- c(
    log((60 * 120) / (40 * 80)),
    log((45 * 120) / (30 * 80)),
    log((20 * 120) / (60 * 80))
  )
  v_a <- 1 / 60 + 1 / 40 + 1 / 80 + 1 / 120
  v_b <- 1 / 45 + 1 / 30 + 1 / 80 + 1 / 120
  v_c <- 1 / 20 + 1 / 60 + 1 / 80 + 1 / 120
  off <- 1 / 80 + 1 / 120
  V <- matrix(off, 3, 3)
  diag(V) <- c(v_a, v_b, v_c)
  cmat <- rbind(c(1, -1, 0), c(0, 1, -1))
  cb <- cmat %*% b
  chisq_cf <- as.numeric(t(cb) %*% solve(cmat %*% V %*% t(cmat)) %*% cb)
  p_cf <- stats::pchisq(chisq_cf, df = 2, lower.tail = FALSE)

  expect_equal(h$homogeneity$statistic, chisq_cf, tolerance = 1e-2)
  expect_equal(h$homogeneity$p.value, p_cf, tolerance = 1e-2)

  # Exact: the same statistic rebuilt from multinom's own 3x3 exposure vcov.
  res <- contrast(fit, type = "or")
  cb2 <- cmat %*% res$estimates$estimate
  chisq_hand <- as.numeric(
    t(cb2) %*% solve(cmat %*% res$vcov %*% t(cmat)) %*% cb2
  )
  expect_equal(h$homogeneity$statistic, chisq_hand, tolerance = 1e-10)
})

# ---- Definition consistency on an adjusted (continuous-confounder) fit -------

test_that("the pooled OR and chi-squared are the GLS / Wald functionals of contrast()", {
  d <- make_polytomous_cc(n = 4000L, seed = 23L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  res <- contrast(fit, type = "or")

  b <- res$estimates$estimate
  V <- res$vcov
  cmat <- matrix(c(1, -1), nrow = 1)
  cb <- cmat %*% b
  chisq_hand <- as.numeric(t(cb) %*% solve(cmat %*% V %*% t(cmat)) %*% cb)
  vinv <- solve(V)
  ones <- c(1, 1)
  bc_hand <- as.numeric(t(ones) %*% vinv %*% b) /
    as.numeric(t(ones) %*% vinv %*% ones)

  expect_equal(h$homogeneity$statistic, chisq_hand, tolerance = 1e-10)
  expect_equal(h$homogeneity$common_or, exp(bc_hand), tolerance = 1e-10)
  # Continuous confounders are handled directly (no constrained refit needed).
  expect_equal(h$n, nrow(d))
  expect_equal(h$reference, "control")
})

# ---- riskclustr external oracle (mlogit engine, independent codebase) --------

test_that("the heterogeneity p-value matches riskclustr::eh_test_subtype", {
  skip_if_not_installed("riskclustr")
  # riskclustr codes the subtype variable as 0 = control, 1..M = subtypes, and
  # fits the multinomial via mlogit (a different engine), so agreement validates
  # the Wald p-value against an independent implementation. A moderate-signal
  # sample keeps the p-value away from the underflow floor for a meaningful
  # comparison.
  d <- make_polytomous_cc(
    n = 1200L,
    beta_a = log(1.6),
    beta_b = log(2.4),
    seed = 31L
  )
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)

  d$subtype <- match(as.character(d$g), c("control", "caseA", "caseB")) - 1L
  oracle <- riskclustr::eh_test_subtype(
    label = "subtype",
    M = 2,
    factors = list("x"),
    data = as.data.frame(d)
  )
  # riskclustr returns `eh_pval` as a one-row data.frame (column `p_het`, row
  # named by the risk factor) and `beta` as a [risk factor x subtype] data.frame.
  # Compare on the chi-squared scale (robust to the p-value magnitude): back out
  # riskclustr's statistic from its heterogeneity p-value (df = M - 1 = 1).
  chisq_oracle <- stats::qchisq(
    oracle$eh_pval["x", "p_het"],
    df = 1,
    lower.tail = FALSE
  )
  expect_equal(h$homogeneity$statistic, chisq_oracle, tolerance = 1e-2)
  # And the per-subtype log-ORs agree with riskclustr's mlogit fit (different
  # engine), confirming the inputs to the test -- not just the test statistic.
  res <- contrast(fit, type = "or")
  expect_equal(
    res$estimates$estimate,
    unname(unlist(oracle$beta["x", ])),
    tolerance = 1e-2
  )
})

# ---- Operating characteristics ----------------------------------------------

test_that("the test has correct size under equal subtype ORs (continuous confounder)", {
  # Under H0 (caseA and caseB share the exposure OR), the Wald p-value should be
  # ~Uniform, so the rejection rate at alpha = 0.05 sits near 0.05. Deterministic
  # via the fixed seed; the band absorbs Monte Carlo noise (R = 200).
  reps <- 200L
  rej <- withr::with_seed(2024L, {
    sum(vapply(
      seq_len(reps),
      function(i) {
        di <- make_polytomous_cc(
          n = 700L,
          beta_a = log(2),
          beta_b = log(2),
          seed = sample.int(1e6, 1)
        )
        fi <- matcha(
          di,
          outcome = "g",
          exposure = "x",
          design = unmatched_cc(),
          confounders = ~age,
          estimator = "polytomous"
        )
        test_homogeneity(fi)$homogeneity$p.value < 0.05
      },
      logical(1)
    ))
  })
  expect_gte(rej / reps, 0.02)
  expect_lte(rej / reps, 0.10)
})

test_that("the test has power against unequal subtype ORs", {
  # caseA OR = 2.5, caseB OR = 0.5 -- strongly heterogeneous.
  d <- make_polytomous_cc(n = 6000L, seed = 23L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  expect_lt(h$homogeneity$p.value, 1e-3)
})

test_that("the pooled OR is more efficient than each subtype OR (Begg & Gray)", {
  # Under a true common OR the GLS pooled estimate has a smaller SE than any
  # single subtype OR.
  d <- make_polytomous_cc(
    n = 8000L,
    beta_a = log(2),
    beta_b = log(2),
    seed = 7L
  )
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  res <- contrast(fit, type = "or")
  z <- stats::qnorm(0.975)
  se_common <- (log(h$homogeneity$ci_upper) - log(h$homogeneity$common_or)) / z
  expect_lt(se_common, min(res$estimates$se))
})

# ---- Factor exposure (one homogeneity test per level) -----------------------

test_that("a factor exposure reports a per-level homogeneity test vs a multinom contrast", {
  d <- make_polytomous_cc(n = 4000L, seed = 12L)
  d$xf <- factor(ifelse(d$x == 1L, "hi", "lo"))
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "xf",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  # "hi" is the alphabetical reference, so the single exposure column is "xflo".
  expect_equal(h$homogeneity$term, "xflo")

  # Hand-built Wald on a raw multinom: the xflo coefficients for caseA/caseB and
  # their vcov sub-block. matcha fits the same formula, so this is exact.
  oracle <- nnet::multinom(g ~ xf + age, data = d, trace = FALSE)
  b <- c(coef(oracle)["caseA", "xflo"], coef(oracle)["caseB", "xflo"])
  vc <- vcov(oracle)
  nm <- c("caseA:xflo", "caseB:xflo")
  V <- vc[nm, nm]
  cmat <- matrix(c(1, -1), nrow = 1)
  cb <- cmat %*% b
  chisq_hand <- as.numeric(t(cb) %*% solve(cmat %*% V %*% t(cmat)) %*% cb)
  expect_equal(h$homogeneity$statistic, chisq_hand, tolerance = 1e-6)
})

test_that("a 3-level factor exposure tests each level separately (multi-column stride)", {
  # Two exposure columns (xf3b, xf3c) drive the subtype-major / column-minor
  # regrouping stride that a single-column exposure never exercises. Level "b"
  # raises both subtypes equally (homogeneous); level "c" raises caseA but lowers
  # caseB (heterogeneous), so the two columns must yield DISTINCT statistics -- a
  # transposed stride would scramble them and fail the per-column oracle.
  d <- withr::with_seed(7L, {
    n <- 5000L
    xf3 <- sample(c("a", "b", "c"), n, replace = TRUE)
    b_a <- c(a = 0, b = 0.6, c = 0.6)[xf3]
    b_b <- c(a = 0, b = 0.6, c = -0.6)[xf3]
    lp_a <- -0.5 + b_a
    lp_b <- -0.5 + b_b
    denom <- 1 + exp(lp_a) + exp(lp_b)
    u <- runif(n)
    g <- ifelse(
      u < 1 / denom,
      "control",
      ifelse(u < (1 + exp(lp_a)) / denom, "caseA", "caseB")
    )
    data.frame(
      g = factor(g, levels = c("control", "caseA", "caseB")),
      xf3 = factor(xf3, levels = c("a", "b", "c"))
    )
  })
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "xf3",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  # "a" is the factor reference, so the two columns are xf3b and xf3c.
  expect_equal(nrow(h$homogeneity), 2L)
  expect_setequal(h$homogeneity$term, c("xf3b", "xf3c"))
  # The heterogeneous level ("c") must give a far larger statistic than the
  # homogeneous one ("b") -- confirms the columns are not swapped.
  stat_b <- h$homogeneity$statistic[h$homogeneity$term == "xf3b"]
  stat_c <- h$homogeneity$statistic[h$homogeneity$term == "xf3c"]
  expect_lt(stat_b, stat_c)

  # Each column's statistic equals a hand-built Wald on the SAME-named multinom
  # coefficients -- pins the per-column regrouping to the exact subtype pair.
  oracle <- nnet::multinom(g ~ xf3, data = d, trace = FALSE)
  vc <- vcov(oracle)
  cmat <- matrix(c(1, -1), nrow = 1)
  for (lv in c("xf3b", "xf3c")) {
    b <- c(coef(oracle)["caseA", lv], coef(oracle)["caseB", lv])
    nm <- c(paste0("caseA:", lv), paste0("caseB:", lv))
    V <- vc[nm, nm]
    cb <- cmat %*% b
    chisq_hand <- as.numeric(t(cb) %*% solve(cmat %*% V %*% t(cmat)) %*% cb)
    expect_equal(
      h$homogeneity$statistic[h$homogeneity$term == lv],
      chisq_hand,
      tolerance = 1e-6
    )
  }
})

# ---- S3 surface -------------------------------------------------------------

test_that("print shows the common OR table and a Wald-test header", {
  d <- make_3group_table()
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  out <- paste(
    utils::capture.output(print(test_homogeneity(fit))),
    collapse = "\n"
  )
  expect_match(out, "Homogeneity of subtype odds ratios \\(Wald\\)")
  expect_match(out, "Reference:  control")
  expect_match(out, "common_or")
})

test_that("tidy returns one row per exposure term with the test columns", {
  d <- make_polytomous_cc(n = 2500L, seed = 9L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  h <- test_homogeneity(fit)
  td <- tidy(h)
  expect_true(data.table::is.data.table(td))
  expect_equal(nrow(td), 1L)
  # Broom column convention, shared with tidy.matchatr_result.
  expect_setequal(
    names(td),
    c(
      "term",
      "estimate",
      "std.error",
      "conf.low",
      "conf.high",
      "statistic",
      "df",
      "p.value"
    )
  )
  expect_identical(td$term, "x")
  # `estimate` is the common OR; `std.error` the OR-scale delta-method SE.
  expect_equal(td$estimate, h$homogeneity$common_or)
  expect_equal(td$std.error, h$homogeneity$se)
})

# ---- Rejection paths --------------------------------------------------------

test_that("homogeneity is rejected for a non-polytomous (binary) fit", {
  d <- data.frame(case = rep(c(1, 0), 100), x = rbinom(200, 1, 0.4))
  fit <- matcha(d, outcome = "case", exposure = "x", design = unmatched_cc())
  expect_snapshot(error = TRUE, test_homogeneity(fit))
})

test_that("homogeneity is rejected for the Mantel-Haenszel and clogit engines", {
  dm <- make_stratified_cc(seed = 4L)
  fm <- matcha(
    dm,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "mh"
  )
  expect_error(test_homogeneity(fm), class = "matchatr_bad_input")

  dc <- make_matched_cc(n_sets = 40L, seed = 2L)
  fc <- matcha(
    dc,
    outcome = "case",
    exposure = "x",
    design = matched_cc(strata = "set"),
    estimator = "clogit"
  )
  expect_error(test_homogeneity(fc), class = "matchatr_bad_input")
})

test_that("homogeneity is rejected for a fit with no estimated model", {
  df <- data.frame(
    case = c(1, 0, 1, 0),
    x = c(1, 0, 1, 0),
    set = c(1, 1, 2, 2),
    t = c(2, 3, 1, 4)
  )
  fit <- matcha(df, "case", "x", counter_matched(strata = "set", time = "t"))
  expect_error(test_homogeneity(fit), class = "matchatr_not_estimated")
})

test_that("a malformed confidence level and a non-fit object are rejected", {
  d <- make_polytomous_cc(n = 600L, seed = 10L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  expect_error(
    test_homogeneity(fit, conf_level = 1.5),
    class = "matchatr_bad_input"
  )
  expect_error(test_homogeneity(list(a = 1)), class = "matchatr_bad_input")
})

test_that("a singular subtype covariance is matchatr_unestimable_exposure", {
  # Direct unit test of the per-term guard: a rank-1 covariance (perfectly
  # correlated subtype log-ORs) makes both the contrast covariance C V C' and V
  # singular, so the pooled OR / test are undefined. Unreachable from a converged
  # multinom (collinearity is rejected at fit time), so it is exercised here at
  # the kernel.
  expect_error(
    homogeneity_one_term(
      beta = c(0.5, 0.5),
      vmat = matrix(c(1, 1, 1, 1), 2, 2),
      term = "x",
      z = 1.96
    ),
    class = "matchatr_unestimable_exposure"
  )
})
