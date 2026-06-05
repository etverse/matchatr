# Polytomous (multinomial) logistic for multiple case / control groups.
#
# Oracles, strongest first:
#  1. Closed-form saturated multinomial == 2x2 Woolf log-OR AND variance. For a
#     binary exposure and a 3-group outcome with no confounder, the
#     baseline-category logit is saturated, so each non-reference equation's
#     exposure coefficient is the log OR of the {reference, group} x {0, 1}
#     subtable and its variance is the Woolf sum 1/a + 1/b + 1/c + 1/d. This
#     pins BOTH the point and the SE independently of nnet::multinom's own vcov
#     (the multinom-fidelity checks below only verify forwarding).
#  2. Truth-based DGP recovery: a cohort drawn from a multinomial with KNOWN
#     per-subtype exposure log-ORs; the estimates fall within a few SE of truth.
#  3. nnet::multinom fidelity: matcha's coef / vcov equal a hand-built multinom.

# A deterministic 3-group case-control table with exact cell counts, so the
# saturated multinomial coefficients and variances have a closed form. Rows:
# control (reference), caseA, caseB; columns: x = 0 / 1.
make_3group_table <- function(
  ctrl1 = 80L,
  ctrl0 = 120L,
  a1 = 60L,
  a0 = 40L,
  b1 = 30L,
  b0 = 70L
) {
  g <- c(
    rep("control", ctrl1 + ctrl0),
    rep("caseA", a1 + a0),
    rep("caseB", b1 + b0)
  )
  x <- c(
    rep(1L, ctrl1),
    rep(0L, ctrl0),
    rep(1L, a1),
    rep(0L, a0),
    rep(1L, b1),
    rep(0L, b0)
  )
  data.frame(
    g = factor(g, levels = c("control", "caseA", "caseB")),
    x = x,
    stringsAsFactors = FALSE
  )
}

# ---- Closed-form 2x2 Woolf oracle (independent point AND variance) ----------

test_that("saturated multinomial OR == closed-form 2x2 Woolf odds ratio", {
  d <- make_3group_table()
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  res <- contrast(fit, type = "or")

  # Closed form per non-reference group k: OR = (n_{k,1} n_{ref,0}) /
  # (n_{k,0} n_{ref,1}); Var(log OR) = 1/n_{k,1} + 1/n_{k,0} + 1/n_{ref,1} +
  # 1/n_{ref,0} (Woolf) over the {control, k} x {0, 1} subtable.
  woolf <- function(k1, k0, r1 = 80, r0 = 120) {
    list(
      or = (k1 * r0) / (k0 * r1),
      se_log = sqrt(1 / k1 + 1 / k0 + 1 / r1 + 1 / r0)
    )
  }
  oa <- woolf(60, 40)
  ob <- woolf(30, 70)

  # contrasts are returned in subtype order (caseA, caseB).
  expect_equal(res$contrasts$comparison, c("caseA: x", "caseB: x"))
  expect_equal(res$contrasts$estimate, c(oa$or, ob$or), tolerance = 1e-4)
  # The log-scale SE lives in $estimates (the reconstructable one); pin it to the
  # Woolf variance, an oracle independent of multinom's information matrix.
  expect_equal(res$estimates$se, c(oa$se_log, ob$se_log), tolerance = 1e-4)

  # The Wald interval is exp(logOR +/- z * se_log); reconstruct it from the
  # closed-form pieces to confirm the bounds are not multinom-derived.
  z <- stats::qnorm(0.975)
  expect_equal(
    res$contrasts$ci_lower,
    c(oa$or, ob$or) * exp(-z * c(oa$se_log, ob$se_log)),
    tolerance = 1e-4
  )
  expect_equal(res$n, nrow(d))
})

# ---- Truth-based DGP recovery -----------------------------------------------

test_that("multinomial recovers known per-subtype exposure log-ORs", {
  d <- make_polytomous_cc(n = 8000L, seed = 23L)
  truth <- attr(d, "truth")
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  res <- contrast(fit, type = "or")
  # estimates carry the log-OR point + SE; each subtype's estimate must fall
  # within ~3.5 SE of its data-generating log-OR (SE-scaled band, not a fixed
  # absolute tolerance).
  est <- res$estimates
  expect_equal(est$term, c("caseA: x", "caseB: x"))
  expect_lt(abs(est$estimate[1] - truth[["caseA.x"]]), 3.5 * est$se[1])
  expect_lt(abs(est$estimate[2] - truth[["caseB.x"]]), 3.5 * est$se[2])
})

# ---- nnet::multinom fidelity (coef / vcov forwarding) -----------------------

test_that("matcha forwards nnet::multinom coefficients and variance exactly", {
  d <- make_polytomous_cc(n = 3000L, seed = 5L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  # Hand-built oracle: the same formula on the same reference-first factor.
  oracle <- nnet::multinom(g ~ x + age, data = d, trace = FALSE)
  expect_equal(unname(coef(fit$model)), unname(coef(oracle)), tolerance = 1e-6)
  expect_equal(unname(vcov(fit$model)), unname(vcov(oracle)), tolerance = 1e-6)

  # The reported per-subtype log-OR / SE are exactly the oracle's exposure
  # coefficients and their vcov-diagonal SEs.
  res <- contrast(fit, type = "or")
  cf <- coef(oracle)
  vc <- vcov(oracle)
  expect_equal(
    res$estimates$estimate,
    c(cf["caseA", "x"], cf["caseB", "x"]),
    tolerance = 1e-8
  )
  expect_equal(
    res$estimates$se,
    c(sqrt(vc["caseA:x", "caseA:x"]), sqrt(vc["caseB:x", "caseB:x"])),
    tolerance = 1e-8
  )
})

test_that("continuous exposure reports a per-unit subtype OR vs multinom", {
  d <- make_polytomous_cc(n = 3000L, seed = 8L)
  # Use age as a continuous exposure, x as the confounder.
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "age",
    design = unmatched_cc(),
    confounders = ~x,
    estimator = "polytomous"
  )
  oracle <- nnet::multinom(g ~ age + x, data = d, trace = FALSE)
  res <- contrast(fit, type = "or")
  cf <- coef(oracle)
  expect_equal(
    res$contrasts$estimate,
    c(exp(cf["caseA", "age"]), exp(cf["caseB", "age"])),
    tolerance = 1e-6
  )
})

test_that("factor exposure reports per-level subtype ORs vs multinom", {
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
  res <- contrast(fit, type = "or")
  # "hi" is the alphabetical reference, so each comparison is the "lo" level.
  expect_equal(res$contrasts$comparison, c("caseA: xflo", "caseB: xflo"))
  oracle <- nnet::multinom(g ~ xf + age, data = d, trace = FALSE)
  cf <- coef(oracle)
  expect_equal(
    res$estimates$estimate,
    c(cf["caseA", "xflo"], cf["caseB", "xflo"]),
    tolerance = 1e-6
  )
})

# ---- Reference handling -----------------------------------------------------

test_that("reference choice sets the baseline and matches a releveled multinom", {
  d <- make_polytomous_cc(n = 3000L, seed = 15L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "caseA"
  )
  # The fitted model's first outcome level is the requested reference.
  expect_identical(fit$model$lev[1], "caseA")
  expect_identical(fit$details$reference, "caseA")
  # Equality with an explicitly releveled multinom (the only difference vs the
  # default-reference fit is which group is the baseline).
  d2 <- d
  d2$g <- stats::relevel(d2$g, ref = "caseA")
  oracle <- nnet::multinom(g ~ x, data = d2, trace = FALSE)
  expect_equal(unname(coef(fit$model)), unname(coef(oracle)), tolerance = 1e-6)
  # The two non-reference equations are now control and caseB.
  expect_setequal(rownames(coef(fit$model)), c("control", "caseB"))
})

test_that("default (NULL) reference is the first level; character outcome sorts", {
  d <- make_polytomous_cc(n = 2000L, seed = 3L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  expect_identical(fit$details$reference, "control") # first declared level

  # A character outcome has no declared order, so factor() sorts: caseA first.
  dc <- d
  dc$g <- as.character(dc$g)
  fitc <- matcha(
    dc,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  expect_identical(fitc$details$reference, "caseA")
})

test_that("an unused outcome level is dropped before counting groups", {
  d <- make_polytomous_cc(n = 1500L, seed = 6L)
  # Declare a fourth level that never occurs; droplevels must remove it so the
  # group count and the fitted equations reflect only observed groups.
  levels(d$g) <- c(levels(d$g), "caseC")
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  expect_setequal(fit$details$group_levels, c("control", "caseA", "caseB"))
  expect_false("caseC" %in% rownames(coef(fit$model)))
})

# ---- tidy / summary / print -------------------------------------------------

test_that("tidy() returns a per-equation table with a y.level column", {
  d <- make_polytomous_cc(n = 2500L, seed = 9L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  td <- tidy(fit, exponentiate = TRUE)
  expect_true(all(c("y.level", "term", "estimate") %in% names(td)))
  # One row per (non-reference level) x (intercept + x + age) = 2 x 3.
  expect_equal(nrow(td), 6L)
  expect_setequal(unique(td$y.level), c("caseA", "caseB"))
  # Exponentiated estimate equals exp(multinom coef) for the exposure rows.
  oracle <- nnet::multinom(g ~ x + age, data = d, trace = FALSE)
  xrow <- td[td$term == "x" & td$y.level == "caseA", ]
  expect_equal(xrow$estimate, exp(coef(oracle)["caseA", "x"]), tolerance = 1e-6)
})

test_that("print shows per-group counts with the reference flagged", {
  d <- make_3group_table()
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  out <- paste(utils::capture.output(print(fit)), collapse = "\n")
  expect_match(out, "groups:")
  expect_match(out, "control")
  expect_match(out, "reference group")
  # summary() prints the exponentiated per-equation table without error.
  expect_output(summary(fit), "engine: multinom")
})

# ---- Rejection paths --------------------------------------------------------

test_that("a two-group outcome is rejected (needs >=3 groups)", {
  d <- make_polytomous_cc(n = 600L, seed = 1L)
  d$g <- factor(ifelse(d$g == "control", "control", "case"))
  expect_snapshot(
    error = TRUE,
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous"
    )
  )
})

test_that("a numeric / logical outcome is rejected by the polytomous estimator", {
  d <- data.frame(g = rep(c(0L, 1L), 100L), x = rbinom(200, 1, 0.4))
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous"
    ),
    class = "matchatr_bad_outcome"
  )
})

test_that("an out-of-range reference is rejected", {
  d <- make_polytomous_cc(n = 600L, seed = 2L)
  expect_snapshot(
    error = TRUE,
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous",
      reference = "nope"
    )
  )
})

test_that("reference on a non-polytomous estimator is rejected", {
  d <- data.frame(case = rep(c(1, 0), 100), x = rbinom(200, 1, 0.4))
  expect_snapshot(
    error = TRUE,
    matcha(
      d,
      outcome = "case",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "logistic",
      reference = "0"
    )
  )
})

test_that("a constant exposure is unestimable", {
  d <- make_polytomous_cc(n = 600L, seed = 4L)
  d$x <- 1L
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous"
    ),
    class = "matchatr_unestimable_exposure"
  )
})

# Regression: an exposure collinear with a confounder. nnet::multinom does not
# alias the redundant column to NA (unlike glm); it splits the coefficient and
# would report a silently halved OR. The design-matrix rank guard
# (reject_collinear_exposure()) must reject it, while confounder-only
# collinearity (the exposure still estimable) must NOT be rejected.
# 2026-06-03 critical-review-loop Issue #1; repro /tmp/matchatr_repro_collinear.R
test_that("an exposure collinear with a confounder is unestimable", {
  d <- make_polytomous_cc(n = 1500L, seed = 4L)
  d$dup <- d$x # confounder identical to the exposure -> rank-deficient design
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      confounders = ~dup,
      estimator = "polytomous"
    ),
    class = "matchatr_unestimable_exposure"
  )

  # Collinearity confined to the confounders leaves the exposure estimable, so
  # it must fit and recover the exposure OR rather than over-reject. The guard
  # rejects only when the exposure itself loses rank.
  z1 <- withr::with_seed(99L, stats::rnorm(nrow(d)))
  d$z1 <- z1
  d$z2 <- z1 # two confounders collinear with each other, but x is orthogonal
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~ z1 + z2,
    estimator = "polytomous"
  )
  res <- contrast(fit, type = "or")
  oracle <- nnet::multinom(g ~ x + z1 + z2, data = d, trace = FALSE)
  cf <- coef(oracle)
  expect_equal(
    res$estimates$estimate,
    c(cf["caseA", "x"], cf["caseB", "x"]),
    tolerance = 1e-6
  )
})

test_that("an ordered-factor exposure is rejected", {
  d <- make_polytomous_cc(n = 600L, seed = 7L)
  d$xo <- ordered(ifelse(d$x == 1L, "hi", "lo"), levels = c("lo", "hi"))
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "xo",
      design = unmatched_cc(),
      estimator = "polytomous"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("RD / RR are unidentified; sandwich / bootstrap variance unsupported", {
  d <- make_polytomous_cc(n = 800L, seed = 10L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    estimator = "polytomous"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
  expect_error(
    contrast(fit, type = "ratio"),
    class = "matchatr_unidentified_estimand"
  )
  expect_error(
    contrast(fit, type = "or", ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
  expect_error(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
  expect_error(
    tidy(fit, robust = TRUE),
    class = "matchatr_unsupported_variance"
  )
})

test_that("effect_modifier is rejected for the polytomous estimator", {
  d <- make_polytomous_cc(n = 600L, seed = 11L)
  d$m <- factor(sample(c("p", "q"), nrow(d), replace = TRUE))
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous",
      effect_modifier = "m"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("polytomous is rejected on a matched design", {
  d <- make_matched_cc(n_sets = 20L)
  d$g <- factor(sample(c("control", "caseA", "caseB"), nrow(d), replace = TRUE))
  expect_error(
    matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = matched_cc(strata = "set"),
      estimator = "polytomous"
    ),
    class = "matchatr_bad_estimator"
  )
})

# ---- Missing data -----------------------------------------------------------

test_that("missing values warn and the fit uses complete cases only", {
  d <- make_polytomous_cc(n = 1200L, seed = 13L)
  d$x[1:7] <- NA
  expect_warning(
    fit <- matcha(
      d,
      outcome = "g",
      exposure = "x",
      design = unmatched_cc(),
      estimator = "polytomous"
    ),
    class = "matchatr_dropped_rows"
  )
  res <- contrast(fit, type = "or")
  expect_equal(res$n, nrow(d) - 7L)
})

# ---- Dispatch ---------------------------------------------------------------

test_that("polytomous routes to the multinom engine with an OR default", {
  routing <- resolve_engine("unmatched_cc", "polytomous")
  expect_identical(routing$engine, "multinom")
  expect_identical(routing$outcome_kind, "polytomous")
  expect_false(routing$conditional)
  expect_identical(default_contrast_type("multinom"), "or")
})

test_that("the result vcov is the selected, symmetric exposure submatrix", {
  d <- make_polytomous_cc(n = 1500L, seed = 14L)
  fit <- matcha(
    d,
    outcome = "g",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous"
  )
  res <- contrast(fit, type = "or")
  V <- res$vcov
  expect_equal(dim(V), c(2L, 2L))
  expect_equal(rownames(V), c("caseA: x", "caseB: x"))
  expect_equal(V, t(V))
  expect_equal(
    unname(sqrt(diag(V))),
    res$estimates$se,
    tolerance = 1e-10
  )
})
