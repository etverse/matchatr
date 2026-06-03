# The matched case-control conditional logistic engine: matcha() fits
# survival::clogit and contrast(type = "or") / tidy() report the conditional OR.
# Oracles: survival::clogit itself (exact pass-through), the handbook §4.4
# induced-abortion ORs on `infert`, and a matched-set DGP with a known
# conditional log-OR built from the conditional likelihood (Breslow & Day 1980).

# --- exact pass-through oracle on infert ---------------------------------

test_that("the clogit engine reproduces survival::clogit exactly on infert", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ induced + spontaneous + strata(stratum),
    data = infert
  )
  # The wrapper must build the same conditional likelihood as a hand-fit clogit.
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  expect_equal(stats::vcov(fit$model), stats::vcov(oracle))

  res <- contrast(fit, type = "or")
  ind <- res$contrasts[res$contrasts$comparison == "induced", ]
  b <- stats::coef(oracle)["induced"]
  se <- sqrt(stats::vcov(oracle)["induced", "induced"])
  z <- stats::qnorm(0.975)
  expect_equal(ind$estimate, unname(exp(b)), tolerance = 1e-8)
  # Wald interval on the log scale, exponentiated.
  expect_equal(
    c(ind$ci_lower, ind$ci_upper),
    unname(exp(b + c(-1, 1) * z * se)),
    tolerance = 1e-8
  )
})

test_that("the conditional ORs match the handbook induced-abortion values", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  td <- tidy(fit, exponentiate = TRUE)
  # Handbook §4.4 CMLE matched-set ORs: induced ~ 4.09 per prior abortion,
  # spontaneous ~ 7.29 (the conditional logistic fit on the matched strata).
  expect_equal(td$estimate[td$term == "induced"], 4.0919, tolerance = 1e-3)
  expect_equal(td$estimate[td$term == "spontaneous"], 7.2854, tolerance = 1e-3)
  # induced enters as a numeric trend, so OR(2+) = OR(1+)^2 ~ 16.7 (the handbook
  # two-or-more value), pinning the per-unit log-linear interpretation.
  b_ind <- log(td$estimate[td$term == "induced"])
  expect_equal(exp(2 * b_ind), 16.74, tolerance = 0.05)
})

# --- truth-based: CMLE recovers the conditional log-OR -------------------

test_that("the conditional OR recovers the matched-set log-OR (truth-based)", {
  df <- make_matched_cc(n_sets = 400L, ratio = 3L, beta_x = log(2.5))
  truth <- attr(df, "truth")
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  # The set-level exposure prevalence is conditioned away, so the CMLE recovers
  # exp(beta_x) up to sampling error (many sets, fixed seed).
  expect_equal(
    log(res$contrasts$estimate),
    unname(truth["beta_x"]),
    tolerance = 0.12
  )
})

test_that("adjusting for a non-matching covariate matches the clogit oracle", {
  df <- make_matched_cc(n_sets = 300L)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    confounders = ~z,
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ x + z + strata(set), data = df)
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  expect_equal(stats::vcov(fit$model), stats::vcov(oracle))
})

# --- factor exposure (per-level OR vs the reference) --------------------

test_that("a factor exposure reports the per-level OR and the reference", {
  df <- infert
  df$ind <- factor(
    ifelse(df$induced >= 1L, "yes", "no"),
    levels = c("no", "yes")
  )
  fit <- matcha(
    df,
    "case",
    "ind",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ ind + strata(stratum), data = df)
  res <- contrast(fit, type = "or")
  expect_identical(res$contrasts$comparison, "indyes")
  expect_identical(res$reference, "no")
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)["indyes"])),
    tolerance = 1e-8
  )
})

# --- frequency matching: several strata columns cross into one factor ----

test_that("multi-column strata cross into a single conditioning factor", {
  df <- make_matched_cc(n_sets = 250L, ratio = 3L)
  # Split the set id into two columns whose crossing uniquely recovers the set
  # (frequency matching on two variables).
  df$set_a <- (df$set - 1L) %/% 5L
  df$set_b <- (df$set - 1L) %% 5L
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = c("set_a", "set_b")),
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ x + strata(set_a, set_b),
    data = df
  )
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
})

# --- result structure / labels ------------------------------------------

test_that("the clogit result is labelled and sized correctly", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "conditional OR")
  expect_identical(res$estimator, "clogit")
  expect_identical(res$engine, "clogit")
  # Analysis n is the rows clogit used (coxph nobs() counts events, not rows).
  expect_identical(res$n, fit$model$n)
  # contrast() with no type defaults to the OR for the clogit engine.
  expect_identical(contrast(fit)$type, "or")
})

test_that("tidy reports the clogit ORs matching the oracle", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ induced + spontaneous + strata(stratum),
    data = infert
  )
  td <- tidy(fit, exponentiate = TRUE)
  # No intercept row (the conditional likelihood has none), one row per term.
  expect_setequal(td$term, c("induced", "spontaneous"))
  expect_equal(
    td$estimate,
    unname(exp(stats::coef(oracle))[td$term]),
    tolerance = 1e-8
  )
})

# --- rejections ----------------------------------------------------------

test_that("RD / RR are rejected as unidentified for clogit", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
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

test_that("sandwich / bootstrap CIs are not available for clogit", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
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

test_that("an exposure with no within-set variation is not estimable", {
  df <- make_matched_cc(n_sets = 50L)
  # Make x constant within every set: it has no conditional contribution, so
  # clogit aliases its coefficient to NA and no OR is identified.
  df$x <- as.integer(df$set %% 2L)
  fit <- suppressWarnings(matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  ))
  expect_error(
    contrast(fit, type = "or"),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("a missing value triggers the dropped-rows warning", {
  df <- infert
  df$induced[1:4] <- NA
  expect_warning(
    matcha(
      df,
      "case",
      "induced",
      matched_cc(strata = "stratum"),
      estimator = "clogit"
    ),
    class = "matchatr_dropped_rows"
  )
})

test_that("clogit rejection messages read clearly", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  expect_snapshot(
    contrast(fit, type = "or", ci_method = "sandwich"),
    error = TRUE
  )
})
