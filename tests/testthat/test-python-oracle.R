# Cross-language oracles: each implemented estimator is cross-checked against an
# independent statsmodels (Python) fit on the SAME committed dataset, so a bug
# shared between matchatr and its R engine cannot hide behind a same-package
# comparison. The Python outputs are pre-computed fixtures under
# fixtures/python/ (see its README); these tests read the fixtures and never
# invoke Python, so CI needs no Python toolchain.

test_that("the unmatched logistic OR matches the statsmodels oracle", {
  res_path <- test_path("fixtures", "python", "logistic_or_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "logistic_or_data.csv"))
  py <- read.csv(res_path)
  # statsmodels labels the intercept "Intercept"; R's model.matrix uses
  # "(Intercept)". Align so the per-term lookup matches.
  py$term[py$term == "Intercept"] <- "(Intercept)"
  rownames(py) <- py$term

  fit <- matcha(
    data,
    "case",
    "x",
    unmatched_cc(),
    confounders = ~age,
    estimator = "logistic"
  )

  # Log-scale coefficients + SEs per term: statsmodels Logit and R's glm solve
  # the same likelihood, so they agree to optimiser tolerance. Read the tidy
  # table as a plain data.frame and index with vector subscripts -- inside
  # data.table's `[`, a bare name matching a column (here `trm` vs the `term`
  # column) is resolved to the column, not the loop variable.
  td <- as.data.frame(tidy(fit))
  for (trm in c("(Intercept)", "x", "age")) {
    expect_equal(
      td$estimate[td$term == trm],
      py[trm, "estimate"],
      tolerance = 1e-4
    )
    expect_equal(
      td$std.error[td$term == trm],
      py[trm, "std_error"],
      tolerance = 1e-4
    )
  }

  # Exposure odds ratio + 95% Wald interval (log scale, exponentiated) on the OR
  # scale -- the quantity contrast() reports.
  cx <- contrast(fit, type = "or")$contrasts
  expect_equal(cx$estimate, py["x", "odds_ratio"], tolerance = 1e-4)
  expect_equal(cx$ci_lower, py["x", "conf_low"], tolerance = 1e-4)
  expect_equal(cx$ci_upper, py["x", "conf_high"], tolerance = 1e-4)
})

test_that("the matched conditional OR matches the statsmodels ConditionalLogit oracle", {
  res_path <- test_path("fixtures", "python", "matched_or_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "matched_or_data.csv"))
  py <- read.csv(res_path)
  rownames(py) <- py$term

  fit <- matcha(
    data,
    "case",
    "x",
    matched_cc(strata = "set"),
    confounders = ~z,
    estimator = "clogit"
  )
  # Conditional log-ORs + SEs per term (no intercept). survival::clogit and
  # statsmodels' ConditionalLogit maximise the same conditional likelihood; two
  # independent partial-likelihood optimisers agree to ~1e-3 (a real bug -- a
  # wrong likelihood or stratification -- would differ by far more).
  tol <- 1e-3
  td <- as.data.frame(tidy(fit))
  for (trm in c("x", "z")) {
    expect_equal(
      td$estimate[td$term == trm],
      py[trm, "estimate"],
      tolerance = tol
    )
    expect_equal(
      td$std.error[td$term == trm],
      py[trm, "std_error"],
      tolerance = tol
    )
  }
  cx <- contrast(fit, type = "or")$contrasts
  expect_equal(cx$estimate, py["x", "odds_ratio"], tolerance = tol)
  expect_equal(cx$ci_lower, py["x", "conf_low"], tolerance = tol)
  expect_equal(cx$ci_upper, py["x", "conf_high"], tolerance = tol)
})

test_that("the nested case-control HR matches the statsmodels ConditionalLogit oracle", {
  res_path <- test_path("fixtures", "python", "nested_hr_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "nested_hr_data.csv"))
  py <- read.csv(res_path)
  rownames(py) <- py$term

  fit <- matcha(
    data,
    "case",
    "x",
    nested_cc(strata = "set", time = "risk_time"),
    confounders = ~z,
    estimator = "clogit"
  )
  # The risk-set conditional partial likelihood; exp(beta) IS the hazard ratio
  # under incidence-density sampling, so the same ConditionalLogit oracle applies.
  # ~1e-3 agreement of two independent partial-likelihood optimisers.
  tol <- 1e-3
  td <- as.data.frame(tidy(fit))
  for (trm in c("x", "z")) {
    expect_equal(
      td$estimate[td$term == trm],
      py[trm, "estimate"],
      tolerance = tol
    )
    expect_equal(
      td$std.error[td$term == trm],
      py[trm, "std_error"],
      tolerance = tol
    )
  }
  cx <- contrast(fit, type = "hr")$contrasts
  expect_equal(cx$estimate, py["x", "hazard_ratio"], tolerance = tol)
  expect_equal(cx$ci_lower, py["x", "conf_low"], tolerance = tol)
  expect_equal(cx$ci_upper, py["x", "conf_high"], tolerance = tol)
})

test_that("the Mantel-Haenszel OR matches the statsmodels StratifiedTable oracle", {
  res_path <- test_path("fixtures", "python", "mh_or_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "mh_or_data.csv"))
  py <- read.csv(res_path)

  fit <- matcha(
    data,
    "case",
    "x",
    unmatched_cc(strata = "agegrp"),
    estimator = "mh"
  )
  # Pooled OR + Robins-Breslow-Greenland interval, the same estimator and
  # variance statsmodels' StratifiedTable computes.
  cx <- contrast(fit, type = "or")$contrasts
  expect_equal(cx$estimate, py$odds_ratio[1], tolerance = 1e-4)
  expect_equal(cx$ci_lower, py$conf_low[1], tolerance = 1e-4)
  expect_equal(cx$ci_upper, py$conf_high[1], tolerance = 1e-4)
})

test_that("the polytomous subtype ORs match the statsmodels MNLogit oracle", {
  res_path <- test_path("fixtures", "python", "polytomous_or_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "polytomous_or_data.csv"))
  py <- read.csv(res_path)

  fit <- matcha(
    data,
    "g",
    "x",
    unmatched_cc(),
    confounders = ~age,
    estimator = "polytomous",
    reference = "control"
  )
  # Each non-reference subtype's exposure OR vs the control reference. nnet's
  # multinom and statsmodels' MNLogit maximise the same baseline-category
  # multinomial likelihood (the confounder `age` MUST match the Python model).
  cx <- as.data.frame(contrast(fit, type = "or")$contrasts)
  for (lev in c("caseA", "caseB")) {
    pr <- py[py$y_level == lev & py$term == "x", ]
    crow <- cx[cx$comparison == paste0(lev, ": x"), ]
    expect_equal(crow$estimate, pr$odds_ratio, tolerance = 1e-4)
    expect_equal(crow$ci_lower, pr$conf_low, tolerance = 1e-4)
    expect_equal(crow$ci_upper, pr$conf_high, tolerance = 1e-4)
  }
})

test_that("the homogeneity Wald test + pooled OR match the statsmodels oracle", {
  res_path <- test_path("fixtures", "python", "homogeneity_results.csv")
  skip_if(!file.exists(res_path), "Python oracle fixture not generated")
  data <- read.csv(test_path("fixtures", "python", "homogeneity_data.csv"))
  py <- read.csv(res_path)

  fit <- matcha(
    data,
    "g",
    "x",
    unmatched_cc(),
    estimator = "polytomous",
    reference = "control"
  )
  hr <- as.data.frame(tidy(test_homogeneity(fit)))
  hr <- hr[hr$term == "x", ]
  # matchatr reconstructs the Wald homogeneity chi-square and the GLS-pooled
  # common OR from nnet's multinomial fit; the Python oracle does the same from
  # statsmodels' MNLogit. Two independent multinomial optimisers agree to ~1e-3.
  tol <- 1e-3
  expect_identical(as.integer(hr$df), as.integer(py$df[1]))
  expect_equal(hr$statistic, py$chisq[1], tolerance = tol)
  expect_equal(hr$estimate, py$pooled_or[1], tolerance = tol)
  expect_equal(hr$conf.low, py$pooled_or_low[1], tolerance = tol)
  expect_equal(hr$conf.high, py$pooled_or_high[1], tolerance = tol)
})
