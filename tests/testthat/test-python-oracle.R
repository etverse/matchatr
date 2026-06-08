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
