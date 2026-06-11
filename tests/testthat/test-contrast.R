# contrast() is the second-step verb. It validates its arguments, then either
# dispatches to the engine's contrast assembly (logistic -> conditional OR) or
# aborts with matchatr_not_estimated for an engine with no wired estimator.
# The logistic OR values / rejections live in test-unconditional.R; here we pin
# the verb's generic contract.

test_that("contrast rejects a non-fit argument", {
  expect_error(contrast(list()), class = "matchatr_bad_input")
  expect_error(contrast(unmatched_cc()), class = "matchatr_bad_input")
})

test_that("contrast on a non-estimated engine aborts with matchatr_not_estimated", {
  df <- make_cc_data(n_sets = 6L)
  # The CCW IPW / AIPW / TMLE engines are not yet wired; their fit carries
  # model = NULL (ccw_gformula is wired, exercised in test-ccw.R).
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(prevalence = 0.05),
    estimator = "ccw_ipw"
  )
  expect_null(fit$model)
  expect_error(contrast(fit), class = "matchatr_not_estimated")
})

test_that("contrast validates the scale, CI method, and confidence level", {
  df <- make_cc_data(n_sets = 6L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  # match.arg rejects out-of-set values.
  expect_error(contrast(fit, type = "bogus"))
  expect_error(contrast(fit, ci_method = "bogus"))
  # conf_level must be a single probability strictly in (0, 1).
  expect_error(
    contrast(fit, type = "or", conf_level = 0),
    class = "matchatr_bad_input"
  )
  expect_error(
    contrast(fit, type = "or", conf_level = 1),
    class = "matchatr_bad_input"
  )
  expect_error(
    contrast(fit, type = "or", conf_level = c(0.9, 0.95)),
    class = "matchatr_bad_input"
  )
})

test_that("the three ci_method choices all resolve for the logistic OR", {
  df <- make_cc_data(n_sets = 20L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_s3_class(
    contrast(fit, type = "or", ci_method = "model"),
    "matchatr_result"
  )
  expect_s3_class(
    contrast(fit, type = "or", ci_method = "sandwich"),
    "matchatr_result"
  )
  # bootstrap is the one CI method the logistic OR declines.
  expect_error(
    contrast(fit, type = "or", ci_method = "bootstrap"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("contrast() defaults to the estimand the engine identifies", {
  df <- make_cc_data(n_sets = 20L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  # No `type` given: the logistic engine identifies the OR, so the default
  # resolves to it rather than the risk difference it would have to reject.
  res <- contrast(fit)
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  # An explicitly requested unidentified estimand is still rejected.
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
})

test_that("the not-estimated message reads clearly", {
  df <- make_cc_data(n_sets = 6L)
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(prevalence = 0.05),
    estimator = "ccw_ipw"
  )
  expect_snapshot(contrast(fit), error = TRUE)
})
