# contrast() is the second-step verb. In this phase no estimation engine has
# run, so it fixes the public signature and aborts with matchatr_not_estimated
# on any fit (model = NULL). These tests pin the verb's contract.

test_that("contrast rejects a non-fit argument", {
  expect_error(contrast(list()), class = "matchatr_bad_input")
  expect_error(contrast(unmatched_cc()), class = "matchatr_bad_input")
})

test_that("contrast on an unestimated fit aborts with matchatr_not_estimated", {
  df <- make_cc_data(n_sets = 6L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_error(contrast(fit), class = "matchatr_not_estimated")
  # Holds across designs / estimators.
  fit2 <- matcha(df, "case", "x", matched_cc(strata = "set"))
  expect_error(contrast(fit2), class = "matchatr_not_estimated")
})

test_that("contrast validates the contrast scale and CI method", {
  df <- make_cc_data(n_sets = 6L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  # match.arg rejects an out-of-set value before the estimation guard fires.
  expect_error(contrast(fit, type = "bogus"))
  expect_error(contrast(fit, ci_method = "bogus"))
})

test_that("the not-estimated message reads clearly", {
  df <- make_cc_data(n_sets = 6L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_snapshot(contrast(fit), error = TRUE)
})
