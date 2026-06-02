# Snapshot the full text of every classed rejection so the wording -- which is
# the user's only guide out of a misuse -- is reviewed and locked. Class-based
# assertions live in the per-component test files; these pin the messages.

test_that("design-constructor rejections read clearly", {
  expect_snapshot(unmatched_cc(prevalence = 1.5), error = TRUE)
  expect_snapshot(matched_cc(strata = "set", ratio = 1.5), error = TRUE)
  expect_snapshot(matched_cc(strata = character(0)), error = TRUE)
  expect_snapshot(nested_cc(strata = "set", time = 5), error = TRUE)
})

test_that("matcha rejections read clearly", {
  df <- make_cc_data(n_sets = 5L)

  expect_snapshot(
    matcha(df, "case", "x", unmatched_cc(), estimator = "bogus"),
    error = TRUE
  )
  expect_snapshot(
    matcha(df, "case", "x", unmatched_cc(), estimator = "ccw_ipw"),
    error = TRUE
  )
  expect_snapshot(
    matcha(df, "case", "missing_col", unmatched_cc()),
    error = TRUE
  )
  expect_snapshot(
    matcha(df, "case", "x", matched_cc(strata = "no_such_set")),
    error = TRUE
  )
  expect_snapshot(
    matcha(df, "case", "case", unmatched_cc()),
    error = TRUE
  )

  bad_y <- df
  bad_y$case <- rnorm(nrow(df))
  expect_snapshot(
    matcha(bad_y, "case", "x", unmatched_cc()),
    error = TRUE
  )

  expect_snapshot(
    matcha(df, "case", "x", "not a design"),
    error = TRUE
  )
})

test_that("the uninformative-stratum warning reads clearly", {
  bad <- make_uninformative_cc()
  expect_snapshot(
    invisible(matcha(bad, "case", "x", matched_cc(strata = "set")))
  )
})
