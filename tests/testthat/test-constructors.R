# The low-level fit / result constructors define the object shapes the rest of
# the package (and later phases) build on. These tests pin the class and slot
# layout so a future refactor cannot silently drop a slot.

test_that("new_matchatr_fit builds the documented slot layout", {
  design <- matched_cc(strata = "set")
  fit <- new_matchatr_fit(
    model = NULL,
    data = data.table::data.table(case = 0:1, x = c(1, 0)),
    outcome = "case",
    exposure = "x",
    confounders = ~age,
    design = design,
    estimator = "clogit",
    engine = "clogit",
    details = list(engine = "clogit", variance_kind = NULL),
    call = quote(matcha())
  )
  expect_s3_class(fit, "matchatr_fit")
  expect_named(
    fit,
    c(
      "model",
      "data",
      "outcome",
      "exposure",
      "confounders",
      "design",
      "estimator",
      "engine",
      "details",
      "call"
    )
  )
  expect_null(fit$model)
  expect_s3_class(fit$design, "matchatr_design")
  # The variance correction is reserved (NULL) until an inference engine runs.
  expect_true("variance_kind" %in% names(fit$details))
  expect_null(fit$details$variance_kind)
})

test_that("new_matchatr_result builds the documented slot layout", {
  res <- new_matchatr_result(
    estimates = data.table::data.table(level = "x", est = 1.2),
    contrasts = data.table::data.table(contrast = "x vs ref", est = 0.3),
    type = "difference",
    estimand = "RD",
    ci_method = "sandwich",
    reference = "ref",
    n = 100L,
    estimator = "ccw_gformula",
    engine = "ccw_gformula",
    vcov = matrix(1, 1, 1),
    call = quote(contrast())
  )
  expect_s3_class(res, "matchatr_result")
  expect_named(
    res,
    c(
      "estimates",
      "contrasts",
      "type",
      "estimand",
      "ci_method",
      "reference",
      "n",
      "estimator",
      "engine",
      "vcov",
      "call"
    )
  )
  expect_identical(res$ci_method, "sandwich")
  expect_identical(res$n, 100L)
})
