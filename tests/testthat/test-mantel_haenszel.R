# The Mantel-Haenszel stratified odds-ratio engine. Oracle: stats::mantelhaen.test
# (correct = FALSE), whose OR and odds-ratio confidence interval use the same
# Robins-Breslow-Greenland variance matchatr computes; plus the closed-form 2x2
# odds ratio for the crude (single-stratum) case.

# Build the [exposure, outcome, stratum] 2x2xK array mantelhaen.test expects.
mh_oracle <- function(df, exposure = "x", strata = "agegrp") {
  arr <- table(
    factor(df[[exposure]], c(1, 0)),
    factor(df$case, c(1, 0)),
    df[[strata]]
  )
  stats::mantelhaen.test(arr, correct = FALSE)
}

test_that("the stratified MH OR and CI match mantelhaen.test (RBG variance)", {
  df <- make_stratified_cc()
  res <- contrast(
    matcha(df, "case", "x", unmatched_cc(strata = "agegrp"), estimator = "mh"),
    type = "or"
  )
  oracle <- mh_oracle(df)

  expect_equal(
    res$contrasts$estimate,
    unname(oracle$estimate),
    tolerance = 1e-6
  )
  # The CI bounds use the Robins-Breslow-Greenland variance (so does the oracle).
  expect_equal(
    c(res$contrasts$ci_lower, res$contrasts$ci_upper),
    as.numeric(oracle$conf.int),
    tolerance = 1e-5
  )
})

test_that("the crude (no-strata) MH equals the closed-form 2x2 odds ratio", {
  df <- make_stratified_cc()
  res <- contrast(
    matcha(df, "case", "x", unmatched_cc(), estimator = "mh"),
    type = "or"
  )
  a <- sum(df$x == 1 & df$case == 1)
  b <- sum(df$x == 0 & df$case == 1)
  cc <- sum(df$x == 1 & df$case == 0)
  d <- sum(df$x == 0 & df$case == 0)
  expect_equal(res$contrasts$estimate, (a * d) / (b * cc), tolerance = 1e-8)
})

test_that("multi-column strata cross into a single stratifying factor", {
  df <- make_stratified_cc()
  res <- contrast(
    matcha(
      df,
      "case",
      "x",
      unmatched_cc(strata = c("agegrp", "sex")),
      estimator = "mh"
    ),
    type = "or"
  )
  # Oracle: mantelhaen.test on the agegrp:sex interaction strata.
  arr <- table(
    factor(df$x, c(1, 0)),
    factor(df$case, c(1, 0)),
    interaction(df$agegrp, df$sex, drop = TRUE)
  )
  oracle <- stats::mantelhaen.test(arr, correct = FALSE)
  expect_equal(
    res$contrasts$estimate,
    unname(oracle$estimate),
    tolerance = 1e-6
  )
})

test_that("a two-level factor exposure gives the same MH OR as 0/1", {
  df <- make_stratified_cc()
  or_num <- contrast(
    matcha(df, "case", "x", unmatched_cc(strata = "agegrp"), estimator = "mh"),
    type = "or"
  )$contrasts$estimate

  dff <- df
  dff$x <- factor(
    ifelse(df$x == 1L, "exposed", "unexposed"),
    levels = c("unexposed", "exposed")
  )
  or_fac <- contrast(
    matcha(dff, "case", "x", unmatched_cc(strata = "agegrp"), estimator = "mh"),
    type = "or"
  )$contrasts$estimate
  expect_equal(or_fac, or_num, tolerance = 1e-10)
})

test_that("the MH result is labelled and sized correctly", {
  df <- make_stratified_cc()
  res <- contrast(
    matcha(df, "case", "x", unmatched_cc(strata = "agegrp"), estimator = "mh"),
    type = "or"
  )
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "Mantel-Haenszel OR")
  expect_identical(res$contrasts$comparison, "x")
  expect_identical(res$n, nrow(df))
  # contrast() with no type defaults to the OR for the MH engine.
  expect_identical(
    contrast(matcha(
      df,
      "case",
      "x",
      unmatched_cc(strata = "agegrp"),
      estimator = "mh"
    ))$type,
    "or"
  )
})

# --- rejections ---------------------------------------------------------

test_that("RD / RR are rejected as unidentified for MH", {
  df <- make_stratified_cc()
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(strata = "agegrp"),
    estimator = "mh"
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

test_that("sandwich / bootstrap CIs are not available for MH (RBG only)", {
  df <- make_stratified_cc()
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(strata = "agegrp"),
    estimator = "mh"
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

test_that("a non-binary exposure is rejected by the MH estimator", {
  df <- make_stratified_cc()
  df$x <- factor(sample(c("lo", "mid", "hi"), nrow(df), replace = TRUE))
  expect_error(
    contrast(matcha(
      df,
      "case",
      "x",
      unmatched_cc(strata = "agegrp"),
      estimator = "mh"
    )),
    class = "matchatr_bad_input"
  )
})

test_that("a zero exposure-outcome margin is not estimable", {
  df <- make_stratified_cc()
  df$x <- 1L # everyone exposed -> no unexposed cell anywhere
  expect_error(
    contrast(matcha(
      df,
      "case",
      "x",
      unmatched_cc(strata = "agegrp"),
      estimator = "mh"
    )),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("MH rejection messages read clearly", {
  df <- make_stratified_cc()
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(strata = "agegrp"),
    estimator = "mh"
  )
  expect_snapshot(
    contrast(fit, type = "or", ci_method = "sandwich"),
    error = TRUE
  )
})
