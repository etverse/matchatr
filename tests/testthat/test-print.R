# Print methods are user-facing; snapshot their output so any wording / layout
# change is reviewed deliberately.

test_that("print.matchatr_design renders each design", {
  expect_snapshot(print(unmatched_cc(prevalence = 0.02)))
  expect_snapshot(print(matched_cc(strata = c("age_grp", "sex"), ratio = 2)))
  expect_snapshot(print(nested_cc(strata = "set", time = "t", ratio = 3)))
  expect_snapshot(print(case_cohort(subcohort = "in_subcohort", time = "t")))
  expect_snapshot(print(two_phase(phase1 = "stratum", phase2 = "in_phase2")))
  expect_snapshot(print(counter_matched(strata = "surrogate", time = "t")))
})

test_that("print.matchatr_fit renders the resolved analysis", {
  df <- make_cc_data(n_sets = 20L, ratio = 2L)

  fit_adj <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(),
    confounders = ~ age + smoke
  )
  expect_snapshot(print(fit_adj))

  fit_clogit <- matcha(df, "case", "x", matched_cc(strata = "set"))
  expect_snapshot(print(fit_clogit))

  fit_ccw <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(prevalence = 0.02),
    confounders = ~age,
    estimator = "ccw_gformula"
  )
  expect_snapshot(print(fit_ccw))
})

test_that("print methods return their input invisibly", {
  d <- unmatched_cc()
  expect_invisible(print(d))
  expect_identical(withVisible(print(d))$value, d)

  fit <- matcha(make_cc_data(), "case", "x", unmatched_cc())
  expect_invisible(print(fit))
})
