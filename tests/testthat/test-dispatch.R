# The dispatch layer maps (design type, estimator) -> engine. These tests pin
# every routed pair and confirm the two orthogonal axes behave: a classical
# estimator is design-specific, while the case-control-weighted family applies
# to any design.

test_that("each design's classical estimator routes to the right engine", {
  expect_identical(
    resolve_engine("unmatched_cc", "logistic")$engine,
    "glm_logistic"
  )
  expect_identical(
    resolve_engine("unmatched_cc", "mh")$engine,
    "mantel_haenszel"
  )
  expect_identical(resolve_engine("matched_cc", "clogit")$engine, "clogit")
  expect_identical(resolve_engine("nested_cc", "clogit")$engine, "clogit")
  expect_identical(resolve_engine("case_cohort", "cch")$engine, "cch")
  expect_identical(
    resolve_engine("two_phase", "survey")$engine,
    "survey_twophase"
  )
  expect_identical(
    resolve_engine("counter_matched", "weighted_cox")$engine,
    "weighted_cox"
  )
})

test_that("the conditional flag is set only for clogit", {
  expect_true(resolve_engine("matched_cc", "clogit")$conditional)
  expect_true(resolve_engine("nested_cc", "clogit")$conditional)
  expect_false(resolve_engine("unmatched_cc", "logistic")$conditional)
  expect_false(resolve_engine("case_cohort", "cch")$conditional)
})

test_that("classical routes are tagged kind = 'classical'", {
  expect_identical(resolve_engine("unmatched_cc", "logistic")$kind, "classical")
  expect_identical(resolve_engine("matched_cc", "clogit")$kind, "classical")
})

test_that("CCW estimators route on any design with kind = 'ccw'", {
  designs <- c(
    "unmatched_cc",
    "matched_cc",
    "nested_cc",
    "case_cohort",
    "two_phase",
    "counter_matched"
  )
  for (d in designs) {
    for (e in ccw_estimators()) {
      r <- resolve_engine(d, e)
      expect_identical(r$kind, "ccw")
      expect_identical(r$engine, e)
      expect_false(r$conditional)
    }
  }
})

test_that("unknown estimators are rejected with matchatr_bad_estimator", {
  expect_error(
    resolve_engine("unmatched_cc", "bogus"),
    class = "matchatr_bad_estimator"
  )
  # An estimator that is valid for a *different* design is still rejected here.
  expect_error(
    resolve_engine("unmatched_cc", "clogit"),
    class = "matchatr_bad_estimator"
  )
  expect_error(
    resolve_engine("matched_cc", "logistic"),
    class = "matchatr_bad_estimator"
  )
  expect_error(
    resolve_engine("case_cohort", "clogit"),
    class = "matchatr_bad_estimator"
  )
})

test_that("default_estimator gives each design its canonical analysis", {
  expect_identical(default_estimator("unmatched_cc"), "logistic")
  expect_identical(default_estimator("matched_cc"), "clogit")
  expect_identical(default_estimator("nested_cc"), "clogit")
  expect_identical(default_estimator("case_cohort"), "cch")
  expect_identical(default_estimator("two_phase"), "survey")
  expect_identical(default_estimator("counter_matched"), "weighted_cox")
})

test_that("design_columns collects every referenced column and dedups", {
  expect_setequal(design_columns(unmatched_cc()), character(0))
  expect_setequal(design_columns(matched_cc(c("a", "b"))), c("a", "b"))
  expect_setequal(design_columns(nested_cc("set", "t")), c("set", "t"))
  expect_setequal(
    design_columns(case_cohort("sc", "t")),
    c("sc", "t")
  )
  expect_setequal(
    design_columns(two_phase(c("s1", "s2"), "p2")),
    c("s1", "s2", "p2")
  )
  # A column reused across slots appears once.
  d <- new_matchatr_design(type = "nested_cc", strata = "g", time = "g")
  expect_identical(design_columns(d), "g")
})
