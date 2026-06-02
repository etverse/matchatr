# Design constructors build the right matchatr_design objects and reject
# malformed structural arguments with classed errors. No data is involved here
# (constructors never see data); column-existence checks live in matcha().

test_that("unmatched_cc carries prevalence and the right weight scheme", {
  d0 <- unmatched_cc()
  expect_s3_class(d0, "matchatr_design")
  expect_identical(d0$type, "unmatched_cc")
  expect_null(d0$prevalence)
  expect_identical(d0$weight_spec$kind, "none")

  d1 <- unmatched_cc(prevalence = 0.02)
  expect_identical(d1$prevalence, 0.02)
  expect_identical(d1$weight_spec$kind, "case_control")
  expect_identical(d1$weight_spec$prevalence, 0.02)
})

test_that("matched_cc stores strata (vector) and ratio", {
  d <- matched_cc(strata = "set")
  expect_identical(d$type, "matched_cc")
  expect_identical(d$strata, "set")
  expect_null(d$ratio)
  expect_identical(d$weight_spec$kind, "none")

  d2 <- matched_cc(strata = c("age_grp", "sex"), ratio = 2)
  expect_identical(d2$strata, c("age_grp", "sex"))
  expect_identical(d2$ratio, 2)
})

test_that("nested_cc requires strata + time and flags inclusion weights", {
  d <- nested_cc(strata = "set", time = "t", ratio = 3)
  expect_identical(d$type, "nested_cc")
  expect_identical(d$strata, "set")
  expect_identical(d$time, "t")
  expect_identical(d$ratio, 3)
  expect_identical(d$weight_spec$kind, "inclusion")
})

test_that("case_cohort stores subcohort + time", {
  d <- case_cohort(subcohort = "in_subcohort", time = "t")
  expect_identical(d$type, "case_cohort")
  expect_identical(d$subcohort, "in_subcohort")
  expect_identical(d$time, "t")
  expect_identical(d$weight_spec$kind, "inclusion")
})

test_that("two_phase stores phase-1 strata and phase-2 selection", {
  d <- two_phase(phase1 = c("stratum", "region"), phase2 = "in_phase2")
  expect_identical(d$type, "two_phase")
  expect_identical(d$phase1, c("stratum", "region"))
  expect_identical(d$phase2, "in_phase2")
  expect_identical(d$weight_spec$kind, "design")
})

test_that("counter_matched stores strata + time", {
  d <- counter_matched(strata = "exposure_surrogate", time = "t")
  expect_identical(d$type, "counter_matched")
  expect_identical(d$strata, "exposure_surrogate")
  expect_identical(d$weight_spec$kind, "counter_match")
})

# --- structural validators ----------------------------------------------

test_that("prevalence must be a single number strictly in (0, 1)", {
  expect_error(unmatched_cc(prevalence = 0), class = "matchatr_bad_prevalence")
  expect_error(unmatched_cc(prevalence = 1), class = "matchatr_bad_prevalence")
  expect_error(
    unmatched_cc(prevalence = -0.1),
    class = "matchatr_bad_prevalence"
  )
  expect_error(
    unmatched_cc(prevalence = 1.5),
    class = "matchatr_bad_prevalence"
  )
  expect_error(
    unmatched_cc(prevalence = "x"),
    class = "matchatr_bad_prevalence"
  )
  expect_error(
    unmatched_cc(prevalence = c(0.1, 0.2)),
    class = "matchatr_bad_prevalence"
  )
  expect_error(
    unmatched_cc(prevalence = NA_real_),
    class = "matchatr_bad_prevalence"
  )
})

test_that("ratio must be a single whole number >= 1", {
  expect_error(matched_cc("set", ratio = 0), class = "matchatr_bad_ratio")
  expect_error(matched_cc("set", ratio = -1), class = "matchatr_bad_ratio")
  expect_error(matched_cc("set", ratio = 1.5), class = "matchatr_bad_ratio")
  expect_error(matched_cc("set", ratio = "2"), class = "matchatr_bad_ratio")
  expect_error(matched_cc("set", ratio = c(1, 2)), class = "matchatr_bad_ratio")
  expect_error(
    matched_cc("set", ratio = NA_integer_),
    class = "matchatr_bad_ratio"
  )
  # An integer-typed whole number is fine.
  expect_no_error(matched_cc("set", ratio = 3L))
})

test_that("a non-finite ratio is rejected with the classed error, not a base crash", {
  # Review 2026-06-02 Issue B1: `Inf %% 1` is NaN, so the old guard hit
  # `if (NA)` and raised an unclassed base error. /tmp/matchatr_repro_ratio_inf.R
  expect_error(matched_cc("set", ratio = Inf), class = "matchatr_bad_ratio")
  expect_error(nested_cc("set", "t", ratio = -Inf), class = "matchatr_bad_ratio")
  expect_error(matched_cc("set", ratio = NaN), class = "matchatr_bad_ratio")
})

test_that("strata must be a non-empty character vector", {
  expect_error(matched_cc(strata = 1), class = "matchatr_bad_strata")
  expect_error(matched_cc(strata = character(0)), class = "matchatr_bad_strata")
  expect_error(matched_cc(strata = c("a", "")), class = "matchatr_bad_strata")
  expect_error(
    matched_cc(strata = NA_character_),
    class = "matchatr_bad_strata"
  )
  expect_error(nested_cc(strata = 5, time = "t"), class = "matchatr_bad_strata")
})

test_that("time and subcohort must be single non-empty strings", {
  expect_error(nested_cc("set", time = 5), class = "matchatr_bad_input")
  expect_error(
    nested_cc("set", time = c("a", "b")),
    class = "matchatr_bad_input"
  )
  expect_error(nested_cc("set", time = ""), class = "matchatr_bad_input")
  expect_error(
    case_cohort(subcohort = 1, time = "t"),
    class = "matchatr_bad_input"
  )
  expect_error(
    case_cohort(subcohort = "sc", time = NA),
    class = "matchatr_bad_input"
  )
})

test_that("every classed design error also carries the matchatr_error parent", {
  expect_error(unmatched_cc(prevalence = 2), class = "matchatr_error")
  expect_error(matched_cc("set", ratio = 0), class = "matchatr_error")
  expect_error(matched_cc(strata = 1), class = "matchatr_error")
})
