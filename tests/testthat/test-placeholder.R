# Placeholder test so the suite is non-empty during the scaffold stage.
# Replaced by real design/estimator tests as the PHASE_*.md docs are implemented
# (starting with PHASE_1 design-object and dispatch tests).

test_that("package can be loaded", {
  expect_true(requireNamespace("matchatr", quietly = TRUE))
})
