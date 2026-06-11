# cc_weights() builds the Rose & van der Laan case-control weights that map a
# case-control sample back to the source population. The defining properties are
# closed-form: the weights sum to n and the weighted case fraction equals q0
# exactly. These are the oracle.

test_that("cc_weights reproduce the source prevalence q0 in the weighted sample", {
  y <- c(rep(1L, 30L), rep(0L, 70L))
  q0 <- 0.05
  w <- cc_weights(q0, y)

  n <- length(y)
  # Weights sum to n (a probability-rescaling property, not a free normalization).
  expect_equal(sum(w), n)
  # The weighted case fraction is exactly q0.
  expect_equal(sum(w[y == 1L]) / sum(w), q0)
  # Each role's weight is q0 / (sample case fraction) and (1 - q0) / (control
  # fraction).
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  expect_equal(unique(w[y == 1L]), q0 / (n1 / n))
  expect_equal(unique(w[y == 0L]), (1 - q0) / (n0 / n))
})

test_that("cc_weights records the prevalence-known flag and defaults to TRUE", {
  y <- c(1L, 1L, 0L, 0L, 0L)
  expect_true(attr(cc_weights(0.1, y), "prevalence_known"))
  expect_false(attr(
    cc_weights(0.1, y, prevalence_known = FALSE),
    "prevalence_known"
  ))
})

test_that("cc_weights carry NA outcomes through as NA and exclude them from counts", {
  y <- c(1L, 1L, 0L, 0L, NA, NA)
  q0 <- 0.2
  w <- cc_weights(q0, y)

  expect_true(all(is.na(w[is.na(y)])))
  # n1 = 2, n0 = 2, n = 4 (the two NAs are excluded), so the non-missing weights
  # still rescale to the q0 margin over the complete cases.
  expect_equal(sum(w[!is.na(y) & y == 1L]) / sum(w, na.rm = TRUE), q0)
  expect_equal(w[1], q0 / (2 / 4))
  expect_equal(w[3], (1 - q0) / (2 / 4))
})

test_that("cc_weights reject a degenerate (single-class) sample", {
  expect_error(
    cc_weights(0.1, rep(1L, 10L)),
    class = "matchatr_bad_outcome"
  )
  expect_error(
    cc_weights(0.1, rep(0L, 10L)),
    class = "matchatr_bad_outcome"
  )
})
