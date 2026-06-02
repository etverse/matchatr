# matcha() is the fit verb: validate -> build/accept design -> resolve dispatch
# -> return a matchatr_fit. No estimation happens in this phase, so model is
# NULL; the tests pin the resolved routing, the validated spec, and every
# data-level rejection path.

test_that("matcha returns a matchatr_fit with the resolved spec", {
  df <- make_cc_data()
  fit <- matcha(
    df,
    outcome = "case",
    exposure = "x",
    design = unmatched_cc(),
    confounders = ~ age + smoke
  )
  expect_s3_class(fit, "matchatr_fit")
  expect_null(fit$model)
  expect_s3_class(fit$data, "data.table")
  expect_identical(fit$outcome, "case")
  expect_identical(fit$exposure, "x")
  expect_identical(fit$estimator, "logistic")
  expect_identical(fit$engine, "glm_logistic")
  expect_identical(fit$details$kind, "classical")
})

test_that("matcha does not mutate the caller's data", {
  df <- make_cc_data()
  snapshot <- df
  fit <- matcha(df, "case", "x", matched_cc(strata = "set"))
  # Original frame unchanged (matcha works on a data.table copy).
  expect_identical(df, snapshot)
  expect_false(data.table::is.data.table(df))
  expect_true(data.table::is.data.table(fit$data))
})

test_that("estimator defaults to the design's canonical analysis", {
  df <- make_cc_data()
  expect_identical(
    matcha(df, "case", "x", unmatched_cc())$estimator,
    "logistic"
  )
  expect_identical(
    matcha(df, "case", "x", matched_cc("set"))$estimator,
    "clogit"
  )
  expect_identical(
    matcha(df, "case", "x", nested_cc("set", "t"))$estimator,
    "clogit"
  )
})

test_that("case / control counts are recorded on the fit", {
  df <- make_cc_data(n_sets = 20L, ratio = 2L)
  fit <- matcha(df, "case", "x", unmatched_cc())
  expect_identical(fit$details$n_cases, 20L)
  expect_identical(fit$details$n_controls, 40L)
})

# --- binary-outcome resolution ------------------------------------------

test_that("binary outcomes are accepted in their three sane encodings", {
  base <- make_cc_data(n_sets = 10L)

  # numeric 0/1
  expect_no_error(matcha(base, "case", "x", unmatched_cc()))

  # logical
  dlog <- base
  dlog$case <- as.logical(dlog$case)
  fit_log <- matcha(dlog, "case", "x", unmatched_cc())
  expect_identical(fit_log$details$n_cases, 10L)

  # two-level factor; second level is the case (glm/clogit event convention)
  dfac <- base
  dfac$case <- factor(
    ifelse(base$case == 1L, "case", "control"),
    levels = c("control", "case")
  )
  fit_fac <- matcha(dfac, "case", "x", unmatched_cc())
  expect_identical(fit_fac$details$n_cases, 10L)
})

test_that("non-binary outcomes are rejected with matchatr_bad_outcome", {
  base <- make_cc_data(n_sets = 10L)

  cont <- base
  cont$case <- rnorm(nrow(cont))
  expect_error(
    matcha(cont, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )

  three <- base
  three$case <- rep(c(0L, 1L, 2L), length.out = nrow(three))
  expect_error(
    matcha(three, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )

  threefac <- base
  threefac$case <- factor(rep(c("a", "b", "c"), length.out = nrow(threefac)))
  expect_error(
    matcha(threefac, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )

  # Degenerate: all cases -- no contrast.
  allcase <- base
  allcase$case <- 1L
  expect_error(
    matcha(allcase, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )
})

test_that("degenerate single-class outcomes are rejected for every encoding", {
  # Review 2026-06-02 Issue R2: the contrast check was numeric-only, so an
  # all-controls logical or 2-level factor slipped through with n_cases = 0.
  # /tmp/matchatr_repro_degenerate_outcome.R
  base <- make_cc_data(n_sets = 6L)

  alllog <- base
  alllog$case <- rep(FALSE, nrow(base))
  expect_error(
    matcha(alllog, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )

  # 2-level factor where only the "control" level actually occurs.
  allfac <- base
  allfac$case <- factor(rep("control", nrow(base)), levels = c("control", "case"))
  expect_error(
    matcha(allfac, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )

  # all-NA outcome: no contrast either.
  allna <- base
  allna$case <- NA_integer_
  expect_error(
    matcha(allna, "case", "x", unmatched_cc()),
    class = "matchatr_bad_outcome"
  )
})

# --- column existence + role checks -------------------------------------

test_that("missing columns are rejected with matchatr_bad_design", {
  df <- make_cc_data()
  expect_error(
    matcha(df, "nope", "x", unmatched_cc()),
    class = "matchatr_bad_design"
  )
  expect_error(
    matcha(df, "case", "nope", unmatched_cc()),
    class = "matchatr_bad_design"
  )
  expect_error(
    matcha(df, "case", "x", unmatched_cc(), confounders = ~ age + nope),
    class = "matchatr_bad_design"
  )
  # design-referenced columns are checked too
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "no_set")),
    class = "matchatr_bad_design"
  )
  expect_error(
    matcha(df, "case", "x", nested_cc(strata = "set", time = "no_t")),
    class = "matchatr_bad_design"
  )
})

test_that("structural argument misuse is rejected with classed errors", {
  df <- make_cc_data()
  expect_error(
    matcha(list(a = 1), "case", "x", unmatched_cc()),
    class = "matchatr_bad_input"
  )
  expect_error(
    matcha(df, "case", "x", "not a design"),
    class = "matchatr_bad_design"
  )
  expect_error(matcha(df, 1, "x", unmatched_cc()), class = "matchatr_bad_input")
  expect_error(
    matcha(df, "case", 1, unmatched_cc()),
    class = "matchatr_bad_input"
  )
  expect_error(
    matcha(df, "case", "x", unmatched_cc(), confounders = "age"),
    class = "matchatr_bad_input"
  )
  expect_error(
    matcha(df, "case", "case", unmatched_cc()),
    class = "matchatr_bad_input"
  )
})

# --- CCW prevalence requirement -----------------------------------------

test_that("CCW estimators require a prevalence on the design", {
  df <- make_cc_data()
  for (e in ccw_estimators()) {
    expect_error(
      matcha(df, "case", "x", unmatched_cc(), estimator = e),
      class = "matchatr_missing_prevalence"
    )
  }
  # With q0 supplied, the same estimators resolve.
  fit <- matcha(
    df,
    "case",
    "x",
    unmatched_cc(prevalence = 0.02),
    confounders = ~age,
    estimator = "ccw_gformula"
  )
  expect_identical(fit$details$kind, "ccw")
  expect_identical(fit$engine, "ccw_gformula")
})

test_that("dispatch accepts CCW on a matched design but still demands q0", {
  # The matched/nested constructors do not expose `prevalence`, so CCW on them
  # always trips the q0 requirement in this layer -- the dispatch accepts the
  # estimator (no bad_estimator), the prevalence gate rejects it.
  df <- make_cc_data()
  expect_error(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "ccw_aipw"),
    class = "matchatr_missing_prevalence"
  )
})

# --- uninformative-stratum warning --------------------------------------

test_that("conditional likelihood warns on uninformative strata", {
  bad <- make_uninformative_cc()
  expect_warning(
    matcha(bad, "case", "x", matched_cc(strata = "set"), estimator = "clogit"),
    class = "matchatr_uninformative_stratum"
  )
})

test_that("a fully informative matched sample raises no stratum warning", {
  df <- make_cc_data()
  expect_no_warning(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "clogit")
  )
})

test_that("a non-conditional analysis never emits the stratum warning", {
  # Same uninformative shape, but an unmatched_cc + logistic analysis is not a
  # conditional likelihood, so the warning gate stays shut.
  bad <- make_uninformative_cc()
  expect_no_warning(matcha(bad, "case", "x", unmatched_cc()))
})

test_that("warn_uninformative_strata flags only sets missing a case or control", {
  # set 1: one case + one control (informative); set 2: two cases (no control).
  strata <- list(c(1L, 1L, 2L, 2L))
  y <- c(1L, 0L, 1L, 1L)
  expect_warning(
    warn_uninformative_strata(strata, y),
    class = "matchatr_uninformative_stratum"
  )

  # Every set has both a case and a control: silent.
  strata_ok <- list(c(1L, 1L, 2L, 2L))
  y_ok <- c(1L, 0L, 1L, 0L)
  expect_no_warning(warn_uninformative_strata(strata_ok, y_ok))

  # Multi-column strata: combine via interaction; (a1,b1) has only a case.
  strata_mc <- list(c("a", "a", "b"), c(1L, 1L, 2L))
  y_mc <- c(1L, 0L, 1L)
  expect_warning(
    warn_uninformative_strata(strata_mc, y_mc),
    class = "matchatr_uninformative_stratum"
  )
})
