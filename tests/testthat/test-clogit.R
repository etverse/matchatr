# The matched case-control conditional logistic engine: matcha() fits
# survival::clogit and contrast(type = "or") / tidy() report the conditional OR.
# Oracles: survival::clogit (exact pass-through), the 1:1 McNemar closed form
# (OR = n10/n01, Var(log OR) = 1/n10 + 1/n01 -- independent of clogit), and a
# matched-set DGP with a known conditional log-OR built from the conditional
# likelihood (Breslow & Day 1980). The infert ORs are a regression pin against
# the canonical survival::clogit example, not an external oracle.

# --- exact pass-through oracle on infert ---------------------------------

test_that("the clogit engine reproduces survival::clogit exactly on infert", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ induced + spontaneous + strata(stratum),
    data = infert
  )
  # The wrapper must build the same conditional likelihood as a hand-fit clogit.
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  expect_equal(stats::vcov(fit$model), stats::vcov(oracle))

  res <- contrast(fit, type = "or")
  ind <- res$contrasts[res$contrasts$comparison == "induced", ]
  b <- stats::coef(oracle)["induced"]
  se <- sqrt(stats::vcov(oracle)["induced", "induced"])
  z <- stats::qnorm(0.975)
  expect_equal(ind$estimate, unname(exp(b)), tolerance = 1e-8)
  # Wald interval on the log scale, exponentiated.
  expect_equal(
    c(ind$ci_lower, ind$ci_upper),
    unname(exp(b + c(-1, 1) * z * se)),
    tolerance = 1e-8
  )
})

test_that("the conditional ORs match the canonical infert clogit example", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  td <- tidy(fit, exponentiate = TRUE)
  # Regression pin against the canonical infert conditional-logistic fit (the
  # `survival::clogit` documented example, `?clogit`): induced ~ 4.09 per prior
  # induced abortion, spontaneous ~ 7.29. These are clogit's own output, so this
  # guards against a gross wrapper break rather than serving as an external
  # oracle (the independent statistical checks are the McNemar and truth tests).
  expect_equal(td$estimate[td$term == "induced"], 4.0919, tolerance = 1e-3)
  expect_equal(td$estimate[td$term == "spontaneous"], 7.2854, tolerance = 1e-3)
  # `induced` enters as a numeric trend, so the predicted OR for two prior
  # abortions is exp(2 * beta) = OR(1)^2 ~ 16.7 (the handbook two-or-more value).
  # This is a model prediction, not a separate assertion -- exp(2 * log(OR)) is
  # algebraically OR^2, so re-asserting it against 16.74 would only re-pin the
  # per-unit OR already checked above; left documented rather than tested.
})

# --- truth-based: CMLE recovers the conditional log-OR -------------------

test_that("the conditional OR recovers the matched-set log-OR (truth-based)", {
  df <- make_matched_cc(n_sets = 400L, ratio = 3L, beta_x = log(2.5))
  truth <- attr(df, "truth")
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  # The set-level exposure prevalence is conditioned away, so the CMLE recovers
  # exp(beta_x) up to sampling error. Use a SELF-SCALING band of a few reported
  # SEs rather than a fixed absolute tolerance: the estimator's Monte-Carlo
  # sampling SD in this DGP is ~0.128 (a 60-seed check), which matches the
  # reported log-OR SE, so a fixed 0.12 tolerance would be under one SD and pass
  # only by luck of the seed. A 3.5-SE band is robust to the seed yet still
  # rejects any bias above ~0.45 (a sign flip sits ~14 SEs away).
  se_log <- res$estimates$se
  expect_lt(
    abs(log(res$contrasts$estimate) - unname(truth["beta_x"])),
    3.5 * se_log
  )
})

# --- 1:1 matching: closed-form McNemar OR AND variance (independent oracle) ---

# Critical-review-loop (2026-06-03, test-audit B1/B2; repro
# /tmp/matchatr_repro_mcnemar.R). The other oracle tests compare against
# survival::clogit run with the same formula, so they validate forwarding, not
# the conditional variance. For 1:1 matching the conditional likelihood reduces
# to McNemar's: among discordant pairs n10 (case exposed, control unexposed) and
# n01 (case unexposed, control exposed), the CMLE is OR = n10/n01 with
# Var(log OR) = 1/n10 + 1/n01 -- both exact and computed here WITHOUT clogit.
test_that("1:1 matching reproduces the closed-form McNemar OR and variance", {
  df <- make_matched_cc(n_sets = 300L, ratio = 1L, beta_x = log(2.5))
  sets <- split(df, df$set)
  n10 <- sum(vapply(
    sets,
    function(s) s$x[s$case == 1L] == 1L && s$x[s$case == 0L] == 0L,
    logical(1)
  ))
  n01 <- sum(vapply(
    sets,
    function(s) s$x[s$case == 1L] == 0L && s$x[s$case == 0L] == 1L,
    logical(1)
  ))

  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  # Point OR = n10 / n01 (the discordant-pair ratio), exact.
  expect_equal(res$contrasts$estimate, n10 / n01, tolerance = 1e-7)
  # Information-matrix Var(log OR) = 1/n10 + 1/n01, exact -- this is the
  # variance check no clogit-vcov comparison can provide.
  expect_equal(res$estimates$se^2, 1 / n10 + 1 / n01, tolerance = 1e-6)
})

test_that("adjusting for a non-matching covariate recovers its known log-OR", {
  # The covariate z now carries a genuine conditional log-OR; a pure-noise z
  # (beta_z = 0) would let a mis-recovered confounder pass unnoticed.
  df <- make_matched_cc(n_sets = 600L, beta_z = log(1.8))
  truth <- attr(df, "truth")
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    confounders = ~z,
    estimator = "clogit"
  )
  # Pass-through against the hand-fit clogit (the wrapper builds the same model).
  oracle <- survival::clogit(case ~ x + z + strata(set), data = df)
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
  expect_equal(stats::vcov(fit$model), stats::vcov(oracle))
  # Truth-based: both the exposure and the adjusted covariate recover their known
  # conditional log-ORs, each within 3.5 of its own reported SE (an exposure /
  # confounder swap or a mis-recovered adjustment would fail this).
  td <- tidy(fit) # log scale, with std.error per term
  for (term in c("x", "z")) {
    est <- td$estimate[td$term == term]
    se <- td$std.error[td$term == term]
    truth_term <- if (term == "x") truth["beta_x"] else truth["beta_z"]
    expect_lt(abs(est - unname(truth_term)), 3.5 * se)
  }
})

# --- factor exposure (per-level OR vs the reference) --------------------

test_that("a factor exposure reports the per-level OR and the reference", {
  df <- infert
  df$ind <- factor(
    ifelse(df$induced >= 1L, "yes", "no"),
    levels = c("no", "yes")
  )
  fit <- matcha(
    df,
    "case",
    "ind",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  oracle <- survival::clogit(case ~ ind + strata(stratum), data = df)
  res <- contrast(fit, type = "or")
  expect_identical(res$contrasts$comparison, "indyes")
  expect_identical(res$reference, "no")
  expect_equal(
    res$contrasts$estimate,
    unname(exp(stats::coef(oracle)["indyes"])),
    tolerance = 1e-8
  )
})

# --- frequency matching: several strata columns cross into one factor ----

test_that("multi-column strata cross into a single conditioning factor", {
  df <- make_matched_cc(n_sets = 250L, ratio = 3L)
  # Split the set id into two columns whose crossing uniquely recovers the set
  # (frequency matching on two variables).
  df$set_a <- (df$set - 1L) %/% 5L
  df$set_b <- (df$set - 1L) %% 5L
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = c("set_a", "set_b")),
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ x + strata(set_a, set_b),
    data = df
  )
  expect_equal(stats::coef(fit$model), stats::coef(oracle))
})

# --- continuous exposure (per-unit OR), truth-based ---------------------

test_that("a continuous exposure recovers the per-unit log-OR (truth-based)", {
  # Matched sets with a CONTINUOUS exposure; the case is drawn from the
  # conditional likelihood (weight exp(x * beta)), so the per-unit CMLE recovers
  # beta. The set-level mean is the matched-away nuisance.
  beta <- 0.7
  df <- withr::with_seed(21L, {
    parts <- lapply(seq_len(500L), function(i) {
      x <- stats::rnorm(4L, stats::rnorm(1L, 0, 1), 1)
      case_idx <- sample.int(4L, 1L, prob = exp(x * beta))
      case <- integer(4L)
      case[case_idx] <- 1L
      data.frame(case = case, x = x, set = i)
    })
    do.call(rbind, parts)
  })
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  expect_lt(abs(log(res$contrasts$estimate) - beta), 3.5 * res$estimates$se)
})

# --- variable matching ratios in one sample -----------------------------

test_that("mixed 1:1 / 1:2 / 1:3 matching ratios recover the log-OR", {
  beta <- log(2.5)
  ratios <- rep(c(1L, 2L, 3L), each = 200L)
  df <- withr::with_seed(31L, {
    parts <- Map(
      function(i, r) {
        m <- r + 1L
        x <- stats::rbinom(m, 1L, stats::plogis(stats::rnorm(1L, 0, 1)))
        case_idx <- sample.int(m, 1L, prob = exp(x * beta))
        case <- integer(m)
        case[case_idx] <- 1L
        data.frame(case = case, x = x, set = i)
      },
      seq_along(ratios),
      ratios
    )
    do.call(rbind, parts)
  })
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  expect_lt(abs(log(res$contrasts$estimate) - beta), 3.5 * res$estimates$se)
})

# --- a concordant-on-exposure set contributes nothing -------------------

test_that("appending a set concordant on exposure leaves the OR unchanged", {
  df <- make_matched_cc(n_sets = 200L)
  or1 <- contrast(
    matcha(df, "case", "x", matched_cc(strata = "set"), estimator = "clogit"),
    type = "or"
  )$contrasts$estimate
  # A set whose members all share the same exposure has no within-set variation,
  # so its conditional-likelihood factor is constant in beta and adds nothing.
  concordant <- data.frame(
    case = c(1L, 0L, 0L, 0L),
    x = c(1L, 1L, 1L, 1L),
    z = 0L,
    set = max(df$set) + 1L
  )
  or2 <- contrast(
    matcha(
      rbind(df, concordant),
      "case",
      "x",
      matched_cc(strata = "set"),
      estimator = "clogit"
    ),
    type = "or"
  )$contrasts$estimate
  expect_equal(or1, or2, tolerance = 1e-6)
})

# --- result structure / labels ------------------------------------------

test_that("the clogit result is labelled and sized correctly", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "conditional OR")
  expect_identical(res$estimator, "clogit")
  expect_identical(res$engine, "clogit")
  # Analysis n is the rows clogit used (coxph nobs() counts events, not rows).
  expect_identical(res$n, fit$model$n)
  # contrast() with no type defaults to the OR for the clogit engine.
  expect_identical(contrast(fit)$type, "or")
})

test_that("tidy reports the clogit ORs matching the oracle", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    confounders = ~spontaneous,
    estimator = "clogit"
  )
  oracle <- survival::clogit(
    case ~ induced + spontaneous + strata(stratum),
    data = infert
  )
  td <- tidy(fit, exponentiate = TRUE)
  # No intercept row (the conditional likelihood has none), one row per term.
  expect_setequal(td$term, c("induced", "spontaneous"))
  expect_equal(
    td$estimate,
    unname(exp(stats::coef(oracle))[td$term]),
    tolerance = 1e-8
  )
})

# --- rejections ----------------------------------------------------------

test_that("RD / RR are rejected as unidentified for clogit", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
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

test_that("sandwich / bootstrap CIs are not available for clogit", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
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

test_that("an exposure with no within-set variation is not estimable", {
  df <- make_matched_cc(n_sets = 50L)
  # Make x constant within every set: it has no conditional contribution, so
  # clogit aliases its coefficient to NA and no OR is identified.
  df$x <- as.integer(df$set %% 2L)
  fit <- suppressWarnings(matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  ))
  expect_error(
    contrast(fit, type = "or"),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("a missing value triggers the dropped-rows warning with the count", {
  df <- infert
  df$induced[1:4] <- NA
  expect_warning(
    matcha(
      df,
      "case",
      "induced",
      matched_cc(strata = "stratum"),
      estimator = "clogit"
    ),
    class = "matchatr_dropped_rows",
    regexp = "4 row"
  )
})

test_that("a dropped uninformative stratum does NOT inflate the dropped count", {
  # clogit drops a stratum with no case or no control from the likelihood but
  # still counts its rows in model$n, so n_dropped = nrow - model$n stays 0 and
  # NO spurious dropped-rows warning fires (the warning is for missing values
  # only). Pin that interaction, which is easy to break by reading the wrong n.
  bad <- make_uninformative_cc()
  fired <- FALSE
  withCallingHandlers(
    matcha(
      bad,
      "case",
      "x",
      matched_cc(strata = "set"),
      estimator = "clogit"
    ),
    matchatr_dropped_rows = function(w) {
      fired <<- TRUE
      invokeRestart("muffleWarning")
    },
    matchatr_uninformative_stratum = function(w) invokeRestart("muffleWarning")
  )
  expect_false(fired)
  fit <- suppressWarnings(matcha(
    bad,
    "case",
    "x",
    matched_cc(strata = "set"),
    estimator = "clogit"
  ))
  expect_identical(fit$model$n, nrow(bad))
})

test_that("clogit rejection messages read clearly", {
  fit <- matcha(
    infert,
    "case",
    "induced",
    matched_cc(strata = "stratum"),
    estimator = "clogit"
  )
  expect_snapshot(
    contrast(fit, type = "or", ci_method = "sandwich"),
    error = TRUE
  )
})
