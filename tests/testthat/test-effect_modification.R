# Effect modification in matched case-control data: matcha(effect_modifier =)
# fits clogit with the exposure x modifier interaction and contrast(type = "or")
# reports the exposure's stratum-specific conditional OR within each modifier
# level. Oracles: the 1:1 McNemar closed form per level (point AND variance,
# independent of clogit), survival::clogit hand-built linear combinations
# (forwarding + labels + multi-level ordering), and a truth DGP with known
# per-level conditional log-ORs.

# --- independent point + variance oracle: McNemar within each level ----------

# When the modifier is constant within set (a matching variable) and matching is
# 1:1, the matched sets split into disjoint groups by level and each level's
# conditional likelihood reduces to McNemar's: OR(level) = n10/n01 with
# Var(log OR) = 1/n10 + 1/n01 over that level's discordant pairs. This validates
# BOTH the per-level point estimate and the linear-combination variance against
# a closed form computed WITHOUT survival::clogit.
test_that("per-level OR and variance match the within-level McNemar closed form", {
  df <- make_matched_cc_em(
    n_sets = 350L,
    ratio = 1L,
    betas = c(a = log(2), b = log(5)),
    within_set = FALSE
  )
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    effect_modifier = "m",
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")

  mcnemar_level <- function(sub) {
    s <- split(sub, sub$set)
    n10 <- sum(vapply(
      s,
      function(p) p$x[p$case == 1L] == 1L && p$x[p$case == 0L] == 0L,
      logical(1)
    ))
    n01 <- sum(vapply(
      s,
      function(p) p$x[p$case == 1L] == 0L && p$x[p$case == 0L] == 1L,
      logical(1)
    ))
    c(or = n10 / n01, var_log = 1 / n10 + 1 / n01)
  }
  for (lev in c("a", "b")) {
    row <- res$contrasts[res$contrasts$comparison == paste0("x | m = ", lev), ]
    est <- res$estimates[res$estimates$term == paste0("x | m = ", lev), ]
    mc <- mcnemar_level(df[df$m == lev, ])
    # Point OR = n10 / n01, exact.
    expect_equal(row$estimate, unname(mc["or"]), tolerance = 1e-7)
    # Var(log OR) = 1/n10 + 1/n01, exact -- the linear-combination variance the
    # clogit vcov alone could not independently confirm.
    expect_equal(est$se^2, unname(mc["var_log"]), tolerance = 1e-7)
  }
})

# --- pass-through oracle: hand-built clogit linear combinations ---------------

# A three-level modifier that VARIES within set (so its main effect is
# estimable, not aliased) plus a non-matching covariate. The per-level OR is the
# linear combination beta_x (+ beta_{x:level}) of the hand-fit clogit, and its
# SE comes from that fit's vcov -- this pins forwarding, the level labelling,
# and the ordering of the interaction columns against the modifier levels.
test_that("stratum-specific ORs match hand-built survival::clogit combinations", {
  df <- make_matched_cc_em(
    n_sets = 600L,
    ratio = 2L,
    betas = c(lo = log(1.5), mid = log(3), hi = log(6)),
    within_set = TRUE,
    seed = 52L
  )
  df$z <- withr::with_seed(99L, stats::rnorm(nrow(df)))
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    confounders = ~z,
    effect_modifier = "m",
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")

  oracle <- survival::clogit(case ~ x * m + z + strata(set), data = df)
  b <- stats::coef(oracle)
  V <- stats::vcov(oracle)
  # The exposure log OR at each modifier level and its variance from the joint
  # vcov: reference level is beta_x; non-reference level adds the interaction.
  level_combo <- list(
    "x | m = lo" = "x",
    "x | m = mid" = c("x", "x:mmid"),
    "x | m = hi" = c("x", "x:mhi")
  )
  for (lab in names(level_combo)) {
    terms_l <- level_combo[[lab]]
    est <- sum(b[terms_l])
    se <- sqrt(sum(V[terms_l, terms_l]))
    z <- stats::qnorm(0.975)
    row <- res$contrasts[res$contrasts$comparison == lab, ]
    er <- res$estimates[res$estimates$term == lab, ]
    expect_equal(row$estimate, unname(exp(est)), tolerance = 1e-7)
    expect_equal(er$se, unname(se), tolerance = 1e-7)
    expect_equal(
      c(row$ci_lower, row$ci_upper),
      unname(exp(est + c(-1, 1) * z * se)),
      tolerance = 1e-7
    )
  }
  expect_identical(res$reference, "lo")
})

# --- truth-based recovery of the per-level conditional log-ORs ---------------

test_that("the per-level ORs recover the known conditional log-ORs (truth)", {
  betas <- c(a = log(2), b = log(4))
  df <- make_matched_cc_em(
    n_sets = 700L,
    ratio = 3L,
    betas = betas,
    within_set = TRUE,
    seed = 63L
  )
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    effect_modifier = "m",
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  # Each level's CMLE recovers its known conditional log-OR up to sampling
  # error; use a self-scaling band of 3.5 reported SEs (mirrors the clogit
  # truth tests) so the check is robust to the seed yet rejects a real bias or a
  # level swap (the two truths sit several SEs apart).
  for (lev in names(betas)) {
    est <- res$estimates[res$estimates$term == paste0("x | m = ", lev), ]
    expect_lt(abs(est$estimate - unname(betas[lev])), 3.5 * est$se)
  }
})

# --- character / logical modifier is coerced to a factor ----------------------

test_that("a character modifier gives the same ORs as the equivalent factor", {
  df <- make_matched_cc_em(n_sets = 300L, ratio = 1L, within_set = FALSE)
  df_chr <- df
  df_chr$m <- as.character(df_chr$m)
  res_f <- contrast(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "m",
      estimator = "clogit"
    ),
    type = "or"
  )
  res_c <- contrast(
    matcha(
      df_chr,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "m",
      estimator = "clogit"
    ),
    type = "or"
  )
  expect_equal(res_c$contrasts$estimate, res_f$contrasts$estimate)
  expect_equal(res_c$estimates$se, res_f$estimates$se)
})

# --- result structure / labels -----------------------------------------------

test_that("the effect-modification result is labelled and shaped correctly", {
  df <- make_matched_cc_em(n_sets = 200L, ratio = 1L, within_set = FALSE)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    effect_modifier = "m",
    estimator = "clogit"
  )
  res <- contrast(fit, type = "or")
  expect_s3_class(res, "matchatr_result")
  expect_identical(res$type, "or")
  expect_identical(res$estimand, "stratum-specific conditional OR")
  expect_identical(res$reference, "a")
  # One row per modifier level, labelled "exposure | modifier = level".
  expect_identical(res$contrasts$comparison, c("x | m = a", "x | m = b"))
  # The stored vcov is the per-level log-OR covariance (square, named by level).
  expect_identical(dim(res$vcov), c(2L, 2L))
  expect_identical(rownames(res$vcov), c("x | m = a", "x | m = b"))
  expect_identical(res$n, fit$model$n)
  # contrast() with no type defaults to the OR for the clogit engine.
  expect_identical(contrast(fit)$type, "or")
  # tidy() renders one row per level.
  expect_identical(nrow(tidy(res)), 2L)
})

# --- rejections --------------------------------------------------------------

test_that("effect_modifier is rejected for a non-clogit engine", {
  df <- make_matched_cc_em(n_sets = 100L, within_set = TRUE)
  expect_error(
    matcha(
      df,
      "case",
      "x",
      unmatched_cc(),
      effect_modifier = "m",
      estimator = "logistic"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("a continuous (numeric) modifier is rejected", {
  df <- make_matched_cc_em(n_sets = 100L, within_set = TRUE)
  df$mnum <- as.numeric(df$m)
  expect_error(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "mnum",
      estimator = "clogit"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("a modifier coinciding with the exposure or outcome is rejected", {
  df <- make_matched_cc_em(n_sets = 100L, within_set = TRUE)
  expect_error(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "x",
      estimator = "clogit"
    ),
    class = "matchatr_bad_input"
  )
  expect_error(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "case",
      estimator = "clogit"
    ),
    class = "matchatr_bad_input"
  )
})

test_that("a 3+-level factor exposure with effect modification is rejected", {
  df <- make_matched_cc_em(n_sets = 150L, within_set = TRUE)
  df$xf <- factor(sample(c("lo", "mid", "hi"), nrow(df), replace = TRUE))
  expect_error(
    matcha(
      df,
      "case",
      "xf",
      matched_cc(strata = "set"),
      effect_modifier = "m",
      estimator = "clogit"
    ),
    class = "matchatr_unsupported_combination"
  )
})

test_that("RD / RR and sandwich / bootstrap remain rejected with a modifier", {
  df <- make_matched_cc_em(n_sets = 200L, within_set = FALSE)
  fit <- matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    effect_modifier = "m",
    estimator = "clogit"
  )
  expect_error(
    contrast(fit, type = "difference"),
    class = "matchatr_unidentified_estimand"
  )
  expect_error(
    contrast(fit, type = "or", ci_method = "sandwich"),
    class = "matchatr_unsupported_variance"
  )
})

test_that("an unestimable interaction level aborts rather than returning NA", {
  # Level "b" has no within-set exposure variation in any of its sets, so the
  # x:mb interaction is aliased to NA: the per-level OR is not identified.
  df_a <- make_matched_cc_em(n_sets = 120L, ratio = 1L, within_set = FALSE)
  df_a <- df_a[df_a$m == "a", ]
  # Append "b" sets that are all concordant on exposure (no discordant pairs).
  bset <- do.call(
    rbind,
    lapply(seq_len(40L), function(i) {
      data.frame(
        case = c(1L, 0L),
        x = c(1L, 1L),
        m = "b",
        set = 10000L + i
      )
    })
  )
  bset$m <- factor(bset$m, levels = c("a", "b"))
  df <- rbind(df_a, bset)
  fit <- suppressWarnings(matcha(
    df,
    "case",
    "x",
    matched_cc(strata = "set"),
    effect_modifier = "m",
    estimator = "clogit"
  ))
  expect_error(
    contrast(fit, type = "or"),
    class = "matchatr_unestimable_exposure"
  )
})

test_that("effect-modification rejection messages read clearly", {
  df <- make_matched_cc_em(n_sets = 100L, within_set = TRUE)
  df$mnum <- as.numeric(df$m)
  expect_snapshot(
    matcha(
      df,
      "case",
      "x",
      matched_cc(strata = "set"),
      effect_modifier = "mnum",
      estimator = "clogit"
    ),
    error = TRUE
  )
})
