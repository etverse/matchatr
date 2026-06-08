#' Test homogeneity of an exposure's odds ratios across disease subtypes
#'
#' @description
#' Given a fitted polytomous (multinomial) case-control model from [matcha()]
#' (`estimator = "polytomous"`), tests whether the exposure acts the same way on
#' every disease subtype — the etiologic-heterogeneity question — and reports the
#' efficient pooled ("common") odds ratio that holds under homogeneity. For each
#' exposure term the null hypothesis is that its log odds ratio is equal across
#' the non-reference outcome groups (H0: beta_1 = beta_2 = ... = beta_M).
#'
#' @details
#' The test is the Wald test of the canonical etiologic-heterogeneity analysis
#' (Begg & Gray, 1984; as implemented in `riskclustr::eh_test_subtype`): with the
#' stacked subtype log odds ratios `b` (length M = number of non-reference
#' groups) and their joint covariance `V` from the multinomial information matrix,
#' and a full-rank contrast matrix `C` (M - 1 rows) that differences consecutive
#' subtypes,
#'
#'   W = (C b)' (C V C')^-1 (C b)  ~  chi-squared with M - 1 degrees of freedom.
#'
#' The common odds ratio is the minimum-variance (generalized-least-squares /
#' inverse-variance) combination of the subtype log odds ratios — the restricted
#' estimator under the equality constraint, asymptotically equivalent to the
#' constrained maximum-likelihood fit:
#'
#'   b_common = (1' V^-1 b) / (1' V^-1 1),   Var(b_common) = 1 / (1' V^-1 1),
#'
#' exponentiated to the odds-ratio scale with a Wald interval on the log scale
#' (so the interval is asymmetric on the OR scale). Because the constraint is
#' imposed on the already-fitted unconstrained model, no refit is needed and the
#' test handles continuous confounders directly. The pooled estimate is more
#' efficient than any single subtype odds ratio (Begg & Gray, 1984): its standard
#' error is smaller than each pooled term's.
#'
#' Each exposure term is tested separately (one "risk factor" per column): a
#' binary or continuous exposure contributes one row, an unordered factor
#' exposure one row per non-reference level. The fit must be the polytomous
#' multinomial engine (three or more outcome groups); any other engine — or a
#' fit that produced no model — is rejected.
#'
#' @param fit A `matchatr_fit` returned by [matcha()] whose `engine` is
#'   `"multinom"` (i.e. `estimator = "polytomous"`).
#' @param conf_level Numeric confidence level for the common-OR interval, a
#'   single number strictly in (0, 1). Defaults to 0.95.
#'
#' @returns A `matchatr_homogeneity` object: a list carrying `homogeneity` (a
#'   `data.table` with one row per exposure term — the term, the common odds
#'   ratio with its Wald bounds, and the homogeneity chi-squared `statistic`,
#'   `df`, and `p.value`), `subtype` (the per-subtype odds ratios it pools),
#'   the baseline `reference` group, the `conf_level`, the analysis size `n`,
#'   and the `estimator` / `engine` labels.
#'
#' @examples
#' set.seed(5)
#' n <- 4000
#' x <- rbinom(n, 1, 0.4)
#' # Subtype A and B share the exposure effect (homogeneity holds).
#' eta <- cbind(control = 0, A = -1 + log(2) * x, B = -1.4 + log(2) * x)
#' prob <- exp(eta) / rowSums(exp(eta))
#' g <- apply(prob, 1, function(p) sample(c("control", "A", "B"), 1, prob = p))
#' d <- data.frame(g = g, x = x)
#' fit <- matcha(d, outcome = "g", exposure = "x",
#'               design = unmatched_cc(), estimator = "polytomous",
#'               reference = "control")
#' test_homogeneity(fit)
#'
#' @seealso [matcha()], [contrast()], [tidy.matchatr_homogeneity()]
#' @references
#' Begg CB, Gray R (1984). Calculation of polytomous logistic regression
#' parameters using individualized regressions. *Biometrika* 71(1), 11-18.
#'
#' Borgan O, Breslow N, Chatterjee N, Gail MH, Scott A, Wild CJ (2018).
#' *Handbook of Statistical Methods for Case-Control Studies*, Chapter 5.
#' @family estimators
#' @export
test_homogeneity <- function(fit, conf_level = 0.95) {
  if (!inherits(fit, "matchatr_fit")) {
    rlang::abort(
      "`fit` must be a `matchatr_fit` object returned by `matcha()`.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  # An engine with no wired estimator leaves `model = NULL`; there is nothing to
  # test until its estimation layer has run.
  require_estimated(fit)
  # Homogeneity across subtypes is only defined for the multinomial fit, which is
  # the only engine with more than one non-reference outcome group to compare.
  if (!identical(fit$engine, "multinom")) {
    rlang::abort(
      c(
        "Homogeneity of subtype odds ratios is only defined for the polytomous estimator.",
        i = paste0(
          "Use `estimator = \"polytomous\"` with a multi-group outcome; ",
          "got engine `",
          fit$engine,
          "`."
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  check_conf_level(conf_level)

  call <- match.call()
  model <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Stacked subtype exposure log ORs and their joint covariance, reusing the
  # contrast engine the polytomous OR result is built from. The flat vectors are
  # subtype-major then exposure-column, so each exposure column's M subtype
  # entries sit at a fixed stride.
  ex <- multinom_exposure_or(
    model,
    fit$exposure,
    conf_level = conf_level,
    call = call
  )
  n_sub <- length(ex$subtypes)
  n_cols <- length(ex$predictors)
  # A multinomial fit always has at least two non-reference groups (>= 3 outcome
  # groups), so the homogeneity contrast has >= 1 degree of freedom; guard
  # defensively so a malformed fit cannot reach the singular contrast below.
  if (n_sub < 2L) {
    rlang::abort(
      c(
        "Homogeneity needs at least two non-reference subtypes to compare.",
        i = "The polytomous fit has fewer than two subtypes."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }

  # One homogeneity test + pooled OR per exposure column (each a separate "risk
  # factor"). Column c's M subtype entries are at flat positions c, c + n_cols,
  # c + 2 n_cols, ... (subtype-major layout).
  rows <- lapply(seq_len(n_cols), function(c) {
    idx <- seq(c, by = n_cols, length.out = n_sub)
    homogeneity_one_term(
      beta = ex$log_or[idx],
      vcov = ex$vcov[idx, idx, drop = FALSE],
      term = ex$predictors[c],
      z = z,
      call = call
    )
  })

  homogeneity <- data.table::rbindlist(lapply(rows, function(r) {
    data.table::data.table(
      term = r$term,
      common_or = r$common_or,
      ci_lower = r$ci_lower,
      ci_upper = r$ci_upper,
      statistic = r$statistic,
      df = r$df,
      p.value = r$p.value
    )
  }))

  # The per-subtype odds ratios the pooled estimate combines, for reference; the
  # Wald interval is on the log scale, exponentiated (asymmetric on the OR
  # scale), matching the polytomous contrast convention.
  subtype <- data.table::data.table(
    comparison = ex$comparison,
    or = exp(ex$log_or),
    ci_lower = exp(ex$log_or - z * ex$se),
    ci_upper = exp(ex$log_or + z * ex$se)
  )

  new_matchatr_homogeneity(
    homogeneity = homogeneity,
    subtype = subtype,
    # The baseline outcome group every subtype OR (and the pooled OR) is taken
    # against (multinom's first outcome level).
    reference = model$lev[1],
    conf_level = conf_level,
    n = nrow(model$residuals),
    estimator = fit$estimator,
    engine = fit$engine,
    call = call
  )
}

#' Homogeneity Wald test and pooled odds ratio for one exposure term
#'
#' The per-exposure-column kernel of [test_homogeneity()]: from the stacked
#' subtype log odds ratios and their joint covariance it forms the Wald
#' chi-squared statistic for equality across subtypes and the
#' generalized-least-squares pooled (common) log odds ratio.
#'
#' @details
#' With `b` the M subtype log odds ratios, `V` their covariance, and a
#' consecutive-difference contrast `C` (`homogeneity_contrast_matrix()`), the
#' statistic is `(C b)' (C V C')^-1 (C b)` on `M - 1` degrees of freedom, and the
#' pooled log odds ratio is `(1' V^-1 b) / (1' V^-1 1)` with variance
#' `1 / (1' V^-1 1)`. A singular contrast covariance or covariance matrix means
#' the subtype odds ratios are not jointly estimable, which aborts with
#' `matchatr_unestimable_exposure` rather than returning a degenerate statistic.
#'
#' @param beta Numeric vector of the M subtype log odds ratios for this exposure
#'   column.
#' @param vcov Their `M x M` covariance sub-matrix (from the multinomial
#'   information matrix).
#' @param term Character label for the exposure column.
#' @param z Numeric Wald critical value for the common-OR interval.
#' @param call Caller environment surfaced in any error.
#' @returns A list with `term`, `common_or`, `ci_lower`, `ci_upper`,
#'   `statistic` (the Wald chi-squared), `df`, and `p.value`.
#' @family estimators
#' @noRd
homogeneity_one_term <- function(
  beta,
  vcov,
  term,
  z,
  call = rlang::caller_env()
) {
  m <- length(beta)
  df <- m - 1L

  unestimable <- function() {
    rlang::abort(
      c(
        paste0(
          "The subtype odds ratios for `",
          term,
          "` are not jointly estimable."
        ),
        i = "Their covariance matrix is singular, so the homogeneity test is undefined."
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }
  # `solve()` only errors on an *exactly* singular matrix; a near-singular one
  # returns large-but-finite garbage. Reject on the reciprocal condition number
  # so a degenerate subtype covariance is caught, not silently inverted.
  invertible <- function(mat) {
    all(is.finite(mat)) && rcond(mat) >= .Machine$double.eps
  }

  # Wald statistic for equality across subtypes: any full-rank contrast spanning
  # the complement of the all-equal direction gives the same W, so the
  # consecutive-difference contrast is used.
  cmat <- homogeneity_contrast_matrix(m)
  cb <- cmat %*% beta
  cvc <- cmat %*% vcov %*% t(cmat)
  if (!invertible(cvc)) {
    unestimable()
  }
  cvc_inv <- solve(cvc)
  statistic <- as.numeric(t(cb) %*% cvc_inv %*% cb)
  p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)

  # GLS / inverse-variance pooled common log OR (the restricted estimator).
  if (!invertible(vcov)) {
    unestimable()
  }
  vinv <- solve(vcov)
  ones <- rep(1, m)
  denom <- as.numeric(t(ones) %*% vinv %*% ones)
  common_log_or <- as.numeric(t(ones) %*% vinv %*% beta) / denom
  se_common <- sqrt(1 / denom)

  list(
    term = term,
    common_or = exp(common_log_or),
    ci_lower = exp(common_log_or - z * se_common),
    ci_upper = exp(common_log_or + z * se_common),
    statistic = statistic,
    df = df,
    p.value = p_value
  )
}

#' Consecutive-difference contrast matrix
#'
#' Builds the `(k - 1) x k` contrast that differences consecutive elements (row
#' `i` is `+1` at column `i` and `-1` at column `i + 1`). Applied to the stacked
#' subtype log odds ratios it spans the orthogonal complement of the all-equal
#' direction, so `C b = 0` exactly when every subtype shares one odds ratio —
#' the homogeneity null.
#'
#' @param k Integer number of subtypes (`>= 2`).
#' @returns A `(k - 1) x k` numeric matrix.
#' @family estimators
#' @noRd
homogeneity_contrast_matrix <- function(k) {
  cmat <- matrix(0, nrow = k - 1L, ncol = k)
  for (i in seq_len(k - 1L)) {
    cmat[i, i] <- 1
    cmat[i, i + 1L] <- -1
  }
  cmat
}

#' Print a matchatr homogeneity test
#'
#' @description
#' Displays a compact summary of a `matchatr_homogeneity` object: the estimator
#' and engine, the baseline reference group, the analysis size, and the table of
#' per-exposure common (pooled) odds ratios with their homogeneity Wald
#' chi-squared statistic, degrees of freedom, and p-value.
#'
#' @param x A `matchatr_homogeneity` object returned by [test_homogeneity()].
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @examples
#' set.seed(5)
#' n <- 2000
#' x <- rbinom(n, 1, 0.4)
#' eta <- cbind(control = 0, A = -1 + log(2) * x, B = -1.4 + log(2) * x)
#' prob <- exp(eta) / rowSums(exp(eta))
#' g <- apply(prob, 1, function(p) sample(c("control", "A", "B"), 1, prob = p))
#' d <- data.frame(g = g, x = x)
#' fit <- matcha(d, outcome = "g", exposure = "x",
#'               design = unmatched_cc(), estimator = "polytomous",
#'               reference = "control")
#' print(test_homogeneity(fit))
#' @seealso [test_homogeneity()]
#' @export
print.matchatr_homogeneity <- function(x, ...) {
  cat("<matchatr_homogeneity>\n")
  cat(
    " Estimator:  ",
    x$estimator,
    "  (engine: ",
    x$engine,
    ")\n",
    sep = ""
  )
  cat(" Test:       Homogeneity of subtype odds ratios (Wald)\n")
  cat(" Reference:  ", x$reference, "\n", sep = "")
  cat(" N:          ", x$n, "\n", sep = "")
  cat(
    "\nCommon (pooled) odds ratio per exposure term and homogeneity test:\n"
  )
  print(x$homogeneity)
  # Show the per-subtype odds ratios the pooled estimate combines, so the print
  # exposes what the common OR is pooling and what the test compares.
  cat("\nPer-subtype odds ratios (pooled):\n")
  print(x$subtype)
  cat(
    "\nA small p-value is evidence the exposure odds ratio differs across ",
    "subtypes.\n",
    sep = ""
  )
  invisible(x)
}

#' Tidy a matchatr homogeneity test
#'
#' Returns the per-exposure homogeneity results carried by a
#' `matchatr_homogeneity` object as a tidy `data.table`, one row per exposure
#' term, with the common (pooled) odds ratio, its confidence bounds, and the
#' homogeneity Wald chi-squared statistic, degrees of freedom, and p-value.
#'
#' @param x A `matchatr_homogeneity` object returned by [test_homogeneity()].
#' @param ... Unused; present for generic consistency.
#' @returns A `data.table` with columns `term`, `common_or`, `conf.low`,
#'   `conf.high`, `statistic`, `df`, and `p.value`.
#' @examples
#' set.seed(5)
#' n <- 2000
#' x <- rbinom(n, 1, 0.4)
#' eta <- cbind(control = 0, A = -1 + log(2) * x, B = -1.4 + log(2) * x)
#' prob <- exp(eta) / rowSums(exp(eta))
#' g <- apply(prob, 1, function(p) sample(c("control", "A", "B"), 1, prob = p))
#' d <- data.frame(g = g, x = x)
#' fit <- matcha(d, outcome = "g", exposure = "x",
#'               design = unmatched_cc(), estimator = "polytomous",
#'               reference = "control")
#' tidy(test_homogeneity(fit))
#' @seealso [test_homogeneity()]
#' @family tidiers
#' @export
tidy.matchatr_homogeneity <- function(x, ...) {
  out <- data.table::data.table(
    term = x$homogeneity$term,
    common_or = x$homogeneity$common_or,
    conf.low = x$homogeneity$ci_lower,
    conf.high = x$homogeneity$ci_upper,
    statistic = x$homogeneity$statistic,
    df = x$homogeneity$df,
    p.value = x$homogeneity$p.value
  )
  out[]
}
