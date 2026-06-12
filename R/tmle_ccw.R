#' Fit a case-control-weighted targeted maximum likelihood estimator
#'
#' Implements CCW-TMLE: the case-control sample is reweighted to the source
#' population with the Rose & van der Laan case-control weights, and a targeted
#' maximum likelihood estimator of the marginal effect is fitted on the weighted
#' sample. This is the one genuinely new estimator in the case-control-weighted
#' family — the etverse has no targeted-learning engine, so the targeting
#' (fluctuation) step is matchatr's own.
#'
#' @details
#' The algorithm (van der Laan & Rubin 2006; van der Laan & Rose 2011), with the
#' case-control weights w as observation weights throughout:
#'
#' 1. Initial outcome model Q̄⁰(A, W) = Ê(Y | A, W) — a weighted logistic fit;
#'    predict Q̄⁰(1, W), Q̄⁰(0, W), Q̄⁰(A, W).
#' 2. Propensity g(W) = P̂(A = 1 | W) — a weighted logistic fit, bounded away
#'    from 0 / 1.
#' 3. Clever covariate H(A, W) = A / g(W) − (1 − A) / (1 − g(W)); so
#'    H(1, W) = 1 / g(W) and H(0, W) = −1 / (1 − g(W)).
#' 4. Fluctuation: weighted logistic regression of Y on H(A, W) with offset
#'    logit Q̄⁰(A, W) and no intercept, giving the tilt ε.
#' 5. Update Q̄*(a, W) = expit(logit Q̄⁰(a, W) + ε H(a, W)).
#' 6. Marginalize: ψ̂₁ = Σ w Q̄*(1, W) / Σ w, ψ̂₀ = Σ w Q̄*(0, W) / Σ w.
#'
#' The targeted Q̄* solves the case-control-weighted efficient-influence-equation
#' score Σ w H (Y − Q̄*) = 0, so the plug-in marginal effect is asymptotically
#' efficient and **doubly robust** (consistent if either Q̄⁰ or g is correct).
#'
#' @param fit A `matchatr_fit` whose `engine` is `"ccw_tmle"`, carrying the
#'   case-control `data`, the `outcome` / `exposure` / `confounders` roles, and a
#'   `design` whose `prevalence` (q0) is set.
#' @returns A `matchatr_ccw_tmle` object holding the targeted treatment-specific
#'   means (`EY1`, `EY0`), the per-subject efficient-influence-function components
#'   (`D1`, `D0`), the case-control `weights`, the sample size `n`, and the tilt
#'   `eps`; [contrast()] turns it into the marginal effect with EIF variance.
#' @family estimators
#' @seealso [matcha()], [contrast()], `ccw_prepare()`, [tmle::tmle()]
#' @noRd
fit_ccw_tmle <- function(fit) {
  prep <- ccw_prepare(fit)
  dt <- prep$dt
  w <- as.numeric(prep$weights)

  # The outcome / propensity model fitter (NULL -> stats::glm); a user may pass
  # mgcv::gam for smooth confounder adjustment of both nuisance models.
  model_fn <- fit$details$model_fn
  if (is.null(model_fn)) {
    model_fn <- stats::glm
  }

  ccw_tmle_target(
    dt = dt,
    outcome = fit$outcome,
    exposure = fit$exposure,
    confounders = fit$confounders,
    weights = w,
    model_fn = model_fn
  )
}

#' Core CCW-TMLE targeting computation
#'
#' Runs the initial outcome / propensity fits, the clever-covariate logistic
#' fluctuation, the update, and the marginalization, returning the targeted
#' treatment-specific means and their efficient-influence-function components.
#'
#' @param dt A `data.table` with the outcome and exposure recoded to 0/1.
#' @param outcome,exposure Character scalars naming the outcome / exposure columns.
#' @param confounders A one-sided formula of confounders.
#' @param weights Numeric case-control weights, one per row of `dt`.
#' @param model_fn The nuisance-model fitter (e.g. `stats::glm`).
#' @param gbound Numeric in (0, 0.5): the propensity is bounded to
#'   `[gbound, 1 - gbound]` so the clever covariate cannot blow up under a
#'   near-positivity violation (the `tmle` package default, 0.025).
#' @param ybound Numeric: the outcome predictions are bounded to
#'   `[ybound, 1 - ybound]` so `logit` is finite.
#' @returns A `matchatr_ccw_tmle` object (see `fit_ccw_tmle()`).
#' @family estimators
#' @seealso `fit_ccw_tmle()`
#' @noRd
ccw_tmle_target <- function(
  dt,
  outcome,
  exposure,
  confounders,
  weights,
  model_fn,
  gbound = 0.025,
  ybound = 1e-5
) {
  y <- as.numeric(dt[[outcome]])
  a <- as.numeric(dt[[exposure]])
  conf_terms <- attr(stats::terms(confounders), "term.labels")

  # Initial outcome model Q̄⁰(A, W): weighted logistic Y ~ A + W. quasibinomial
  # fits the same mean model as binomial but is silent on the non-integer
  # "successes" the fractional case-control weights produce.
  q_formula <- stats::reformulate(c(exposure, conf_terms), response = outcome)
  qfit <- model_fn(
    q_formula,
    data = dt,
    family = stats::quasibinomial(),
    weights = weights
  )
  dt1 <- data.table::copy(dt)
  dt1[[exposure]] <- 1
  dt0 <- data.table::copy(dt)
  dt0[[exposure]] <- 0
  q_aw <- bound01(stats::predict(qfit, type = "response"), ybound)
  q_1w <- bound01(
    stats::predict(qfit, newdata = dt1, type = "response"),
    ybound
  )
  q_0w <- bound01(
    stats::predict(qfit, newdata = dt0, type = "response"),
    ybound
  )

  # Propensity g(W) = P(A = 1 | W): weighted logistic A ~ W, bounded.
  g_formula <- stats::reformulate(conf_terms, response = exposure)
  gfit <- model_fn(
    g_formula,
    data = dt,
    family = stats::quasibinomial(),
    weights = weights
  )
  gw <- bound01(stats::predict(gfit, type = "response"), gbound)

  # Clever covariate H(A, W) = A/g − (1−A)/(1−g); H(1,W) = 1/g, H(0,W) = −1/(1−g).
  h_aw <- a / gw - (1 - a) / (1 - gw)
  h_1w <- 1 / gw
  h_0w <- -1 / (1 - gw)

  # Fluctuation: weighted logistic Y ~ -1 + H with offset logit Q̄⁰(A,W). The
  # single coefficient eps is the tilt that solves the weighted EIF score
  # Σ w H (Y − Q̄*) = 0. A degenerate (non-finite) tilt means no update.
  flucfit <- stats::glm(
    y ~ -1 + h_aw,
    offset = stats::qlogis(q_aw),
    family = stats::quasibinomial(),
    weights = weights
  )
  eps <- unname(stats::coef(flucfit))
  if (!is.finite(eps)) {
    # A non-finite tilt means the fluctuation did not converge (e.g. a degenerate
    # clever covariate); fall back to the untargeted initial fit rather than
    # propagate NaN, and surface it so the estimate is not silently un-targeted.
    rlang::warn(
      c(
        "The CCW-TMLE fluctuation did not converge; using the untargeted initial fit.",
        i = "The estimate reverts to the initial outcome model (no targeting step)."
      ),
      class = c("matchatr_tmle_convergence", "matchatr_warning")
    )
    eps <- 0
  }

  # Update Q̄* = expit(logit Q̄⁰ + eps H).
  qs_aw <- stats::plogis(stats::qlogis(q_aw) + eps * h_aw)
  qs_1w <- stats::plogis(stats::qlogis(q_1w) + eps * h_1w)
  qs_0w <- stats::plogis(stats::qlogis(q_0w) + eps * h_0w)

  # Marginal targeted treatment-specific means (weighted to the population).
  ey1 <- stats::weighted.mean(qs_1w, weights)
  ey0 <- stats::weighted.mean(qs_0w, weights)

  # Efficient-influence-function components for each mean:
  #   D₁ = (A/g)(Y − Q̄*(A,W)) + Q̄*(1,W) − ψ₁,
  #   D₀ = ((1−A)/(1−g))(Y − Q̄*(A,W)) + Q̄*(0,W) − ψ₀.
  d1 <- (a / gw) * (y - qs_aw) + qs_1w - ey1
  d0 <- ((1 - a) / (1 - gw)) * (y - qs_aw) + qs_0w - ey0

  structure(
    list(
      EY1 = ey1,
      EY0 = ey0,
      D1 = d1,
      D0 = d0,
      weights = weights,
      n = length(y),
      eps = eps
    ),
    class = "matchatr_ccw_tmle"
  )
}

#' Assemble the CCW-TMLE marginal contrast
#'
#' Turns a targeted CCW-TMLE fit into a `matchatr_result` reporting the marginal
#' risk difference (`type = "difference"`), risk ratio (`type = "ratio"`), or
#' marginal odds ratio (`type = "or"`) with the efficient-influence-function
#' variance weighted by the case-control weights.
#'
#' @details
#' With the case-control weights treated as fixed (known q0), the variance of a
#' marginal mean is the weighted EIF variance Var(ψ̂) = Σ (wᵢ Dᵢ)² / n², where n
#' is the sample size and Σ wᵢ = n. The risk difference uses D = D₁ − D₀; the
#' risk ratio and odds ratio use the delta-method log-scale influence functions
#' D₁/ψ₁ − D₀/ψ₀ and D₁/(ψ₁(1−ψ₁)) − D₀/(ψ₀(1−ψ₀)), with the interval formed on
#' the log scale and exponentiated. A bootstrap interval (which must resample
#' within the case / control strata and recompute the weights) is deferred.
#'
#' @param fit A `matchatr_fit` whose `model` is a `matchatr_ccw_tmle` object.
#' @param type Character contrast scale: `"difference"`, `"ratio"`, or `"or"`.
#' @param ci_method Character variance source; `"model"` / `"sandwich"` both use
#'   the EIF variance (recorded as `"sandwich"`); `"bootstrap"` is rejected.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` carrying the targeted intervention means and the
#'   marginal contrast with EIF variance.
#' @family estimators
#' @seealso [contrast()], `fit_ccw_tmle()`
#' @noRd
contrast_ccw_tmle <- function(
  fit,
  type,
  ci_method,
  conf_level,
  call = rlang::caller_env()
) {
  if (!type %in% c("difference", "ratio", "or")) {
    rlang::abort(
      c(
        paste0(
          "A case-control-weighted estimator reports a marginal effect, not `type = \"",
          type,
          "\"`."
        ),
        i = paste0(
          'Use `type = "difference"` (risk difference), `"ratio"` (risk ',
          'ratio), or `"or"` (marginal odds ratio).'
        )
      ),
      class = c("matchatr_unidentified_estimand", "matchatr_error"),
      call = call
    )
  }
  if (identical(ci_method, "bootstrap")) {
    rlang::abort(
      c(
        '`ci_method = "bootstrap"` is not available for the case-control-weighted estimators.',
        i = paste0(
          'Use `ci_method = "model"` or `ci_method = "sandwich"` (the efficient ',
          "influence-function variance)."
        )
      ),
      class = c("matchatr_unsupported_variance", "matchatr_error"),
      call = call
    )
  }

  m <- fit$model
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  ey1 <- m$EY1
  ey0 <- m$EY0
  # EIF SE of a weighted marginal mean (or any linear combination of the D's):
  # sqrt(Σ (w D)²) / n, with Σ w = n (the case-control weights treated as fixed).
  se_eif <- function(d) sqrt(sum((m$weights * d)^2)) / m$n

  se_y1 <- se_eif(m$D1)
  se_y0 <- se_eif(m$D0)
  estimates <- data.table::data.table(
    intervention = c("treated", "control"),
    estimate = c(ey1, ey0),
    se = c(se_y1, se_y0),
    ci_lower = c(ey1 - z * se_y1, ey0 - z * se_y0),
    ci_upper = c(ey1 + z * se_y1, ey0 + z * se_y0)
  )

  if (identical(type, "difference")) {
    est <- ey1 - ey0
    se <- se_eif(m$D1 - m$D0)
    lower <- est - z * se
    upper <- est + z * se
    estimand <- "marginal risk difference"
  } else if (identical(type, "ratio")) {
    log_est <- log(ey1) - log(ey0)
    se_log <- se_eif(m$D1 / ey1 - m$D0 / ey0)
    est <- exp(log_est)
    lower <- exp(log_est - z * se_log)
    upper <- exp(log_est + z * se_log)
    # OR/RR-scale `se` is the delta-method value (RR * SE(log RR)); the interval
    # is the log-scale Wald exponentiated, so `se` does not reconstruct it.
    se <- est * se_log
    estimand <- "marginal risk ratio"
  } else {
    log_est <- stats::qlogis(ey1) - stats::qlogis(ey0)
    se_log <- se_eif(m$D1 / (ey1 * (1 - ey1)) - m$D0 / (ey0 * (1 - ey0)))
    est <- exp(log_est)
    lower <- exp(log_est - z * se_log)
    upper <- exp(log_est + z * se_log)
    se <- est * se_log
    estimand <- "marginal odds ratio"
  }

  new_matchatr_result(
    estimates = estimates,
    contrasts = data.table::data.table(
      comparison = "treated vs control",
      estimate = est,
      se = se,
      ci_lower = lower,
      ci_upper = upper
    ),
    type = type,
    estimand = estimand,
    # The EIF plug-in variance is the influence-function / sandwich variance, the
    # same family the other CCW engines report.
    ci_method = "sandwich",
    reference = "control",
    n = m$n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = NULL,
    call = call
  )
}
