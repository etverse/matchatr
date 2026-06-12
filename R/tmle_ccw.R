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
#'   `[gbound, 1 - gbound]` so the clever covariate stays finite (the `tmle`
#'   package default, 0.025).
#' @param ybound Numeric: outcome predictions are bounded to `[ybound, 1 - ybound]`
#'   so `logit` is finite.
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
