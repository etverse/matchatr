#' Absolute risk from a sampled-cohort survival fit
#'
#' @description
#' Estimates the cumulative incidence F_x(t) from a fitted survival design. Three
#' engines are supported: the case-cohort pseudo-likelihood (`"cch"`, Borgan &
#' Liestøl 1990) and the IPW nested case-control weighted Cox (`"ipw_cox"`,
#' Samuelsen-weighted Breslow over the reused controls) both build
#' F_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t)) from an inverse-probability-weighted (IPW)
#' Breslow cumulative baseline hazard Λ̂₀(t); the IPW nested case-control weighted
#' Weibull AFT (`"ipw_aft"`) builds F_x(t) = 1 − exp(−exp((log t − η)/σ̂))
#' directly from the fitted parametric survival curve (η is the AFT linear
#' predictor, σ̂ the scale). Pointwise confidence intervals use the delta method on
#' the complementary log-log scale (the "log-log" CI for the survival function,
#' inverted to the risk scale).
#'
#' @param fit A `matchatr_fit` object returned by [matcha()]. The case-cohort
#'   (`"cch"`), IPW nested case-control weighted Cox (`"ipw_cox"`), and IPW nested
#'   case-control Weibull AFT (`"ipw_aft"`) engines are supported.
#' @param ... Unused; present for S3 consistency.
#' @returns A `matchatr_absolute_risk` object. See `absolute_risk.matchatr_fit`
#'   for details on the return structure.
#' @examples
#' \dontrun{
#' fit <- matcha(cohort, outcome = "d", exposure = "x",
#'               design = case_cohort(subcohort = "sc", time = "t"),
#'               confounders = ~z, estimator = "cch")
#' absolute_risk(fit, newdata = data.frame(x = 1, z = 0), times = c(1, 2, 3))
#' }
#' @family contrasts
#' @seealso [matcha()], [contrast()], [case_cohort()], [sample_ncc()]
#' @export
absolute_risk <- function(fit, ...) UseMethod("absolute_risk")

#' @rdname absolute_risk
#'
#' @param newdata A data frame of covariate values at which to evaluate F_x(t).
#'   Must contain columns matching the exposure and confounders used in `fit`.
#'   Each row is one covariate pattern; the result table carries a `row` column
#'   that indexes back to `newdata`.
#' @param times Non-empty numeric vector of evaluation times. Duplicates are
#'   dropped; times are sorted before evaluation. Times before the first event
#'   return F̂ = 0; times after the last event return the last Breslow value
#'   (step-function extrapolation).
#' @param conf_level Numeric in (0, 1). Confidence level for the pointwise
#'   intervals. Default `0.95`.
#'
#' @returns A `matchatr_absolute_risk` list with elements:
#'   - `$estimates`: a `data.table` with columns `row` (newdata row index),
#'     `time`, `estimate` (F̂_x(t)), `ci_lower`, `ci_upper` (delta-method CI
#'     on the probability scale).
#'   - `$times`: the sorted evaluation times.
#'   - `$newdata`: the supplied covariate patterns.
#'   - `$conf_level`, `$ci_method`, `$engine`, `$estimator`, `$method`.
#'
#' @export
absolute_risk.matchatr_fit <- function(
  fit,
  newdata,
  times,
  conf_level = 0.95,
  ...
) {
  # The cumulative-incidence estimator needs a fitted survival model with a
  # well-defined survival curve S(t | x). The Cox-type engines (case-cohort
  # pseudo-likelihood, IPW nested case-control weighted Cox) supply a cumulative
  # baseline hazard; the IPW AFT supplies a parametric Weibull curve.
  supported <- c("cch", "ipw_cox", "ipw_aft")
  if (!fit$engine %in% supported) {
    rlang::abort(
      c(
        paste0(
          "`absolute_risk()` is not implemented for the `",
          fit$engine,
          "` engine."
        ),
        i = paste0(
          "Supported engines: the case-cohort (`cch`), IPW nested ",
          "case-control weighted Cox (`ipw_cox`), and IPW nested case-control ",
          "Weibull AFT (`ipw_aft`) survival fits."
        )
      ),
      class = c("matchatr_not_implemented", "matchatr_error")
    )
  }
  if (is.null(fit$model)) {
    rlang::abort(
      "The fit has no estimated model. Run `matcha()` before `absolute_risk()`.",
      class = c("matchatr_not_estimated", "matchatr_error")
    )
  }
  if (!is.data.frame(newdata) || nrow(newdata) == 0L) {
    rlang::abort(
      "`newdata` must be a non-empty data frame.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  if (!is.numeric(times) || length(times) == 0L || any(!is.finite(times))) {
    rlang::abort(
      "`times` must be a non-empty numeric vector of finite values.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  switch(
    fit$engine,
    cch = absolute_risk_cch(
      fit,
      newdata = newdata,
      times = times,
      conf_level = conf_level
    ),
    ipw_cox = absolute_risk_ncc(
      fit,
      newdata = newdata,
      times = times,
      conf_level = conf_level
    ),
    ipw_aft = absolute_risk_aft(
      fit,
      newdata = newdata,
      times = times,
      conf_level = conf_level
    )
  )
}

#' Assemble an absolute-risk result from a Breslow step function
#'
#' Shared core for the case-cohort and IPW nested case-control absolute-risk
#' engines: given a (cumulative baseline hazard, log-hazard variance) step
#' function and the fitted coefficients, it evaluates
#' F_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t)) for each row of `newdata` at each requested
#' time and attaches a delta-method complementary-log-log CI. Both engines differ
#' only in how Λ̂₀(t) is built (`ipw_breslow_cch()` vs `ipw_breslow_ncc()`); the
#' covariate linear predictor, the delta-method variance, and the result
#' structure are identical, so they live here.
#'
#' @details
#' The CI is on the log(−log(1 − F)) = log(Λ_x(t)) scale, whose variance is the
#' covariate part `x' V_β x` plus the baseline part `Var(log Λ̂₀(t))`. The
#' interval is exponentiated back to a survival and inverted to the risk scale,
#' so it is asymmetric on the probability scale and stays within the unit interval.
#'
#' @param fit A `matchatr_fit` supplying the `exposure` / `confounders` / `data`
#'   columns and the `engine` / `estimator` labels recorded on the result.
#' @param newdata A data frame of covariate patterns.
#' @param times Numeric vector of evaluation times (sorted, de-duplicated here).
#' @param beta Named numeric vector of fitted coefficients (`coef(fit$model)`).
#' @param vcov_beta Variance-covariance matrix of `beta` (the engine's reported
#'   variance: the cch pseudo-likelihood variance, or the IPW robust sandwich).
#' @param breslow A list from an `ipw_breslow_*()` helper with `$times`,
#'   `$cumhaz`, and `$var_log_cumhaz` (each starting with a `t = 0` fence post).
#' @param conf_level Numeric confidence level.
#' @param method Character label recorded on the result (the pseudo-likelihood
#'   variant for cch, or `"IPW"` for the nested case-control weighted Cox).
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @noRd
assemble_absolute_risk <- function(
  fit,
  newdata,
  times,
  beta,
  vcov_beta,
  breslow,
  conf_level,
  method
) {
  times <- sort(unique(times))
  z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  # Interpolate Breslow and its variance to the requested evaluation times.
  # method = "constant" with rule = 2 gives the step-function convention: before
  # the first event -> cumhaz = 0, after the last event -> last cumhaz value.
  cumhaz_t <- stats::approx(
    breslow$times,
    breslow$cumhaz,
    xout = times,
    method = "constant",
    f = 0,
    rule = 2
  )$y
  var_log_t <- stats::approx(
    breslow$times,
    breslow$var_log_cumhaz,
    xout = times,
    method = "constant",
    f = 0,
    rule = 2
  )$y
  # Before the first event, the cumhaz is 0 so the variance is also 0.
  cumhaz_t[is.na(cumhaz_t)] <- 0
  var_log_t[is.na(var_log_t)] <- 0

  # Linear predictor for each newdata row, using the fitted coefficient names.
  lp_info <- ar_lp_from_newdata(fit, newdata, beta)

  rows <- vector("list", nrow(newdata) * length(times))
  k <- 0L
  for (r in seq_len(nrow(newdata))) {
    lp_r <- lp_info$lp[r]
    # x_r' V_beta x_r — the covariate contribution to Var(log Lambda_x(t)).
    # Computed once per row, not per time (times share the same beta variance).
    x_r <- lp_info$mm[r, , drop = FALSE]
    var_beta_r <- as.numeric(x_r %*% vcov_beta %*% t(x_r))

    for (j in seq_along(times)) {
      # Lambda_x(t) = exp(beta'x) * Lambda_0(t)
      lambda_x <- exp(lp_r) * cumhaz_t[j]

      if (!is.finite(lambda_x) || lambda_x <= 0) {
        # cumhaz_t = 0 before the first event -> F = 0 exactly, no CI needed
        rows[[k <- k + 1L]] <- list(
          row = r,
          time = times[j],
          estimate = 0,
          ci_lower = 0,
          ci_upper = 0
        )
        next
      }

      # Delta method on log(-log(1 - F(t))) = log(Lambda_x(t)) = xi
      # Var(xi) = x' V_beta x  +  Var(log Lambda_0(t))
      se_xi <- sqrt(var_beta_r + var_log_t[j])
      ci <- cloglog_risk_ci(log(lambda_x), se_xi, z_crit)

      rows[[k <- k + 1L]] <- list(
        row = r,
        time = times[j],
        estimate = ci$estimate,
        ci_lower = ci$ci_lower,
        ci_upper = ci$ci_upper
      )
    }
  }

  new_matchatr_absolute_risk(
    estimates = data.table::rbindlist(rows),
    times = times,
    newdata = newdata,
    conf_level = conf_level,
    ci_method = "delta (log-log)",
    engine = fit$engine,
    estimator = fit$estimator,
    method = method
  )
}

#' Risk estimate and complementary-log-log CI from a cloglog linear predictor
#'
#' Shared inversion for every absolute-risk engine. Given ξ = log(−log S(t | x))
#' (the cumulative incidence on the complementary log-log scale) and its standard
#' error, returns the cumulative incidence F = 1 − exp(−exp(ξ)) and a pointwise CI
#' formed by inverting the symmetric Wald interval ξ ± z·SE(ξ) back to the risk
#' scale. Because higher ξ means lower survival, ξ − z·SE maps to the lower risk
#' bound. The CI bounds are clamped to `[0, 1]` to defend against numerical overflow
#' at extreme times. The Cox-type engines pass ξ = log Λ_x(t); the AFT engine
#' passes ξ = (log t − η)/σ̂ — the same scale, different ξ.
#'
#' @param xi Numeric scalar ξ = log(−log S(t | x)).
#' @param se_xi Numeric scalar standard error of ξ (delta method).
#' @param z_crit Numeric critical value (e.g. `qnorm(0.975)`).
#' @returns A list with `$estimate`, `$ci_lower`, `$ci_upper` on the probability
#'   scale.
#' @family contrasts
#' @noRd
cloglog_risk_ci <- function(xi, se_xi, z_crit) {
  f_est <- 1 - exp(-exp(xi))
  # Higher xi -> lower S -> higher F, so xi - z*se gives the lower risk bound.
  ci_lo <- 1 - exp(-exp(xi - z_crit * se_xi))
  ci_hi <- 1 - exp(-exp(xi + z_crit * se_xi))
  list(
    estimate = f_est,
    ci_lower = max(0, min(1, ci_lo)),
    ci_upper = max(0, min(1, ci_hi))
  )
}

#' Construct a matchatr_absolute_risk object
#'
#' Wraps the estimates table and evaluation metadata in the `matchatr_absolute_risk`
#' S3 structure shared by the case-cohort, IPW Cox, and IPW AFT engines, so the
#' class and field layout live in one place.
#'
#' @param estimates A `data.table` with columns `row`, `time`, `estimate`,
#'   `ci_lower`, `ci_upper`.
#' @param times Sorted numeric evaluation times.
#' @param newdata The covariate patterns supplied to [absolute_risk()].
#' @param conf_level Numeric confidence level.
#' @param ci_method Character label for the interval method.
#' @param engine,estimator Character labels carried from the fit.
#' @param method Character label for the baseline-hazard / parametric variant.
#' @returns A `matchatr_absolute_risk` object.
#' @family contrasts
#' @noRd
new_matchatr_absolute_risk <- function(
  estimates,
  times,
  newdata,
  conf_level,
  ci_method,
  engine,
  estimator,
  method
) {
  structure(
    list(
      estimates = estimates,
      times = times,
      newdata = newdata,
      conf_level = conf_level,
      ci_method = ci_method,
      engine = engine,
      estimator = estimator,
      method = method
    ),
    class = c("matchatr_absolute_risk", "matchatr")
  )
}

#' Prepend a t = 0 fence post to a Breslow step function
#'
#' Adds a `t = 0` knot (cumulative hazard 0, variance 0) so that evaluation times
#' before the first event map to `F = 0` rather than extrapolating to the first
#' event's value. The fence is skipped when an event already occurs at `t = 0`
#' (possible with rounded / discrete times): prepending another `0` would
#' duplicate the knot and make `stats::approx()` collapse the two, silently
#' dropping the real time-0 increment.
#'
#' @param t_events Sorted numeric vector of unique event times.
#' @param cumhaz Cumulative baseline hazard at each event time.
#' @param var_log Var(log Λ̂₀) at each event time.
#' @returns A list with `$times`, `$cumhaz`, `$var_log_cumhaz`, strictly
#'   increasing in `$times`.
#' @family contrasts
#' @noRd
breslow_step_with_fence <- function(t_events, cumhaz, var_log) {
  if (length(t_events) == 0L || t_events[1] > 0) {
    list(
      times = c(0, t_events),
      cumhaz = c(0, cumhaz),
      var_log_cumhaz = c(0, var_log)
    )
  } else {
    # An event already sits at t = 0; the series starts there, no fence needed.
    list(times = t_events, cumhaz = cumhaz, var_log_cumhaz = var_log)
  }
}

#' Linear predictor for covariate patterns from a survival fit
#'
#' Builds the model matrix for `data` from the **fitted model's own terms**,
#' selects the columns corresponding to the fitted coefficients, and returns both
#' the model matrix and the linear predictor. Used by all three absolute-risk
#' engines (case-cohort `cch`, IPW nested case-control weighted Cox `ipw_cox`, and
#' IPW nested case-control Weibull AFT `ipw_aft`), for both `newdata` patterns and
#' the analysis sample itself.
#'
#' @details
#' For the model-terms engines (`ipw_cox` coxph, `ipw_aft` survreg) the design is
#' built with `model.matrix(delete.response(terms(model)), model.frame(..., xlev =
#' model$xlevels))` rather than re-deriving the formula from `term.labels`. This
#' is essential for any **data-dependent** term — `poly()`, `ns()` / `bs()`,
#' `scale()`, `cut()` — whose basis depends on the rows it is computed over: the
#' fitted `terms` object carries the original basis in its `"predvars"`
#' attribute, so applying it to new rows reproduces the basis the model was
#' fitted with. Re-deriving the formula and calling `model.matrix` on the new
#' rows alone would recompute the basis from those rows, silently producing a
#' different design with the *same* coefficient names (so a name-based guard
#' would not catch it) and a wrong linear predictor.
#'
#' The case-cohort (`cch`) engine takes the original-formula path instead:
#' `survival::cch` rewrites its model formula into a non-standard internal form,
#' so its fitted `terms` object does not carry a `predvars`/contrasts map aligned
#' with the reported coefficient names. The design is rebuilt from the
#' user-facing `exposure + confounders` formula (whose standard R contrasts match
#' the `cch` coefficient names); data-dependent confounder transforms are not
#' reproduced for `cch`, which its designs do not use.
#'
#' @param fit A `matchatr_fit` with engine `"cch"`, `"ipw_cox"`, or `"ipw_aft"`.
#' @param data A data frame of covariate patterns (a `newdata` grid or the
#'   analysis sample). Must carry the columns the fitted model's terms reference.
#' @param beta Named numeric vector of fitted coefficients.
#' @returns A list with `$mm` (model matrix, p columns aligned to `names(beta)`)
#'   and `$lp` (numeric vector of length `nrow(data)`).
#' @family contrasts
#' @noRd
ar_lp_from_newdata <- function(fit, data, beta) {
  mm_new <- tryCatch(
    if (identical(fit$engine, "cch")) {
      # cch: rebuild from the user-facing formula (see @details).
      conf_terms <- if (is.null(fit$confounders)) {
        character(0)
      } else {
        attr(stats::terms(fit$confounders), "term.labels")
      }
      stats::model.matrix(
        stats::reformulate(c(fit$exposure, conf_terms)),
        data = data
      )
    } else {
      # coxph / survreg: reuse the fitted terms (predvars basis) and factor levels
      # so any data-dependent transform is reproduced, not recomputed.
      mt <- stats::delete.response(stats::terms(fit$model))
      stats::model.matrix(
        mt,
        stats::model.frame(mt, data = data, xlev = fit$model$xlevels)
      )
    },
    error = function(e) {
      rlang::abort(
        c(
          "Could not build a model matrix from the supplied covariate data.",
          i = conditionMessage(e)
        ),
        class = c("matchatr_bad_input", "matchatr_error")
      )
    }
  )
  # Guard: all coefficient names must appear in the model matrix
  missing_cols <- setdiff(names(beta), colnames(mm_new))
  if (length(missing_cols) > 0L) {
    rlang::abort(
      c(
        "Columns in the covariate data do not match the fitted model.",
        i = paste0(
          "Missing or mismatched columns: ",
          paste(missing_cols, collapse = ", ")
        )
      ),
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  mm_sel <- mm_new[, names(beta), drop = FALSE]
  list(mm = mm_sel, lp = as.vector(mm_sel %*% beta))
}

#' Print a matchatr_absolute_risk object
#'
#' Displays a compact summary of the absolute-risk result: the engine,
#' pseudo-likelihood / weighting method, CI method, the evaluation times, and the
#' number of covariate patterns, followed by the estimates table.
#'
#' @param x A `matchatr_absolute_risk` object returned by [absolute_risk()].
#' @param ... Unused; present for S3 consistency.
#' @returns Invisibly returns `x`.
#' @family contrasts
#' @seealso [absolute_risk()]
#' @export
print.matchatr_absolute_risk <- function(x, ...) {
  cat("<matchatr_absolute_risk>\n")
  cat(" Engine:     ", x$engine, "\n", sep = "")
  cat(" Method:     ", x$method, "\n", sep = "")
  cat(" CI method:  ", x$ci_method, "\n", sep = "")
  cat(
    " Times:      ",
    length(x$times),
    " evaluation time",
    if (length(x$times) == 1L) "" else "s",
    " (",
    paste(
      round(x$times[seq_len(min(5L, length(x$times)))], 2),
      collapse = ", "
    ),
    if (length(x$times) > 5L) ", ..." else "",
    ")\n",
    sep = ""
  )
  cat(
    " Patterns:   ",
    nrow(x$newdata),
    " covariate pattern",
    if (nrow(x$newdata) == 1L) "" else "s",
    "\n",
    sep = ""
  )
  cat("\nAbsolute risk estimates:\n")
  print(x$estimates)
  invisible(x)
}

#' Tidy a matchatr_absolute_risk object
#'
#' Returns the estimates table as a `data.table` with one row per (covariate
#' pattern row, evaluation time), with columns `row`, `time`, `estimate`,
#' `ci_lower`, `ci_upper`.
#'
#' @param x A `matchatr_absolute_risk` object.
#' @param ... Unused.
#' @returns A `data.table`.
#' @family contrasts
#' @seealso [absolute_risk()]
#' @export
tidy.matchatr_absolute_risk <- function(x, ...) {
  x$estimates
}
