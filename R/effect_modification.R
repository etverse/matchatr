#' Assemble stratum-specific odds ratios from a clogit interaction fit
#'
#' Turns a matched case-control conditional logistic fit that carries an
#' `exposure x modifier` interaction into a `matchatr_result` reporting the
#' exposure's conditional odds ratio within each modifier level — the
#' stratum-specific OR. The exposure log OR at the modifier's reference level is
#' the exposure main coefficient `beta_x`; at each non-reference level it is the
#' linear combination `beta_x + beta_{x:level}`. Each combination's variance is
#' read from the joint partial-likelihood variance, so the per-level Wald
#' intervals account for the covariance between the main and interaction terms.
#'
#' @details
#' For a binary exposure and a modifier that is constant within each matched set
#' (the canonical case where the modifier is itself a matching variable), the
#' matched sets split into disjoint groups by modifier level and the conditional
#' likelihood factorises across them. The per-level estimate and variance then
#' reduce exactly to the 1:1 McNemar closed form on that level's discordant
#' pairs (OR = n10/n01, Var(log OR) = 1/n10 + 1/n01), which is the independent
#' oracle for both the point estimate and the variance.
#'
#' The interval is Wald on the log-odds scale and exponentiated, so it is
#' asymmetric on the OR scale; the OR-scale `se` is the delta-method value
#' OR * SE(log OR), kept for reference, while the reconstructable log-scale
#' estimate and SE live in the result's `estimates`. The exposure must
#' contribute a single coefficient (binary, continuous, or two-level factor);
#' a multi-level factor exposure is rejected upstream and again here as a
#' defensive backstop. An aliased (non-identified) main or interaction
#' coefficient aborts with `matchatr_unestimable_exposure` rather than returning
#' a silent `NA` OR.
#'
#' @param fit A `matchatr_fit` whose `model` is a fitted `survival::clogit` built
#'   with the `exposure * effect_modifier` interaction.
#' @param model The fitted `survival::clogit` object.
#' @param conf_level Numeric confidence level in (0, 1).
#' @param ci_method Character variance source recorded on the result (only
#'   `"model"` reaches here; the caller rejects the others).
#' @param n Integer analysis sample size (the clogit row count, `model$n`).
#' @param call Caller environment surfaced in any error.
#' @returns A `matchatr_result` with one stratum-specific odds-ratio row per
#'   modifier level, the reference level recorded in `reference`, and the
#'   per-level log-OR variance-covariance matrix in `vcov`.
#' @family estimators
#' @seealso [contrast()], `fit_clogit()`, `conditional_or_result()`
#' @noRd
stratum_specific_or_result <- function(
  fit,
  model,
  conf_level,
  ci_method,
  n,
  call = rlang::caller_env()
) {
  exposure <- fit$exposure
  em <- fit$effect_modifier
  beta <- stats::coef(model)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)

  # The exposure main-effect coefficient. A single-coefficient exposure is
  # required (enforced at the entry point); guard again so a multi-coefficient
  # exposure cannot silently mis-pair with the interaction columns.
  main_idx <- exposure_coef_index(model, exposure, call = call)
  int_idx <- interaction_coef_index(model, exposure, em, call = call)
  if (length(main_idx) != 1L) {
    rlang::abort(
      c(
        paste0(
          "Effect modification needs a single-coefficient exposure, but `",
          exposure,
          "` contributes ",
          length(main_idx),
          " coefficients."
        ),
        i = "Use a binary, continuous, or two-level factor exposure."
      ),
      class = c("matchatr_unsupported_combination", "matchatr_error"),
      call = call
    )
  }

  # Modifier levels in factor order (reference first). The interaction
  # coefficients map one-to-one to the non-reference levels in that same order
  # (a single-coefficient exposure contributes exactly one interaction column
  # per non-reference modifier level).
  levs <- model$xlevels[[em]]
  if (is.null(levs) || length(int_idx) != length(levs) - 1L) {
    rlang::abort(
      c(
        paste0(
          "The `",
          exposure,
          " * ",
          em,
          "` interaction does not map cleanly to the modifier levels."
        ),
        i = "Ensure the effect modifier is a categorical column with at least two observed levels."
      ),
      class = c("matchatr_unsupported_combination", "matchatr_error"),
      call = call
    )
  }

  # Each modifier level's exposure log OR is a linear combination of the
  # coefficients: beta_x at the reference level, beta_x + beta_{x:level}
  # elsewhere. An aliased (NA) main or interaction coefficient means the per-
  # level OR is not identified -- refuse rather than emit a silent NA.
  level_pos <- c(
    list(main_idx),
    lapply(int_idx, function(j) c(main_idx, j))
  )
  if (anyNA(beta[c(main_idx, int_idx)])) {
    rlang::abort(
      c(
        paste0(
          "The stratum-specific odds ratio for `",
          exposure,
          "` is not estimable."
        ),
        i = paste0(
          "The exposure main effect or an `",
          exposure,
          ":",
          em,
          "` interaction is collinear with the strata / confounders."
        )
      ),
      class = c("matchatr_unestimable_exposure", "matchatr_error"),
      call = call
    )
  }

  # Variance over the estimable coefficients (model partial-likelihood
  # information). A contrast matrix C with a 1 on each coefficient summed for a
  # level turns the joint variance into the per-level log-OR covariance
  # C V C'; its diagonal is each level's Var(log OR), accounting for the
  # main-interaction covariance.
  ev <- estimable_vcov(model, robust = FALSE)
  comparison <- paste0(exposure, " | ", em, " = ", levs)
  contrast_mat <- matrix(
    0,
    nrow = length(levs),
    ncol = length(ev$est_pos),
    dimnames = list(comparison, NULL)
  )
  for (i in seq_along(level_pos)) {
    contrast_mat[i, match(level_pos[[i]], ev$est_pos)] <- 1
  }
  log_or <- as.numeric(contrast_mat %*% beta[ev$est_pos])
  vcov_levels <- contrast_mat %*% ev$vcov %*% t(contrast_mat)
  dimnames(vcov_levels) <- list(comparison, comparison)
  se <- sqrt(diag(vcov_levels))

  log_lower <- log_or - z * se
  log_upper <- log_or + z * se

  estimates <- data.table::data.table(
    term = comparison,
    estimate = log_or, # log OR per modifier level
    se = se,
    ci_lower = log_lower,
    ci_upper = log_upper
  )
  contrasts <- data.table::data.table(
    comparison = comparison,
    estimate = exp(log_or), # OR per modifier level
    se = exp(log_or) * se, # delta-method SE on the OR scale
    ci_lower = exp(log_lower),
    ci_upper = exp(log_upper)
  )

  new_matchatr_result(
    estimates = estimates,
    contrasts = contrasts,
    type = "or",
    estimand = "stratum-specific conditional OR",
    ci_method = ci_method,
    reference = levs[1],
    n = n,
    estimator = fit$estimator,
    engine = fit$engine,
    vcov = vcov_levels,
    call = call
  )
}

#' Locate the coefficient(s) of an exposure-by-modifier interaction term
#'
#' Maps the `exposure:modifier` interaction to its coefficient position(s) by
#' term *position* via the parametric term assignment (`term_assign()`), so it
#' is collision-free even when `glm`-style fitters produce non-unique
#' coefficient names. The interaction term is the one whose component variables
#' are exactly the exposure and the modifier (order-independent), so it is found
#' regardless of which side `terms()` lists first.
#'
#' @param model A fitted model carrying an `exposure * modifier` interaction.
#' @param exposure Character scalar exposure column name.
#' @param modifier Character scalar effect-modifier column name.
#' @param call Caller environment surfaced in the error.
#' @returns An integer vector of coefficient positions (in `coef(model)`) for the
#'   interaction term, ordered to match the modifier's non-reference levels;
#'   aborts with `matchatr_bad_input` if the interaction term is absent.
#' @family estimators
#' @noRd
interaction_coef_index <- function(
  model,
  exposure,
  modifier,
  call = rlang::caller_env()
) {
  ta <- term_assign(model)
  # A two-way interaction label is "a:b"; split on ":" and compare the variable
  # set so "exposure:modifier" and "modifier:exposure" both match.
  parts <- strsplit(ta$labels, ":", fixed = TRUE)
  is_int <- vapply(
    parts,
    function(v) length(v) == 2L && setequal(v, c(exposure, modifier)),
    logical(1)
  )
  pos <- which(is_int)
  if (length(pos) != 1L) {
    rlang::abort(
      paste0(
        "No `",
        exposure,
        " : ",
        modifier,
        "` interaction term is present in the fitted model."
      ),
      class = c("matchatr_bad_input", "matchatr_error"),
      call = call
    )
  }
  which(ta$assign == pos)
}
