#' Construct a `matchatr_fit` object
#'
#' Low-level constructor for the object returned by [matcha()]. It bundles the
#' validated analysis specification — data, outcome/exposure roles, the
#' sampling `design`, the chosen `estimator`, and the resolved `engine` key —
#' together with a `details` list that carries the design's weighting scheme
#' and reserves the slots later inference fills in. The `model` slot holds the
#' fitted estimation object; it is `NULL` until an estimation engine is run,
#' mirroring the deferred-fit pattern used elsewhere in the etverse.
#'
#' @param model Fitted estimation object, or `NULL` when no engine has run.
#' @param data A data.table of the analysis data.
#' @param outcome Character scalar naming the case-status column.
#' @param exposure Character scalar naming the exposure column.
#' @param confounders A one-sided formula of confounders, or `NULL`.
#' @param design The `matchatr_design` object describing the sampling
#'   structure.
#' @param estimator Character scalar estimator name.
#' @param engine Character scalar engine key the (design, estimator) pair
#'   resolved to.
#' @param details Named list of estimator/design metadata. Reserves
#'   `variance_kind`, `cc_weights`, and `design_weights` — the case-control
#'   weights and design (inclusion-probability) weights are kept in distinct
#'   slots because their variance consequences differ.
#' @param call The original [matcha()] call, for printing.
#' @returns A list with class `"matchatr_fit"`.
#' @family constructors
#' @noRd
new_matchatr_fit <- function(
  model,
  data,
  outcome,
  exposure,
  confounders,
  design,
  estimator,
  engine,
  details = list(),
  call = NULL
) {
  structure(
    list(
      model = model,
      data = data,
      outcome = outcome,
      exposure = exposure,
      confounders = confounders,
      design = design,
      estimator = estimator,
      engine = engine,
      details = details,
      call = call
    ),
    class = "matchatr_fit"
  )
}

#' Construct a `matchatr_result` object
#'
#' Low-level constructor for the object returned by the (reused) `contrast()`
#' verb. It holds the intervention-specific estimates and the pairwise causal
#' contrasts together with the variance-covariance matrix and the metadata a
#' print / tidy method needs. Populated by the causal-contrast layer; the
#' constructor exists here so the result class is defined alongside the fit.
#'
#' @param estimates A data.table of design / intervention-specific estimates.
#' @param contrasts A data.table of pairwise contrasts with SEs and CIs.
#' @param type Character contrast type (e.g. `"difference"`, `"ratio"`,
#'   `"or"`).
#' @param estimand Character estimand label.
#' @param ci_method Character CI method (e.g. `"sandwich"`, `"bootstrap"`).
#' @param reference Character name of the reference level, or `NULL`.
#' @param n Integer analysis sample size.
#' @param estimator Character estimator name.
#' @param engine Character engine key.
#' @param vcov Variance-covariance matrix of the estimates, or `NULL`.
#' @param call The original `contrast()` call, for printing.
#' @returns A list with class `"matchatr_result"`.
#' @family constructors
#' @noRd
new_matchatr_result <- function(
  estimates,
  contrasts,
  type,
  estimand,
  ci_method,
  reference = NULL,
  n,
  estimator,
  engine,
  vcov = NULL,
  call = NULL
) {
  structure(
    list(
      estimates = estimates,
      contrasts = contrasts,
      type = type,
      estimand = estimand,
      ci_method = ci_method,
      reference = reference,
      n = n,
      estimator = estimator,
      engine = engine,
      vcov = vcov,
      call = call
    ),
    class = "matchatr_result"
  )
}
