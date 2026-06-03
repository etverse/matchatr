#' Fit a case-control-type design
#'
#' @description
#' The fit verb for matchatr (mirroring `causatr::causat()` and
#' `survatr::surv_fit()`). `matcha()` takes the analysis data, the
#' outcome/exposure roles, a sampling `design` object, and an `estimator`, then
#' validates the request and resolves it to an estimation engine. The two
#' arguments are orthogonal: `design` selects the sampling structure (strata,
#' time, prevalence q0, inclusion weights) and `estimator` selects the analysis
#' (conditional vs marginal; odds ratio vs hazard ratio vs risk difference).
#'
#' @details
#' Weights are never read from or written to `data`. The design's `weight_spec`
#' records the intended scheme; the case-control weights (q0-based, Rose &
#' van der Laan) and design / inclusion-probability weights (Samuelsen, Borgan)
#' are kept in distinct slots on the fit (`details$cc_weights`,
#' `details$design_weights`) because their variance consequences differ.
#'
#' The resolved engine is run as part of the fit: an implemented estimator (the
#' unmatched case-control logistic regression) populates the `model` slot, while
#' an engine with no wired estimator leaves it `NULL`. `details$engine` records
#' the engine the (design, estimator) pair resolved to.
#'
#' @param data A data.frame or data.table. Not mutated; a data.table copy is
#'   stored on the fit.
#' @param outcome Character scalar naming the binary case-status column (a
#'   logical, two-level factor, or numeric 0/1 column).
#' @param exposure Character scalar naming the exposure column.
#' @param design A `matchatr_design` object from one of the design
#'   constructors ([unmatched_cc()], [matched_cc()], [nested_cc()],
#'   [case_cohort()], [two_phase()], [counter_matched()]).
#' @param confounders A one-sided formula of confounders (e.g. `~ age + smoke`),
#'   or `NULL` for an unadjusted analysis.
#' @param estimator Character scalar naming the analysis, or `NULL` to use the
#'   design's canonical default. Classical choices are design-specific
#'   (`"logistic"` / `"mh"` for unmatched CC, `"clogit"` for matched CC / NCC,
#'   `"cch"` for case-cohort); the case-control-weighted causal estimators
#'   `"ccw_gformula"`, `"ccw_ipw"`, `"ccw_aipw"`, `"ccw_tmle"` apply to any
#'   design but require a prevalence q0 on the design.
#' @param model_fn Optional model-fitting function for the unmatched
#'   case-control logistic engine, with a `(formula, family, data)` interface.
#'   Defaults to [stats::glm()]; pass e.g. `mgcv::gam` to adjust for a confounder
#'   with a smooth term (`confounders = ~ s(age)`) while keeping the exposure
#'   parametric. Ignored by the other engines.
#' @param effect_modifier `NULL` or a character scalar naming a categorical
#'   (logical / character / factor) column whose levels modify the exposure
#'   effect. When supplied, the conditional logistic engine fits
#'   `outcome ~ exposure * effect_modifier + confounders + strata(set)` and
#'   `contrast(type = "or")` reports the stratum-specific odds ratio of the
#'   exposure within each modifier level (one OR per level, with a Wald interval
#'   from the joint partial-likelihood variance). Supported only for
#'   `estimator = "clogit"` with a single-coefficient exposure (binary,
#'   continuous, or two-level factor); the modifier may coincide with a matching
#'   variable. Defaults to `NULL` (no effect modification).
#'
#' @returns A `matchatr_fit` object: a list with the validated specification
#'   (`data`, `outcome`, `exposure`, `confounders`, `design`, `estimator`,
#'   `engine`, `effect_modifier`), a `details` list (resolved engine, weighting
#'   scheme, reserved
#'   variance / weight slots, case and control counts), and the originating
#'   `call`. The `model` slot holds the fitted estimation object for an
#'   implemented engine, or `NULL` otherwise.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   case  = rep(c(1, 0), each = 50),
#'   x     = rbinom(100, 1, 0.4),
#'   age   = rnorm(100, 50, 10),
#'   smoke = rbinom(100, 1, 0.3),
#'   set   = rep(1:50, times = 2)
#' )
#'
#' # Unmatched case-control -> conditional odds ratio (logistic)
#' matcha(df, outcome = "case", exposure = "x",
#'        design = unmatched_cc(), confounders = ~ age + smoke)
#'
#' # Matched case-control -> conditional logistic regression
#' matcha(df, outcome = "case", exposure = "x",
#'        design = matched_cc(strata = "set"), estimator = "clogit")
#'
#' @seealso [unmatched_cc()], [matched_cc()], [nested_cc()], [case_cohort()]
#' @export
matcha <- function(
  data,
  outcome,
  exposure,
  design,
  confounders = NULL,
  estimator = NULL,
  model_fn = NULL,
  effect_modifier = NULL
) {
  # Record the user's call so the fit can echo it in print().
  call <- match.call()

  # Structural input validation, cheapest checks first.
  if (!is.data.frame(data)) {
    rlang::abort(
      "`data` must be a data.frame or data.table.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }
  if (!inherits(design, "matchatr_design")) {
    rlang::abort(
      c(
        "`design` must be a `matchatr_design` object.",
        i = "Build one with e.g. `unmatched_cc()`, `matched_cc()`, or `nested_cc()`."
      ),
      class = c("matchatr_bad_design", "matchatr_error")
    )
  }
  check_string(outcome)
  check_string(exposure)
  if (!is.null(confounders)) {
    check_formula(confounders)
  }
  # `model_fn` is the pluggable logistic fitter (e.g. `stats::glm`, `mgcv::gam`);
  # it must be a function that accepts a `family` argument (it is always called
  # as `model_fn(formula, family = binomial(), data = )`).
  if (!is.null(model_fn)) {
    if (!is.function(model_fn)) {
      rlang::abort(
        "`model_fn` must be a model-fitting function (e.g. `stats::glm`) or NULL.",
        class = c("matchatr_bad_input", "matchatr_error")
      )
    }
    fmls <- tryCatch(names(formals(args(model_fn))), error = function(e) NULL)
    if (length(fmls) > 0L && !any(c("family", "...") %in% fmls)) {
      rlang::abort(
        c(
          "`model_fn` must accept a `family` argument (or `...`).",
          i = "It is called as `model_fn(formula, family = binomial(), data = )`; use e.g. `stats::glm` or `mgcv::gam`."
        ),
        class = c("matchatr_bad_input", "matchatr_error")
      )
    }
  }
  if (identical(outcome, exposure)) {
    rlang::abort(
      "`outcome` and `exposure` must be different columns.",
      class = c("matchatr_bad_input", "matchatr_error")
    )
  }

  # Work on a data.table copy: user-facing functions return data.table and the
  # caller's frame must never be mutated.
  dt <- data.table::as.data.table(data)

  # Duplicated column names make `[[` ambiguous (it returns the first match), so
  # reject them before any column is looked up by name.
  check_unique_colnames(dt)

  # Every named column must exist before we touch its values.
  check_cols_exist(dt, outcome, arg = "outcome")
  check_cols_exist(dt, exposure, arg = "exposure")
  # matchatr reports per-level / single-variable effects; an ordered-factor
  # exposure (polynomial contrasts) is rejected up front.
  check_exposure_not_ordered(dt, exposure)
  if (!is.null(confounders)) {
    # all.vars() strips transforms (`I(age^2)` -> `age`) to plain names.
    check_cols_exist(dt, all.vars(confounders), arg = "confounders")
  }
  design_cols <- design_columns(design)
  if (length(design_cols) > 0L) {
    check_cols_exist(dt, design_cols, arg = "design")
  }

  # No column may serve two incompatible roles (outcome / exposure vs covariate
  # / design). Confounders and design columns are allowed to overlap.
  confounder_vars <- if (is.null(confounders)) {
    character(0)
  } else {
    all.vars(confounders)
  }
  check_role_collisions(outcome, exposure, confounder_vars, design_cols)

  # Resolve the analysis: design default when unspecified, then map the
  # (design, estimator) pair to an engine (rejects unknown estimators).
  if (is.null(estimator)) {
    estimator <- default_estimator(design$type)
  }
  check_string(estimator)
  routing <- resolve_engine(design$type, estimator)

  # Effect modification is the `exposure x modifier` interaction in the
  # conditional logistic model: validate the modifier (categorical, distinct
  # from outcome/exposure, single-coefficient exposure) and restrict it to the
  # conditional logistic engine, where the matched-set conditioning makes the
  # per-level odds ratio the design-faithful contrast.
  if (!is.null(effect_modifier)) {
    check_effect_modifier(dt, effect_modifier, outcome, exposure)
    if (!identical(routing$engine, "clogit")) {
      rlang::abort(
        c(
          "`effect_modifier` is supported only for the conditional logistic estimator.",
          i = paste0(
            "Use `estimator = \"clogit\"` with a matched (or nested) ",
            "case-control design; got engine `",
            routing$engine,
            "`."
          )
        ),
        class = c("matchatr_bad_input", "matchatr_error")
      )
    }
  }

  # Case-control weighting reweights the sample to the source population, which
  # is impossible without the marginal prevalence q0.
  if (identical(routing$kind, "ccw") && is.null(design$prevalence)) {
    rlang::abort(
      c(
        paste0(
          "Estimator `",
          estimator,
          "` needs the source-population prevalence q0 to reweight the sample."
        ),
        i = "Supply it on the design, e.g. `unmatched_cc(prevalence = 0.02)`."
      ),
      class = c("matchatr_missing_prevalence", "matchatr_error")
    )
  }

  # Every design in this layer carries a binary case indicator; resolving it
  # both validates the outcome and yields the 0/1 vector for the strata checks.
  y01 <- resolve_binary_outcome(dt, outcome)

  # Conditional partial likelihood drops sets without both a case and a
  # control; warn so a mis-sampled design is not silently thinned.
  if (isTRUE(routing$conditional) && !is.null(design$strata)) {
    strata_list <- lapply(design$strata, function(col) dt[[col]])
    warn_uninformative_strata(strata_list, y01)
  }

  details <- list(
    engine = routing$engine,
    kind = routing$kind,
    conditional = routing$conditional,
    # The design declares how weights are to be computed; the realised weight
    # vectors are filled in by the weighting layer into the distinct slots
    # below. `variance_kind` records which sampling-variance correction applies.
    weight_spec = design$weight_spec,
    cc_weights = NULL,
    design_weights = NULL,
    variance_kind = NULL,
    # Pluggable logistic fitter (NULL -> stats::glm), used by the glm_logistic
    # engine; e.g. mgcv::gam for smooth confounder adjustment.
    model_fn = model_fn,
    n_cases = sum(y01 == 1L, na.rm = TRUE),
    n_controls = sum(y01 == 0L, na.rm = TRUE)
  )

  fit <- new_matchatr_fit(
    model = NULL,
    data = dt,
    outcome = outcome,
    exposure = exposure,
    confounders = confounders,
    design = design,
    estimator = estimator,
    engine = routing$engine,
    effect_modifier = effect_modifier,
    details = details,
    call = call
  )

  # Run the resolved engine: an implemented estimator populates `model`, while
  # an engine with no wired estimator leaves it NULL.
  fit$model <- run_engine(fit)
  fit
}
