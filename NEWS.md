# matchatr (development version)

## 2026-06-03 — McNemar 1:1 matched-pair odds ratio (PHASE_3 Chunk 2)

Adds the closed-form 1:1 estimator alongside the conditional logistic engine.
`matcha(design = matched_cc(strata = ...), estimator = "mcnemar")` computes the
matched-pair odds ratio directly from the discordant-pair counts — `OR = n10/n01`
with `Var(log OR) = 1/n10 + 1/n01` (McNemar 1947; Breslow & Day 1980) — without
fitting `survival::clogit`, and `contrast(type = "or")` reports it with a Wald
interval. Pairs concordant on exposure carry no information and cancel.

- The estimator applies only to genuine 1:1 pairs: a matched set with more than
  one case or more than one control is M:1 (or richer) matching with no two-cell
  closed form and is rejected with `matchatr_not_one_to_one`, pointing to
  `estimator = "clogit"`. A one-sided (or empty) set of discordant pairs gives a
  boundary OR of 0 / ∞ and aborts with `matchatr_unestimable_exposure`.
- The exposure must be binary (logical / two-level factor / numeric 0/1); a
  non-binary exposure is declined (`matchatr_bad_input`) toward `clogit`. RD / RR
  remain unidentified and `ci_method = "sandwich"` / `"bootstrap"` are declined
  (`matchatr_unsupported_variance`); a missing pair member drops to complete
  pairs with a `matchatr_dropped_rows` warning.
- Validated three ways: exact agreement with `survival::clogit` on the same 1:1
  binary data (the conditional likelihood reduces to McNemar's), the OR and
  variance against the hand-counted closed form (independent of clogit), and a
  matched-pair DGP with a known log-OR (CMLE recovers β within a self-scaling SE
  band). A dedicated test pins the OR²-bias invariant: the unconditional 1:1 MLE
  with a parameter per pair is exactly twice the conditional estimate.

## 2026-06-03 — Matched case-control conditional logistic regression (PHASE_3 Chunk 1)

Opens the matched case-control layer with the conditional maximum-likelihood
odds ratio. `matcha(design = matched_cc(strata = ...), estimator = "clogit")`
(the design's default estimator) fits `survival::clogit` —
`outcome ~ exposure + confounders + strata(set)`, each matched set a stratum —
and `contrast(type = "or")` reports the exposure's conditional odds ratio with a
partial-likelihood-information Wald interval. Conditioning on the matched-set
totals removes the matching-variable nuisance parameters, so only the exposure /
adjustment ORs are reported; the matching variables are controlled implicitly.

- The conditional likelihood is the correctness invariant: unconditional logistic
  regression on matched-set indicators biases the OR (for 1:1 matching its MLE
  converges to the squared OR; Pike et al. 1980, Breslow & Day 1980).
- Several matching columns cross into one `strata()` term (frequency matching);
  a factor exposure reports one OR per level versus its reference; non-matching
  covariates adjust via `confounders`.
- RD / RR are rejected as unidentified (shared with the logistic engine). The
  conditional fit reports the information-matrix interval only, so
  `ci_method = "sandwich"` / `"bootstrap"` are declined
  (`matchatr_unsupported_variance`); cluster-robust variance for reused controls
  is deferred to the risk-set designs. An exposure with no within-stratum
  variation aborts with `matchatr_unestimable_exposure`.
- Validated against an independent closed-form oracle for 1:1 matching (McNemar:
  OR = n10/n01 and Var(log OR) = 1/n10 + 1/n01, exact and computed without
  clogit), a matched-set DGP built from the conditional likelihood with a known
  log-OR (CMLE recovers β for binary, continuous, mixed-ratio, and adjusted
  cases), `survival::clogit` pass-through, and a regression pin against the
  canonical `infert` clogit example (induced ≈ 4.09, spontaneous ≈ 7.29).
- The shared `conditional_or_result()` assembly (exposure coefficient by term
  position, Wald interval on the log scale, exponentiated) now backs both the
  unmatched logistic and the matched conditional logistic engines.

## 2026-06-03 — Mantel-Haenszel stratified odds ratio (PHASE_2 Chunk 3)

Completes the unmatched case-control layer with the closed-form Mantel-Haenszel
estimator. `unmatched_cc()` gains a `strata` argument (the stratifying
variable(s); several are crossed into one factor, and none gives the crude
single-table OR). `matcha(estimator = "mh")` computes the summary odds ratio
over the per-stratum 2×2 tables, and `contrast(type = "or")` reports it with a
**Robins-Breslow-Greenland** (1986) Wald confidence interval — valid in both the
sparse-data and large-strata limits.

- The exposure must be binary (a 2×2 table per stratum); a categorical (k>2) or
  continuous exposure is rejected (`matchatr_bad_input`) with a pointer to
  `estimator = "logistic"`.
- RD / RR are rejected as unidentified (shared with the logistic engine). The MH
  variance is the closed-form RBG estimator, so `ci_method = "sandwich"` /
  `"bootstrap"` are declined (`matchatr_unsupported_variance`); a zero
  exposure-outcome margin aborts with `matchatr_unestimable_exposure`.
- Validated against `stats::mantelhaen.test(correct = FALSE)`, whose odds-ratio
  estimate **and** confidence interval use the same RBG variance — the matchatr
  OR and CI match it exactly — plus the closed-form 2×2 odds ratio for the crude
  case.
- The fitter-agnostic coefficient-extraction helpers (`term_assign`,
  `estimable_vcov`, `exposure_coef_index`, `parametric_positions`) moved to
  `R/coef_extract.R`, shared by the logistic and Mantel-Haenszel engines.

## 2026-06-03 — PHASE_2 Chunk 2 critical-review fixes

Follow-up fixes from an adversarial review of the categorical/GAM exposure work.

- Exposure-coefficient extraction is now by term *position* (collision-free)
  rather than by reconstructed coefficient name. `glm` permits non-unique
  coefficient names from `factor x level` concatenation — e.g. exposure `ses`
  level `low` and confounder `se` level `slow` both yield `"seslow"` — and the
  previous name-based selection could grab the confounder's coefficient as a
  spurious second exposure OR, while name-based SE indexing returned the wrong
  (first-match) or `NA` standard error. Variance is now indexed by position via
  a unified `term_assign()` (glm `model.matrix` assign / gam `$assign` + `$pterms`)
  and a position-based `estimable_vcov()`, robust to both name collisions and
  aliasing.
- An ordered-factor exposure is now rejected up front in `matcha()`
  (`matchatr_bad_input`) instead of at `contrast()` time, so no model is fit and
  no spurious missing-data warning is emitted for an analysis that cannot yield a
  per-level OR.
- `tidy()` / `summary()` on a GAM fit now report only the parametric
  coefficients. Previously the smooth-basis terms (`s(age).1`, ...) were listed
  as rows and, with `exponentiate = TRUE`, their `exp()` was reported as an
  "odds ratio" with a Wald p-value — penalized basis weights are neither.
- `model_fn` is validated more strictly: it must accept a `family` argument (or
  `...`) — caught at `matcha()` with a classed `matchatr_bad_input` instead of a
  raw "unused argument" base error — and the fitted object must be a binomial
  fit. A fitter that ignores `family` and returns, say, an OLS `lm` now aborts
  with `matchatr_bad_model_fit` rather than silently exponentiating a
  linear-probability slope into a bogus "odds ratio".
- The recorded reference level for a factor exposure is now read from the
  fitted model's `xlevels` (the baseline actually used), not the factor's first
  *declared* level. A factor with an unused first level (e.g. declared
  `absent < med < high` but only `med`/`high` observed) previously mislabeled the
  contrast as "vs absent"; it now correctly reports "vs med".

## 2026-06-03 — Categorical / ordinal / GAM exposures (PHASE_2 Chunk 2)

Extends the unmatched case-control logistic engine to every exposure type.

- **Categorical (k>2) exposure**: a factor exposure now yields one odds ratio per
  non-reference level; `contrast()` returns k-1 rows and the result records the
  factor's reference level (`result$reference`).
- **Ordinal trend**: an ordinal exposure entered as a numeric score yields a
  single per-step trend OR.
- **Pluggable fitter**: `matcha(..., model_fn = )` selects the logistic fitter
  (default `stats::glm`); e.g. `model_fn = mgcv::gam` with `confounders = ~ s(age)`
  adjusts for a confounder with a smooth term while the exposure stays parametric
  with an interpretable OR. Exposure-coefficient extraction is now by name, so it
  works across `glm` and `gam`.
- **Ordered-factor exposure rejected** (`matchatr_bad_input`): `glm`/`gam` fit it
  with polynomial contrasts (`.L`, `.Q`, ...), which are not per-level ORs; the
  error points to a numeric score (trend) or an unordered factor (per-level).
- Validated against `stats::glm` for every exposure type, and against the
  Ille-et-Vilaine **`esoph`** case-control data (handbook Ch3): the categorical
  alcohol ORs reproduce the canonical monotone dose-response.

## 2026-06-02 — PHASE_2 Chunk 1 critical-review fixes

Follow-up fixes from an adversarial review of the unmatched case-control
logistic engine.

- `tidy()` / `summary()` with `robust = TRUE`, and `contrast(type = "or",
  ci_method = "sandwich")`, returned misaligned standard errors when the fitted
  model had an aliased (rank-deficient / collinear) term. `sandwich::sandwich()`
  drops aliased columns while `stats::coef()` keeps them, so the shorter SE
  vector was recycled against the coefficient vector and every term at or after
  the first aliased one got the wrong SE / CI / p-value. Standard errors are now
  aligned to the coefficient vector by name — aliased terms get `NA` — via a
  shared `estimable_vcov()` helper that reduces both variance sources to the
  estimable coefficient set.
- A constant or collinear exposure aliases to `NA` in `glm`; `contrast(type =
  "or")` previously returned an `NA` odds ratio silently. It now aborts with the
  classed `matchatr_unestimable_exposure`, mirroring the degenerate-outcome
  rejection.
- `matchatr_result$n` reported `nrow(data)`, overcounting when `glm` dropped rows
  with missing covariates. It now reports `stats::nobs(model)` — the complete-case
  count actually used, matching `causatr` — and a `matchatr_dropped_rows` warning
  reports how many rows were dropped. Actual missing-data handling (multiple
  imputation) stays delegated to `causatr::causat_mice` / the future imputation
  phase, not reimplemented here.
- `contrast()` now defaults `type` to the estimand the design identifies: `"or"`
  for the classical odds-ratio engines, `"difference"` otherwise. `contrast(fit)`
  on an unmatched case-control logistic fit therefore returns the conditional OR
  instead of erroring on the previously hard-coded risk-difference default.

## 2026-06-02 — Unmatched case-control logistic conditional OR (PHASE_2 Chunk 1)

First estimator engine: the classical unmatched case-control analysis. `matcha()`
now *runs* the resolved engine as part of the fit — the `"logistic"` estimator
fits `stats::glm(family = binomial)` for `outcome ~ exposure + confounders` and
stores it in the `model` slot (engines without a wired estimator still leave
`model = NULL`). Only the slope coefficients carry over from the cohort model;
the case-control intercept is offset by the log sampling-fraction ratio
(Prentice & Pyke, 1979) and is never reported as a baseline risk.

- `contrast(fit, type = "or")` returns the exposure's conditional odds ratio as a
  `matchatr_result`, with a Wald interval from the model information matrix
  (`ci_method = "model"`, the new default) or the Huber-White sandwich
  (`ci_method = "sandwich"`). `contrast()` gained a `conf_level` argument.
- The risk difference and risk ratio are **not identified** from an unmatched
  case-control sample without the source-population prevalence q0:
  `contrast(type = "difference" / "ratio")` aborts with the classed
  `matchatr_unidentified_estimand`, pointing to `type = "or"` or to a
  case-control-weighted estimator with `prevalence`. Bootstrap CIs for the
  conditional OR abort with `matchatr_unsupported_variance` (a Wald interval is
  reported instead).
- New S3 surface: `tidy.matchatr_fit()` (broom-style coefficient / OR table on the
  log-odds or, with `exponentiate = TRUE`, the OR scale; model-based or `robust`
  sandwich SE), `tidy.matchatr_result()`, `summary.matchatr_fit()`, and
  `print.matchatr_result()`.
- Tested against three oracles: `stats::glm` pass-through (point estimate, Wald and
  sandwich CI), the closed-form 2×2 odds ratio with Woolf variance, and a cohort
  data-generating process with a known log-OR (the conditional OR recovers the
  cohort slope, since case-control sampling shifts only the intercept).

## 2026-06-02 — PHASE_1 validation hardening (critical review)

Follow-up fixes from an adversarial review of the PHASE_1 foundation, closing
input-validation gaps before the estimation layers depend on them.

- `ratio = Inf` (or `-Inf`) now raises the classed `matchatr_bad_ratio` instead
  of an unclassed base-R error. `Inf %% 1` is `NaN` and `NaN != 0` is `NA`, so
  the old guard reached `if (NA)`; a finiteness check now runs before the modulo.
- A degenerate single-class outcome (all cases, all controls, or all `NA`) is
  now rejected with `matchatr_bad_outcome` for every encoding. The contrast
  check previously fired only for numeric columns, so an all-controls logical or
  two-level factor slipped through with `n_cases = 0`; `resolve_binary_outcome()`
  now coerces to 0/1 first and applies one uniform both-classes-present check.
- `matcha()` rejects `data` with duplicated column names (`matchatr_bad_input`).
  `[[` resolves a duplicated name to its first match, so a duplicated outcome /
  exposure / strata column would otherwise be silently chosen and the
  existence check could even blame the wrong column.
- `matcha()` rejects a column assigned to two incompatible roles
  (`matchatr_bad_input`): the outcome or exposure also appearing as a confounder
  or a design column (e.g. the exposure double-entered in `confounders`, or the
  outcome used as the matched-set id). Confounders and design columns may still
  overlap, so frequency-matching on a variable and adjusting for it remains valid.

## 2026-06-02 — Design taxonomy, data model, and two-step API

First implemented layer: the sampling-design objects and the `matcha()` fit
verb that the rest of the package builds on. No estimator runs yet — this phase
is plumbing, validation, and dispatch.

- Six design constructors returning a unified `matchatr_design` S3 object:
  `unmatched_cc()`, `matched_cc()`, `nested_cc()`, `case_cohort()`,
  `two_phase()`, and `counter_matched()`. Each carries its sampling structure
  (strata, time, matching ratio, prevalence q0, subcohort / phase columns) and a
  `weight_spec` declaring the intended weighting scheme. Structural arguments are
  validated at construction (q0 strictly in (0, 1); ratio a whole number >= 1;
  strata a non-empty character vector) with classed `matchatr_*` errors.
- `matcha(data, outcome, exposure, design, confounders, estimator)` validates
  the request, resolves the orthogonal `(design, estimator)` axes to an
  estimation engine via an internal dispatch table, and returns a
  `matchatr_fit`. The case-control-weighted family (`"ccw_gformula"`,
  `"ccw_ipw"`, `"ccw_aipw"`, `"ccw_tmle"`) routes on any design but requires a
  prevalence q0; the classical estimators are design-specific (`"logistic"` /
  `"mh"`, `"clogit"`, `"cch"`, ...). The `model` slot is `NULL` until an
  estimation engine is run.
- Weights are never a data column: the fit reserves distinct
  `details$cc_weights` and `details$design_weights` slots (case-control weights
  and inclusion-probability weights have different variance consequences) plus a
  `details$variance_kind` slot for the eventual sampling-variance correction.
- Classed rejection paths, all matching on a `matchatr_error` parent:
  `matchatr_bad_estimator` (unknown / design-incompatible estimator),
  `matchatr_bad_outcome` (non-binary case indicator),
  `matchatr_missing_prevalence` (CCW without q0),
  `matchatr_bad_prevalence` / `matchatr_bad_ratio` / `matchatr_bad_strata`
  (malformed construction), `matchatr_bad_design` (missing columns / wrong
  object), and a `matchatr_uninformative_stratum` warning when a conditional
  likelihood would drop a matched set with no case or no control.
- `print` methods for `matchatr_design` and `matchatr_fit`.
- `contrast()`, the second-step verb, is defined as a validated skeleton: it
  fixes the public signature (`type`, `ci_method`) and the `matchatr_result`
  return contract, and aborts with `matchatr_not_estimated` on a fit whose
  estimation engine has not run (`model = NULL`). The estimation body is filled
  in by the causal-contrast phases.

## 2026-06-02 — Package scaffold and design roadmap

Initial bootstrap of matchatr: causal inference for (matched) case-control,
nested case-control, and case-cohort study designs, as part of the etverse
ecosystem.

- Package scaffold: DESCRIPTION (Imports `causatr`, `survatr`, `survival`),
  MIT license, testthat 3 edition, Air formatting, `altdoc` + Quarto website
  (matching the other etverse packages), GitHub Actions (R-CMD-check,
  test-coverage / Codecov, format-check, altdoc), Makefile, and Claude Code
  configuration (`CLAUDE.md`, `.claude/hard-rules.md`, symlinked etverse skills,
  posit-dev skills bundle).
- Full design roadmap in `PHASE_1`–`PHASE_20` at the repository root, mapping the
  *Handbook of Statistical Methods for Case-Control Studies* (Borgan et al., 2018)
  to an implementation plan: design taxonomy and two-step API (PHASE_1); classical
  estimators (unmatched / matched / multiple-group, PHASE_2–4); time-to-event
  sampling designs (nested case-control, case-cohort, IPW-NCC, PHASE_5–7); the
  causal layer (strategy + case-control weighting g-formula / IPW / AIPW / TMLE +
  design-weighted causal survival, PHASE_8–10); and efficiency / advanced /
  extension phases (PHASE_11–20).

No estimator code is implemented yet — every phase is at `Status: DESIGN`.
