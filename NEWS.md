# matchatr (development version)

## 2026-06-08 — Counter-matched NCC sampler and weighted partial likelihood (PHASE_5 Chunk 3)

Adds `sample_ncc_counter_matched()` and the `weighted_cox` engine for the
Langholz-Borgan (1995) counter-matched partial likelihood.

- **`sample_ncc_counter_matched(cohort, time, event, surrogate, m, match, entry)`**
  generates a counter-matched NCC dataset: at each event time the case is matched
  to `m` controls drawn from the *opposite* surrogate stratum. The output appends
  `set`, `case`, `risk_time`, and `log_w` (the log-sampling-weight for the
  weighted partial likelihood). Optional population-stratum matching (`match =`) and
  delayed entry (`entry =`) are supported, mirroring `sample_ncc()`. A case with no
  eligible control in the opposite stratum aborts with `matchatr_empty_risk_set`.
- **`matcha(design = counter_matched(strata, time, weights = "log_w"), estimator =
  "weighted_cox")`** fits `survival::coxph` with the log-weights as a Cox offset
  (`outcome ~ exposure + confounders + strata(set) + offset(log_w)`) and
  `contrast()` reports the hazard ratio (`type = "hr"`). `survival::coxph` is used
  directly (not the `clogit()` wrapper) because `clogit` does not pass `offset`
  through to its internal coxph call.
- The log-weight formula: the case represents its entire same-stratum risk set
  (log_w = log(n_same + 1)); each sampled control represents the opposite stratum
  divided by the controls drawn (log_w = log(n_opp / m_take)). This is the key
  correction that makes the weighted partial likelihood consistent for the Cox HR,
  while the unweighted clogit on counter-matched data is biased.
- `type = "or"`, `"difference"`, `"ratio"` and `ci_method = "sandwich"`,
  `"bootstrap"` are rejected for this estimator.
- New functions: `sample_ncc_counter_matched()`, `resolve_surrogate()` (in
  `R/risk_set_sampling.R`); `fit_weighted_cox()`, `contrast_weighted_cox()` (in
  `R/weighted_cox.R`). Wired in `R/dispatch.R` and `R/contrast.R`.
- Tests in `tests/testthat/test-weighted_cox.R`: structural invariants, exact
  log-weight formula check, truth recovery within 3.5 SE, full-cohort coxph
  oracle, and rejection paths.

## 2026-06-08 — Python cross-language oracles for the classical estimators

Adds an independent **Python (`statsmodels`) cross-language oracle** for every
implemented classical estimator, so a bug shared between matchatr and its R
engine cannot hide behind a same-package comparison. Each oracle reads the same
committed dataset both languages share and compares against a committed Python
result CSV; tests never invoke Python, so CI needs no Python toolchain (each is
guarded with `skip_if(!file.exists(...))`). The fixtures and a regeneration
recipe live in `tests/testthat/fixtures/python/`; the comparisons are in
`test-python-oracle.R`.

- Covered: the unmatched logistic OR (vs `Logit`), the matched conditional OR and
  the nested case-control HR (vs `ConditionalLogit`), the Mantel–Haenszel summary
  OR with the Robins–Breslow–Greenland interval (vs `StratifiedTable`), the
  polytomous subtype ORs (vs `MNLogit`), and the homogeneity Wald χ² + GLS pooled
  common OR (hand-built from `MNLogit`'s coefficients and covariance).
- `statsmodels` anchors these maximum-likelihood estimators; `delicatessen`
  (M-estimation + sandwich) is reserved for the causal / sandwich estimands of
  the later case-control-weighting phases.
- Test-only change: no user-facing behaviour is affected.

## 2026-06-08 — Risk-set control sampling: `sample_ncc()` (PHASE_5 Chunk 2)

Adds the exported `sample_ncc()`, which generates a nested case-control dataset
from a cohort by risk-set (incidence-density) control sampling. Each event
anchors a matched set holding the case and `m` controls drawn without replacement
from the subjects at risk at that failure time; the result is an analysis-ready
`data.table` (the cohort columns plus `set`, the per-set `case` indicator, and
`risk_time`) that feeds straight into
`matcha(design = nested_cc(strata = "set", time = "risk_time"))`. This closes the
generate-then-analyse loop opened by Chunk 1.

- **`sample_ncc(cohort, time, event, m = 1, match = NULL, entry = NULL)`.**
  `match = ~ s1 + s2` confines each case's controls to its own population
  stratum; `entry =` supplies a delayed-entry (left-truncation) column so a
  subject counts as at risk only after entering follow-up. Sampling uses the
  ambient random stream (wrap in `withr::with_seed()` / precede with
  `set.seed()`), and the input cohort is never mutated.
- **Native sampler, `Epi::ccwc` as an oracle.** Risk-set sampling is implemented
  natively (base R + data.table) so it is always available and deterministically
  seedable — consistent with the package's pattern of delegating only to
  Imports-tier *estimation* engines and hand-rolling sampling / closed forms.
  `Epi::ccwc` is an external cross-check in the tests (every sampled control must
  lie in the eligible risk-set pool), not a runtime dependency.
- **A case with no eligible control aborts with `matchatr_empty_risk_set`** — in
  the generation path, producing such a set is a sampling failure (a misspecified
  time origin/scale, an entry/exit mismatch, or over-fine `match` strata), unlike
  an uninformative analysis stratum, which `clogit` merely drops with a warning.
  A late failure time with fewer than `m` eligible controls keeps all of them (a
  smaller set), which is correct, not an error.
- **`match` strata are crossed with collision-proof `interaction()` codes**, and a
  missing value in a `match` column is rejected up front (`matchatr_bad_input`):
  an undefined stratum must not silently merge with the literal string `"NA"` nor
  match anyone.
- Validated by structural invariants (one case per set, controls genuinely at
  risk, no within-set reuse), the `Epi::ccwc` risk-set-definition cross-check, a
  cohort DGP whose known Cox log-HR the sampled subsample recovers within 3.5 SE
  and agrees with the full-cohort `survival::coxph` fit, and the classical
  `(m+1)/m` null efficiency. The test-only `sample_ncc_riskset()` fixture now
  delegates to `sample_ncc()`, so the Chunk 1 analysis tests also exercise the
  exported sampler.

## 2026-06-08 — Nested case-control hazard ratio (PHASE_5 Chunk 1)

Opens the time-to-event sampling designs with the classical nested case-control
(NCC) analysis. `matcha(design = nested_cc(strata = ..., time = ...), estimator =
"clogit")` now fits the risk-set conditional partial likelihood and `contrast()`
reports the exposure's **hazard ratio** (`type = "hr"`). The NCC analysis reuses
the matched-design `clogit` engine unchanged — a sampled risk set and a matched
set are the same stratum construction (`outcome ~ exposure + confounders +
strata(set)`) — so the value is the design-faithful estimand, not new estimation
machinery.

- **OR = HR exactly under risk-set (incidence-density) sampling** (Prentice &
  Breslow 1978), with no rare-disease assumption. The conditional partial
  likelihood gives `exp(beta)`; the sampling design fixes its meaning, so a
  matched case-control design reports the conditional odds ratio and a nested
  case-control design the hazard ratio. The shared `conditional_or_result()` /
  `stratum_specific_or_result()` assemblies gained a `type` argument that carries
  the scale label (`"or"` / `"hr"`); the arithmetic is identical.
- **A new contrast scale `type = "hr"`** (hazard ratio). `contrast()` defaults to
  it for the nested design (`default_contrast_type()` is now design-aware) and to
  `"or"` for the matched design. Each conditional design identifies exactly one
  scale: requesting an odds ratio from a risk-set design, or a hazard ratio from a
  matched design, names an estimand the design does not target and aborts with
  `matchatr_unidentified_estimand`. `summary()` and `print()` label the NCC table
  / contrast as a hazard ratio.
- The design's `time` column records how controls were sampled; the risk-set
  membership is read from `strata`, so the conditional likelihood does not enter
  `time` here (it feeds the later inclusion-weight / weighted-Cox designs). Risk
  set reuse — and hence a cluster-robust variance — belongs to those designs, so
  `ci_method = "sandwich"` / `"bootstrap"` stay unsupported; a risk set with no
  control is an uninformative stratum (`matchatr_uninformative_stratum`, dropped
  by `clogit`), consistent with the matched design.
- Validated against a cohort DGP with a known Cox log-HR (`make_ncc_cohort()`)
  from which an incidence-density NCC sample is drawn (`sample_ncc_riskset()`):
  the conditional likelihood recovers the cohort β within 3.5 SE and agrees with
  the full-cohort `survival::coxph` β (the OR = HR check); exact pass-through
  against a hand-fit `survival::clogit`; and the relative efficiency at the null
  matches the classical `(m+1)/m` variance ratio (Goldstein & Langholz 1992).
  Non-binary outcomes are rejected (`matchatr_bad_outcome`) and a within-set
  constant exposure is `matchatr_unestimable_exposure`.

## 2026-06-08 — Homogeneity test + pooled common OR for disease subtypes (PHASE_4 Chunk 2)

Completes the multiple-case/control-group layer with `test_homogeneity()`, which
answers the etiologic-heterogeneity question the polytomous design is for: does
the exposure act the same way on every disease subtype? Given an unconstrained
polytomous fit (`estimator = "polytomous"`), for each exposure term it tests the
null that the exposure odds ratio is constant across the non-reference subtypes
(H0: beta_1 = ... = beta_M) and reports the efficient pooled ("common") odds
ratio that holds under homogeneity. A binary or continuous exposure yields one
test; an unordered factor exposure one test per level. The new `R/homogeneity.R`
(`test_homogeneity()`, `homogeneity_one_term()`, `print` / `tidy` methods for the
`matchatr_homogeneity` class) carries it; `new_matchatr_homogeneity()` joins the
constructors.

- The homogeneity test is the canonical **Wald** test of etiologic heterogeneity
  (Begg & Gray 1984; `riskclustr::eh_test_subtype`):
  `W = (C b)' (C V C')^-1 (C b) ~ chi-squared(M - 1)`, where `b` is the stacked
  subtype log odds ratios and `V` their multinomial-information covariance
  (reused from `multinom_exposure_or()`, now also returning the structured
  subtype / predictor pieces). The common odds ratio is the minimum-variance
  (GLS / inverse-variance) restricted estimator `(1' V^-1 b)/(1' V^-1 1)`,
  asymptotically equal to the constrained maximum-likelihood fit. Both are
  computed on the **already-fitted** unconstrained model — no constrained refit —
  so there is no new dependency and continuous confounders are handled directly
  (chosen over a Poisson-surrogate or `VGAM` likelihood-ratio test, which would
  need discrete covariate patterns or a heavy new dependency). It mirrors the
  `C V C'` contrast construction already used for stratum-specific odds ratios.
- A non-polytomous fit (`logistic` / `mh` / `clogit`) or a non-`matchatr_fit`
  object is `matchatr_bad_input`; a fit whose engine produced no model is
  `matchatr_not_estimated`; a malformed `conf_level` is `matchatr_bad_input`.
- Validated against three oracles: the saturated 3-group / binary-exposure case
  reproduces the closed-form 2x2 Woolf chi-squared **and** pooled odds ratio
  (independent of `multinom`'s `vcov()`); the chi-squared / pooled OR equal the
  exact `C V C'` / GLS functionals of `contrast()` on an adjusted
  (continuous-confounder) fit; a truth DGP gives correct size under equal subtype
  ORs and power under unequal ones, with the pooled SE smaller than each subtype
  SE (Begg & Gray efficiency); and the heterogeneity p-value cross-checks against
  `riskclustr::eh_test_subtype` (the `mlogit` engine), which joins `Suggests`.

## 2026-06-03 — Polytomous logistic for multiple case / control groups (PHASE_4 Chunk 1)

Opens the multiple-case/control-group layer with unconstrained polytomous
(multinomial) logistic regression. `matcha(design = unmatched_cc(), estimator =
"polytomous")` fits a baseline-category multinomial model
(`outcome ~ exposure + confounders`) via `nnet::multinom` for an outcome with
three or more groups — multiple disease subtypes, or several control groups. The
new `reference =` argument selects the baseline group (releveled to the front);
each non-reference equation's exposure coefficient is that subtype's log odds
ratio versus the reference, so `contrast(type = "or")` reports one OR per
(subtype, exposure-coefficient) with an information-matrix Wald interval, and
`tidy()` renders the full per-equation table with a `y.level` column. The new
`R/polytomous.R` (`fit_polytomous()`, `contrast_polytomous()`,
`multinom_exposure_or()`, `tidy_multinom()`) and `resolve_polytomous_outcome()`
carry the engine; the dispatch grew an `outcome_kind` axis so `matcha()` picks
the multi-group resolver instead of the binary one.

- A two-group, numeric, or logical outcome is `matchatr_bad_outcome` (routed back
  to the binary `logistic` / `mh` / `clogit` estimators); an out-of-range
  `reference`, a `reference` on a non-polytomous estimator, an ordered-factor
  exposure, and `effect_modifier` are each `matchatr_bad_input`; a constant
  exposure — or one collinear with a confounder — is
  `matchatr_unestimable_exposure`. Because `nnet::multinom` (unlike `stats::glm`)
  does not alias a rank-deficient column to `NA` but splits the coefficient
  across the dependent columns (a silently attenuated OR),
  `reject_collinear_exposure()` catches it by the design-matrix rank, mirroring
  glm's NA-aliasing; collinearity confined to the confounders leaves the
  exposure estimable and is not rejected. `polytomous` on a matched design is
  `matchatr_bad_estimator`. RD / RR stay `matchatr_unidentified_estimand` and
  sandwich / bootstrap variances are `matchatr_unsupported_variance` (the engine
  reports the multinomial information matrix only).
- Validated against two oracles: the saturated 3-group / binary-exposure
  multinomial reproduces the closed-form 2×2 Woolf log-OR **and** its variance
  per subtype (independent of multinom's own `vcov()`), and a truth-based cohort
  DGP recovers the known per-subtype exposure log-ORs within 3.5 SE, with
  `nnet::multinom` coefficient / variance equality pinning the adjusted,
  continuous, and factor-exposure fits. `nnet` is a recommended R package, so it
  joins `Imports` without adding an install-time dependency.

## 2026-06-03 — Effect modification in matched case-control (PHASE_3 Chunk 3)

Completes the matched case-control layer with stratum-specific odds ratios.
Passing `effect_modifier = "m"` to `matcha(estimator = "clogit")` fits
`outcome ~ exposure * m + confounders + strata(set)`, and `contrast(type = "or")`
reports the exposure's conditional odds ratio within each modifier level — the
exposure main coefficient `beta_x` at the modifier's reference level and the
linear combination `beta_x + beta_{x:level}` elsewhere — with Wald intervals
from the joint partial-likelihood variance (so the per-level SE accounts for the
main–interaction covariance). The new `R/effect_modification.R`
(`stratum_specific_or_result()`) assembles the per-level ORs via a contrast
matrix `C V C'`; the modifier may coincide with a matching variable, which is
the canonical use (does the exposure effect differ across the matching factor?).

- Supported for a single-coefficient exposure (binary, continuous, or two-level
  factor) crossed with a categorical (logical / character / factor) modifier; a
  character / logical modifier is coerced to a factor. A 3+-level factor exposure
  is `matchatr_unsupported_combination`, a continuous (numeric) modifier or a
  non-`clogit` engine is `matchatr_bad_input`, a modifier coinciding with the
  outcome / exposure is `matchatr_bad_input`, and an aliased per-level
  interaction is `matchatr_unestimable_exposure`. RD / RR stay unidentified and
  sandwich / bootstrap intervals remain declined.
- Validated three ways: the within-level 1:1 McNemar closed form pins both the
  per-level point estimate (`OR = n10/n01`) and the linear-combination variance
  (`Var(log OR) = 1/n10 + 1/n01`) independently of `survival::clogit`;
  hand-built `survival::clogit` linear combinations check forwarding, the level
  labelling, and the interaction-column ordering for a three-level modifier with
  covariate adjustment; and a truth DGP with known per-level conditional log-ORs
  confirms recovery within a self-scaling SE band.
- M:1-with-covariates and variable-ratio matching needed no new code (the
  conditional likelihood handles any matched-set composition); their truth-based
  cells were already covered by Chunk 1's `clogit` engine.
- The effect modifier is `droplevels()`'d before fitting: a modifier factor
  carrying an unused (zero-observation) level previously added an all-zero
  interaction column aliased to `NA` and aborted the whole stratum-specific OR
  as unestimable. Unused levels are now dropped (preserving the order of the
  remaining declared levels, so a user-set reference is kept).

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
