# matchatr (development version)

## 2026-06-12 — CCW-TMLE: case-control-weighted targeted learning (PHASE_9 Chunk 3)

`matcha(estimator = "ccw_tmle")` adds targeted maximum likelihood estimation to
the case-control-weighting family — the one genuinely new engine, since the
etverse has no targeted-learning code to delegate to. It reports the same marginal
risk difference / risk ratio / marginal odds ratio as the other CCW estimators,
and like CCW-AIPW it is **doubly robust** (consistent if either the outcome or the
propensity model is correct).

The targeting step (`R/tmle_ccw.R`, van der Laan & Rubin 2006; van der Laan & Rose
2011) runs entirely on the case-control-weighted sample: an initial weighted
logistic outcome model Q̄⁰(A, W); a weighted propensity g(W) bounded away from
0/1; the clever covariate H(A, W) = A/g(W) − (1 − A)/(1 − g(W)); a weighted
logistic fluctuation of Y on H with offset logit Q̄⁰ giving the tilt ε; the update
Q̄*(a, W) = expit(logit Q̄⁰(a, W) + ε H(a, W)); and the marginalized
treatment-specific means. The variance is the efficient influence function
weighted by the case-control weights (delta-method log-scale intervals for the RR
and OR). The shared `ccw_prepare()` (factored out of `fit_ccw()`) builds the
weights and 0/1-coded sample for every CCW engine, and `tidy()` / `summary()` now
branch on the engine being a CCW estimator (CCW-TMLE's model is a
`matchatr_ccw_tmle`, not a `causatr_fit`).

Validated in `test-tmle_ccw.R` against `tmle::tmle(obsWeights = )` on the same
case-control weights (plain-glm initial fit, matching gbound): the targeted risk
difference and its SE agree essentially exactly, the risk / odds ratios within
~1%. A truth DGP recovers the analytical marginal RD / RR / mOR, and a
double-robustness DGP confirms CCW-TMLE recovers the marginal truth whether the
outcome or the propensity working model is misspecified. `tmle` is added to
`Suggests`.

**Unified missing-data handling for the CCW family.** `ccw_prepare()` now
complete-cases the analysis sample for every CCW engine — rows with a missing
outcome, exposure, or confounder are dropped with a classed `matchatr_dropped_rows`
warning, and the case-control weights are computed on the complete-case sample so
the weighted case fraction still equals q0. This is the matchatr-wide convention
(it matches the classical engines and causatr's default `na.action`), and it
unifies what was inconsistent behaviour: CCW-g-formula silently complete-cased,
CCW-IPW / CCW-AIPW errored on a confounder NA, and CCW-TMLE crashed (its
hand-rolled prediction / clever-covariate vectors misaligned against the
full-length data). Multiple imputation (the recommended approach for missing
confounders; Dashti et al. 2024) and an outcome-missingness / IPCW extended-TMLE
(for missing outcomes) are the deferred principled alternatives (PHASE_13). A
non-converging TMLE fluctuation now warns (`matchatr_tmle_convergence`) and falls
back to the untargeted initial fit rather than propagating `NaN`.

## 2026-06-11 — CCW-IPW and doubly-robust CCW-AIPW (PHASE_9 Chunk 2)

Two more case-control-weighted marginal estimators: `matcha(estimator =
"ccw_ipw")` (inverse-probability weighting) and `matcha(estimator = "ccw_aipw")`
(augmented IPW). Both report the same marginal risk difference / risk ratio /
marginal odds ratio as `ccw_gformula`, on a case-control sample with a known
prevalence q0. **CCW-AIPW is doubly robust** — consistent for the marginal effect
if **either** the outcome model or the propensity model is correctly specified
(Rose & van der Laan 2014, *Biometrics* 70(1)).

The implementation is delegation-first: `fit_ccw()` (`R/ccw.R`) is parameterized
over the estimator, mapping `ccw_gformula` / `ccw_ipw` / `ccw_aipw` to
`causatr::causat(estimator = "gcomp" / "ipw" / "aipw", weights = cc_weights)` —
the case-control weights enter as causatr's observation weights — and reuses the
same `contrast_ccw()` standardization unchanged. The IPW / AIPW propensity fitter
is named explicitly (`propensity_model_fn = stats::glm`) so causatr does not warn
about defaulting it. Point estimate and influence-function/sandwich variance are
causatr's; matchatr owns only the weighting. The missing-q0 / non-binary-exposure
/ missing-confounders / off-scale / bootstrap rejections are shared across the CCW
family.

Validated in `test-ccw.R`: the exact pseudo-cohort `causatr` oracle and the
marginal-truth recovery now cover all three estimators, plus a **double-robustness**
test — a cohort where exactly one of the `~ w`-linear working models is correct
(the other omits a quadratic term); CCW-AIPW recovers the analytical marginal
risk difference whether the outcome or the propensity model is the misspecified
one, while the corresponding singly-robust estimator (g-formula or IPW) is biased.

## 2026-06-11 — Case-control-weighted marginal causal contrasts (PHASE_9 Chunk 1)

First **marginal** causal effect from a case-control sample. `matcha(estimator =
"ccw_gformula")` on an unmatched case-control design carrying a known source
prevalence q0 (`unmatched_cc(prevalence = q0)`) reports the marginal risk
difference (`contrast(type = "difference")`, the default), risk ratio
(`type = "ratio"`), or marginal odds ratio (`type = "or"`) — the quantities a
plain conditional odds ratio cannot deliver (non-collapsible, no baseline risk).

The estimator is the Rose & van der Laan case-control-weighted g-formula: the
new `cc_weights()` (`R/weights_cc.R`) computes the weights q0 / (n1/n) for cases
and (1 − q0) / (n0/n) for controls, which reweight the sample's outcome margin to
the source population so the weighted empirical distribution mimics the cohort;
`fit_ccw()` (`R/ccw.R`) then fits a weighted g-computation via
`causatr::causat(estimator = "gcomp")` and `contrast()` standardizes it to the
marginal effect over the treat-all / treat-none static interventions by
forwarding to `causatr::contrast()`. Point estimate and influence-function /
sandwich variance are delegated to causatr; matchatr owns only the weighting
layer. A non-binary exposure (`matchatr_bad_input`), a missing q0
(`matchatr_missing_prevalence`), an off-scale contrast
(`matchatr_unidentified_estimand`), and a bootstrap interval
(`matchatr_unsupported_variance`) are each rejected. A marginal g-formula
contrast has no information-matrix variance distinct from the influence-function
one, so `ci_method = "model"` and `"sandwich"` both yield causatr's sandwich
interval and the result records `"sandwich"`; `tidy()` and `summary()` on a ccw
fit surface the marginal contrast (the fitted model is a `causatr_fit` with no
conditional coefficient table).

Validated in `test-ccw.R`: an **exact** pseudo-cohort oracle (the hand-weighted
`causatr::causat()` + `contrast()` agrees to machine precision, confirming the
weighting and intervention plumbing) and a **truth DGP** (a cohort with an
analytical marginal RD / RR / mOR; the case-control-weighted g-formula recovers
the marginal truth, which lands on a tolerance band disjoint from the conditional
odds ratio a logistic fit reports — the non-collapsibility pin). `cc_weights()`
is checked against its closed form (weighted case fraction == q0) in
`test-weights_cc.R`.

CCW-IPW / CCW-AIPW (Chunk 2), CCW-TMLE (Chunk 3, the one genuinely new targeting
engine), and the estimated-q0 variance correction with matched / nested CC
support and within-stratum bootstrap (Chunk 4) remain pending.

## 2026-06-10 — Time-varying additive excess risk for IPW NCC (PHASE_7 follow-up)

New exported verb `excess_risk(fit, times)` for an `ipw_aalen` fit reports the
**time-varying** Aalen cumulative regression functions
B_j(t) = ∫₀ᵗ β_j(s) ds — the cumulative excess hazard for each covariate, the
additive analogue of `absolute_risk()`. Where `contrast(type = "excess")` reports
one time-*constant* excess hazard per covariate (the Lin-Ying model),
`excess_risk()` relaxes the constant-effect assumption and returns the full
cumulative regression function over time, completing the Phase-7 deferred item
"time-varying additive effects / cumulative regression B(t)".

The weighted least-squares estimator (`aalen_cumulative()`, `R/excess_risk.R`)
accumulates dB̂(t_i) = (X̃ᵀWX̃)⁻¹ X̃ᵀW dN(t_i) over the event times of the
deduplicated, Samuelsen-weighted analysis sample, with the Aalen martingale
pointwise variance (X̃ᵀWX̃)⁻¹{Σ w² x̃x̃ᵀ}(X̃ᵀWX̃)⁻¹; the interval is a symmetric Wald
band on the linear scale (B can be negative). The estimator truncates with a
`matchatr_truncated_excess` warning if the weighted design becomes singular in a
sparse late risk set. `print` / `tidy` methods render the per-(covariate, time)
table; non-`ipw_aalen` engines are rejected with `matchatr_not_implemented`.

Validated in `test-excess_risk.R`: B̂_j(t) and the pointwise SE match
`timereg::aalen` (without `const()`) `cum` and `var.cum` to machine precision
(1e-8), including a three-level factor exposure; a known constant-excess-hazard
truth DGP recovers B_x(t) = β_x·t within a SE band.

## 2026-06-10 — Split R/weighted_cox.R (internal refactor)

The Samuelsen IPW weighted Cox engine and its shared helpers
(`fit_ipw_cox()` / `contrast_ipw_cox()` / `ncc_ipw_analysis_data()` /
`require_ipw_ncc_columns()`) move from `R/weighted_cox.R` into a new
`R/ipw_cox.R`; `R/weighted_cox.R` now holds only the counter-matched
`fit_weighted_cox()` / `contrast_weighted_cox()`. Pure code move, no behaviour
change — both files are now under the ~300-line guideline. No user-visible
change.

## 2026-06-10 — Non-Weibull AFT distributions for IPW NCC (PHASE_7 follow-up)

`matcha(estimator = "ipw_aft")` gains a `dist` argument selecting the
accelerated-failure-time baseline — `"weibull"` (default), `"exponential"`,
`"lognormal"`, or `"loglogistic"` — completing the Phase-7 deferred item
"non-Weibull AFT distributions". All four are log-location-scale AFT models, so
`contrast(type = "af")` reports the same time-ratio estimand exp(β) under each;
they differ in the baseline error distribution (and therefore the survival-curve
shape). `"exponential"` is the one-parameter Weibull (it fixes σ = 1).

`absolute_risk()` follows the distribution: the cumulative incidence is
F̂_x(t) = G((log t − η̂)/σ̂), where G is the baseline error CDF — extreme-value
(complementary log-log) for weibull/exponential, Φ for lognormal, plogis for
loglogistic. The delta-method CI is the Wald interval on the standardised
residual mapped through the monotone G; for the extreme-value baselines this is
exactly the cloglog inversion the Cox-type engines share. The exponential's fixed
scale carries no log-scale parameter, so the scale term correctly drops out of
the gradient.

`dist` joins `model_fn` / `effect_modifier` / `reference` as an estimator-specific
`matcha()` argument: supplying it for a non-AFT estimator, or naming an
unsupported `survreg` distribution, is `matchatr_bad_input`.

Validated in `test-aft_ncc.R` / `test-absolute_risk_aft.R`: each distribution's
coefficient/SE matches an independent `multipleNCC::KMprob` + `survreg`
reconstruction (1e-6); each survival curve round-trips through
`predict.survreg(type = "quantile")` and matches a numDeriv reconstruction
through its error CDF (1e-7).

## 2026-06-10 — AFT survival-curve absolute risk for IPW NCC (PHASE_7 follow-up)

`absolute_risk()` gains an `ipw_aft` engine path, completing the Phase-7 deferred
item "AFT acceleration-factor absolute risk". The fitted weighted Weibull
accelerated failure time model is a parametric survival curve, so the cumulative
incidence is read directly off the coefficients,

  F̂_x(t) = 1 − exp(−exp((log t − η̂) / σ̂)),

where η̂ is the AFT linear predictor and σ̂ the scale — no Breslow step function.
Pointwise intervals use the delta method on the complementary log-log scale over
θ = (β, log σ), with the gradient ∂ξ/∂β = −x̃/σ, ∂ξ/∂(log σ) = −ξ and the robust
Lin-Wei sandwich `survival::survreg(robust = TRUE)` stores in `vcov()`; the
interval is inverted to the risk scale by the shared cloglog inversion the
Cox-type engines already use (`cloglog_risk_ci()` / `new_matchatr_absolute_risk()`,
factored out of `assemble_absolute_risk()` in `R/absolute_risk.R`; the engine is
`R/absolute_risk_aft.R`). The result is weight-agnostic — KM (Samuelsen) and
GLM/GAM working-model weights both feed it through the same fit.

Validated in `test-absolute_risk_aft.R`: F̂_x(t) round-trips through
`predict.survreg(type = "quantile")` (survival's own inverse CDF) to 1e-7; the
estimate and CI match an independent `numDeriv` reconstruction of the ξ(θ)
gradient (including factor contrasts); the NCC subsample recovers the full-cohort
`survreg` curve within sampling tolerance; and a Weibull truth DGP's analytical
F_x(t) is covered by the CI. `absolute_risk()` on the additive (`ipw_aalen`)
engine — which has no survival-curve verb — is rejected with
`matchatr_not_implemented`.

## 2026-06-10 — Additive-hazards and AFT models for IPW NCC (PHASE_7 Chunk 5)

Two non-Cox alternative models on the deduplicated Samuelsen-weighted NCC sample
(Handbook Ch19 §19.5), completing Phase 7.

- **Accelerated failure time** (`estimator = "ipw_aft"`, `R/aft_ncc.R`): a
  weighted Weibull AFT via `survival::survreg(weights = ipw_weight, robust = TRUE)`;
  `contrast(type = "af")` reports the time ratio exp(β) — the acceleration factor,
  the factor by which a unit of exposure multiplies survival time (Kang, Lu & Liu
  2017, Biometrics 73(1)). The robust sandwich is `survreg`'s; the time ratio
  shares the exponentiated-coefficient interval shape of the odds / hazard ratios.
- **Additive hazards** (`estimator = "ipw_aalen"`, `R/additive_ncc.R` +
  `R/lin_ying.R`): the weighted constant additive-hazards model (Lin & Ying 1994,
  Biometrika 81(1)), implemented in matchatr rather than delegated — the point
  estimate γ̂ = A⁻¹B is a closed form and the robust variance is the
  martingale-residual sandwich A⁻¹(Σ η̂_iη̂_iᵀ)A⁻¹. `contrast(type = "excess")`
  reports the excess hazard γ (additive rate difference; Borgan & Langholz 1997,
  Biometrics 53(2)), a linear-scale, possibly-negative estimand whose Wald
  interval is symmetric (not exponentiated). A constant or collinear covariate
  makes the weighted design matrix singular and aborts with
  `matchatr_unestimable_exposure` rather than a raw LAPACK error. `timereg` is a
  test-only oracle, not a runtime dependency.

Two new contrast scales (`type = "af"`, `type = "excess"`); each engine identifies
exactly one and rejects the others, bootstrap/sandwich variance, non-`incl_prob`
data, and non-nested designs with classed errors.

Validated in `test-aft_ncc.R` / `test-additive_ncc.R`: AFT recovers the
full-cohort `survreg` coefficient (3.5-SE) and matches an independent `KMprob` +
`survreg` reconstruction to machine precision; the additive estimator recovers a
known excess hazard (3.5-SE, binary and three-level-factor exposure) and matches
`timereg::aalen` exactly on the point estimate (full coefficient vector, including
a complex continuous-exposure / factor-confounder set with heavy ties) and within
5% on the robust SE.

## 2026-06-10 — Multiple endpoints from one reused NCC control set (PHASE_7 Chunk 4)

Reuses a single nested case-control control set to estimate hazard ratios for
more than one endpoint, the IPW reformulation's payoff over the matched analysis
(Samuelsen 1997; Støer & Samuelsen 2012). Two modes feed the existing `ipw_cox`
weighted Cox:

- **Combined-event reuse.** Drawing the NCC on the union "any-failure" event
  ascertains every endpoint's cases at once, so each cause-specific endpoint is
  analysed directly: `matcha(ncc, outcome = "<cause>", estimator = "ipw_cox")`.
  This required generalising the analysis-sample builder (`ncc_ipw_analysis_data()`,
  `R/weighted_cox.R`) to keep weight 1 for any subject that is a case of the
  analysed endpoint **or** the failing subject of some sampled risk set — a
  competing-endpoint case is ascertained by the sampling and must not revert to
  the control weight 1/π_j on a row where it happened to be drawn as a control.
  For the single-endpoint analysis the two clauses coincide, so this is a no-op
  generalisation (the `multipleNCC::wpl` exact-agreement tests still hold).
- **Cohort-augmented reuse.** New exported `reuse_ncc_endpoint(ncc, cohort, time,
  event)` (`R/multi_endpoint.R`) reuses a control set sampled for a *primary*
  endpoint to fit a *secondary* one: the controls keep their primary Samuelsen
  inclusion weights 1/π_j (the inclusion probability is a property of the
  sampling, not the endpoint), and the secondary endpoint's cases that were not
  sampled are augmented from the Phase-1 cohort with weight 1. Requires the
  cohort; `matchatr_missing_phase1` / `matchatr_bad_input` / `matchatr_bad_outcome`
  guard the inputs.

Tested in `test-multi_endpoint.R` against: `multipleNCC::wpl` (exact agreement on
log-HR and SE for each endpoint of a combined-event NCC, since `wpl` ascertains
all cases and computes π over all event times); an independent `KMprob` +
`survival::coxph` reconstruction of the cohort-augmented fit (machine precision);
and a competing-risks truth DGP with known cause-specific Cox log-HRs (each mode
recovers the full-cohort HR within a 3.5-SE band). Additive/AFT models
(Ch19 §19.5) remain deferred.

`reuse_ncc_endpoint()` rejects two malformed inputs with `matchatr_bad_input`
rather than producing degenerate output: an empty NCC sample (whose augmented
set ids would come from `max(integer(0)) = -Inf`), and an augmented secondary
case with a missing / non-finite event time (which the downstream weighted Cox
would silently drop) — mirroring `sample_ncc()`'s rejection of cases whose risk
set cannot be formed.

## 2026-06-09 — IPW absolute risk correctness for complicated designs (critical review)

Three correctness fixes to the IPW nested case-control absolute risk, all in the
"complicated case" regime the initial tests (plain factors, continuous times) did
not exercise.

- **Data-dependent confounder bases** (`poly()`, `ns()`/`bs()`, `scale()`):
  `ar_lp_from_newdata()` rebuilt the design from `term.labels` and called
  `model.matrix()` on `newdata` alone, recomputing the basis from those rows — a
  *different* basis than the fit, with the *same* coefficient names (so the
  name-based guard passed) and a silently wrong linear predictor (F off by ~0.4).
  It now builds the design from the fitted model's `terms` (reusing the
  `predvars` basis) and `xlevels` for the `ipw_cox`/coxph engine. The `cch` engine
  keeps the original-formula path (its non-standard internal formula has no usable
  `predvars`).
- **Tied event times**: `fit_ipw_cox()` fitted the weighted Cox with the coxph
  default Efron ties, but `ipw_breslow_ncc()` uses the plain Breslow baseline, so
  the two disagreed at tied event times (F off by ~0.03). The fit now uses
  `ties = "breslow"` so the partial-likelihood coefficients and the Breslow
  cumulative baseline hazard are mutually consistent; for continuous failure times
  (the usual NCC setting) this is identical to Efron.
- **Event at the time origin**: the `t = 0` fence post prepended on top of a real
  `t = 0` event time (possible with rounded times) duplicated the knot and made
  `approx()` collapse it. `breslow_step_with_fence()` now skips the fence when an
  event already sits at `t = 0`.
- The hand-rolled Breslow now agrees with `survival::survfit` to machine precision
  across all three regimes; regression tests in `test-absolute_risk_ncc.R` cover a
  tied-time design, a `poly(z, 2)` confounder, and an event at `t = 0`.

## 2026-06-09 — IPW Breslow absolute risk for nested case-control (PHASE_7 Chunk 3)

Extends `absolute_risk(fit, newdata, times)` to the IPW nested case-control
(`"ipw_cox"`) engine, alongside the existing case-cohort (`"cch"`) path. Returns
`F̂_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t))` for each covariate pattern at each time.

- **Native IPW Breslow cumulative baseline hazard** (`ipw_breslow_ncc()`,
  `R/absolute_risk_ncc.R`) over the deduplicated, Samuelsen-weighted NCC analysis
  sample (cases at weight 1, each unique control at 1/π_j): the increment
  `dΛ̂₀(t_k) = (Σ events) / (Σ_{at risk} w_j exp(β̂ᵀ x_j))` is the Horvitz-Thompson
  weighted Breslow step, so the upweighted controls stand in for the unsampled
  cohort. The hand-rolled step function agrees with `survival::survfit` on the
  same weighted Cox to machine precision — across KM and GLM/GAM weights, factor
  and continuous confounders, and multiple covariate patterns.
- **Pointwise CIs** via the delta method on the complementary log-log scale,
  `Var(log Λ_x(t)) = x' V_β x + Var(log Λ̂₀(t))` with the IPW robust sandwich for
  the β part. Monte Carlo coverage is conservative (errs wide, never
  anti-conservative); the baseline-hazard variance is the within-sample
  Nelson-Aalen approximation, as in the case-cohort path.
- **Refactor**: the per-time / per-pattern F_x(t) assembly and the delta-method CI
  are shared by both engines via `assemble_absolute_risk()` (`R/absolute_risk.R`);
  the case-cohort Breslow moved to `R/absolute_risk_cch.R`. The deduplication of
  the NCC analysis sample is shared between the weighted Cox fit and the Breslow
  via `ncc_ipw_analysis_data()` (`R/weighted_cox.R`).
- A non-`cch`/`ipw_cox` engine still aborts with `matchatr_not_implemented`
  (message updated to name both supported engines).
- Validated in `tests/testthat/test-absolute_risk_ncc.R`: exact `survival::survfit`
  agreement (1e-8) including a GLM-weighted factor-confounder design, full-cohort
  survfit cross-check, and analytical exponential truth recovery with CI coverage.

## 2026-06-09 — Working-model inclusion-probability weights for NCC (PHASE_7 Chunk 2)

Adds `compute_ncc_weights()` for GLM and GAM working-model inclusion probabilities.

- **`compute_ncc_weights(ncc, cohort, method, selection_formula, time, entry)`**
  (new, `R/weights_design.R`) takes an NCC dataset from
  `sample_ncc(incl_prob = TRUE)` and the full Phase-1 cohort, builds the
  augmented (eligible subject × event time) selection dataset, fits a logistic
  regression (`method = "glm"`, [stats::glm()]) or generalised additive model
  (`method = "gam"`, [mgcv::gam()]) for the binary control-selection indicator,
  and replaces the `ipw_weight` column with the resulting working-model inverse
  inclusion probabilities. Cases are forced to weight 1.
- **`matchatr_missing_phase1`** fires when `cohort = NULL` or the `time` column
  is absent from `cohort`: both signals that the Phase-1 event times needed to
  reconstruct risk sets are unavailable.
- **Default `selection_formula = ~ risk_time`** — a time-only logistic model
  matching the simplest GLM specification of Borgan, Samuelsen & Aastveit
  (2003). Users can extend with cohort covariates.
- **Oracle**: `multipleNCC::wpl(weight.method = "glm")` — log-HR agreement
  within 2e-2 (minor formula difference in default time term between the two
  implementations).

## 2026-06-09 — IPW NCC vignette DGP and test quality fixes (critical review)

Two correctness issues in the IPW NCC vignette and test suite.

- **Vignette DGP**: `vignettes/ipw-ncc.qmd` generated the observed time `t` and
  event indicator `d` from two independent `rexp(n, rate)` draws, producing
  inconsistent censoring (658 subjects with `t < 5` but `d = 0`, 675 with
  `t = 5` but `d = 1`). The biased DGP recovered HR ≈ 2.39 instead of the
  true 2.0. Fixed to a single draw: `tt <- rexp(n, rate); t = pmin(tt, 5);
  d = as.integer(tt <= 5)`. Also removed a dead first cohort block that was
  immediately overwritten.
- **Rejection test outcome column**: the `ipw_weight`-missing rejection test
  passed `"case"` (the per-set NCC indicator) as the outcome; changed to `"d"`
  (the cohort event indicator) to match documented usage. The test was
  functionally correct (it errors before fitting) but modelled the wrong API
  call.

## 2026-06-09 — Samuelsen KM IPW weights + `ipw_cox` engine for NCC data (PHASE_7 Chunk 1)

Implements the IPW reformulation of nested case-control data: break the
matching, weight each unique control by the inverse of its Samuelsen
Kaplan-Meier inclusion probability, and fit a standard weighted Cox model.

- **`sample_ncc(incl_prob = TRUE)`** gains two new output columns. `.cohort_row`
  records the original cohort row index for each NCC row (used to deduplicate
  controls that were sampled into multiple risk sets). `ipw_weight` holds the
  Samuelsen (1997) KM inverse inclusion probability 1/π_j, where
  π_j = 1 − prod(1 − m_i/n_elig_i) over all event times where j was eligible
  (not just those where j was actually sampled); cohort cases are forced to
  weight 1. The formula is computed in O(n × K) via the internal
  `samuelsen_km_weights()` helper.
- **`matcha(design = nested_cc(...), estimator = "ipw_cox")`** deduplicates the
  NCC data by `.cohort_row`, fits `survival::coxph(weights = ipw_weight, robust
  = TRUE)`, and `contrast()` reports the exposure's hazard ratio with the
  Lin-Wei robust sandwich variance. The default `type = "hr"` is inherited from
  the nested design; `type = "or"` and `ci_method = "bootstrap"` are rejected.
- Missing `ipw_weight` or `.cohort_row` columns abort with
  `matchatr_missing_ipw_weights` / `matchatr_bad_input`, pointing to
  `sample_ncc(incl_prob = TRUE)`.
- Validated: KM weights match the closed-form KM probability formula to 1e-8;
  truth-based DGP recovers the full-cohort Cox HR within 3.5 SE; exact
  agreement with `multipleNCC::wpl(weight.method = "KM")` on log-HR and SE
  (tolerance 1e-6). New vignette `vignettes/ipw-ncc.qmd` demonstrates the full
  pipeline alongside the classical NCC and full-cohort Cox comparisons.

## 2026-06-09 — IPW Breslow absolute risk from case-cohort fits (PHASE_6 Chunk 3)

New exported verb `absolute_risk(fit, newdata, times)` for `matchatr_fit` objects
fitted with the `"cch"` engine. Returns `F̂_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t))`
at each evaluation time for each row of `newdata`.

- **IPW Breslow cumulative baseline hazard** `Λ̂₀(t)`: at each event time t_k the
  denominator is scaled by `N / n_sub` (cohort size over subcohort size) so the
  subcohort-only risk set correctly represents the full cohort. For Borgan I/II
  stratified methods the per-stratum weight `N_s / n_sub_s` is used instead.
- **Pointwise CIs** via the delta method on the complementary log-log scale
  (`log(−log(1 − F))`) — the standard log-log CI for the survival function inverted
  to the cumulative-incidence scale. Variance = `x′ V_β x + Σ (dΛ̂₀)² / Λ̂₀²`.
- **All five `cch` methods** (Prentice, SelfPrentice, LinYing, I.Borgan, II.Borgan)
  are handled; Borgan I/II use per-stratum IPW weights.
- Times before the first event return `F̂ = 0` exactly (fence-post at `t = 0` in
  the Breslow step function). CIs are clamped to `[0, 1]`.
- `tidy.matchatr_absolute_risk()` returns the long-form estimates `data.table`.
- Non-`cch` engines abort with `matchatr_not_implemented`.

## 2026-06-09 — Remove dead n_dropped warning in fit_cch (critical review)

`fit_cch()` computed `n_dropped = nrow(subset_dt) - model$n`, but `model$n` for
a `survival::cch` object equals `subcohort_size + n_events` — it double-counts
subjects who are both subcohort members and cases. `n_dropped` was therefore
always negative and the warning could never fire. Removed the dead code; NA
handling is left to `survival::cch`'s own `na.action`.

## 2026-06-09 — Borgan I/II IPW case-cohort estimators (PHASE_6 Chunk 2)

Extends the `cch` engine with the Borgan I/II IPW estimators for stratified
subcohort sampling.

- **`case_cohort(..., method = "I.Borgan", stratum = "<col>")`** and `"II.Borgan"`
  now work end-to-end. `stratum` names the column that defines the subcohort
  sampling strata; `survival::cch` weights each subject by the inverse of its
  stratum-specific subcohort sampling fraction (Borgan et al. 2000).
- **`cohort.size` is computed per stratum** from the full cohort (not the
  cases+subcohort subset), so the denominator correctly reflects the total stratum
  size. A named integer vector is passed to `survival::cch`.
- `method = "I.Borgan"` or `"II.Borgan"` without `stratum` aborts with
  `matchatr_bad_design`; a missing stratum column in the data aborts similarly.
- `case_cohort()` gains the `stratum` parameter; `print.matchatr_design()` shows
  the stratum column when present.
- Tests in `test-case_cohort.R`: nwtco oracle for both Borgan methods (exact
  coefficient and SE match vs direct `survival::cch`), truth DGP with two sampling
  strata (`make_stratified_case_cohort_data()`), and rejection path snapshot.

## 2026-06-08 — Case-cohort Cox pseudo-likelihood (PHASE_6 Chunk 1)

Adds the `case_cohort()` design constructor and the `cch` engine wrapping
`survival::cch()` for the Prentice, Self-Prentice, and Lin-Ying pseudo-likelihood
hazard ratio.

- **`case_cohort(subcohort, time, method, id)`** declares a case-cohort sampling
  structure. `method` selects the pseudo-likelihood variant (`"Prentice"`,
  `"SelfPrentice"`, `"LinYing"`, `"I.Borgan"`, `"II.Borgan"`); `id` names the
  subject-identifier column so `survival::cch` correctly pairs subjects that
  appear as both subcohort member and case. Default: `method = "Prentice"`.
- **`matcha(design = case_cohort(subcohort, time, method), estimator = "cch")`**
  builds the `Surv(time, status) ~ exposure + confounders` call, subsets to
  cases + subcohort members (censored non-members contribute nothing to the
  pseudo-likelihood), and delegates to `survival::cch` with `cohort.size =
  nrow(data)`. `contrast()` reports the exposure's **hazard ratio** (`type = "hr"`)
  using the variance that `survival::cch` returns for the chosen method — the
  correct asymptotic / pseudo-likelihood variance, not the naive information matrix.
- The pseudo-likelihood's dependent score factors mean its controls (subcohort
  members) are reused across failure times, so the naive information-matrix SE is
  never used. Self-Prentice and Prentice share the same asymptotic variance; LinYing
  uses an independent robust variance estimator.
- `type = "or"`, `"difference"`, `"ratio"` and `ci_method = "sandwich"`,
  `"bootstrap"` are rejected for this estimator.
- New functions: `fit_cch()`, `contrast_cch()`, `cch_exposure_coef_names()` in
  `R/case_cohort.R`. Wired in `R/dispatch.R`, `R/contrast.R`, `R/cc_design.R`.
- Tests in `tests/testthat/test-case_cohort.R`: nwtco oracle (3 methods vs direct
  `survival::cch` call), truth recovery within 3.5 SE, full-cohort `coxph`
  agreement, structural checks, and all rejection paths.

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
