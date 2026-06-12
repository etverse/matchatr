# Phase 9 — Case-Control-Weighted Marginal Causal Contrasts

> **Status: COMPLETE — Chunks 1–4 done.**
> `matcha(estimator = "ccw_gformula" | "ccw_ipw" | "ccw_aipw" | "ccw_tmle")` reports
> the marginal RD / RR / marginal OR from an unmatched **or matched** case-control
> sample with a known (or cohort-estimated) q0, via `cc_weights()` + `causatr`
> g-computation / IPW / AIPW (`R/ccw.R`) or matchatr's own targeting engine
> (`R/tmle_ccw.R`); CCW-AIPW and CCW-TMLE are doubly robust. Variance is the causatr
> sandwich / TMLE EIF, optionally widened for an estimated q0 (Chunk 4b) or replaced
> by the design-preserving within-stratum bootstrap (Chunk 4a). A nested (risk-set)
> CC design is rejected toward `ipw_cox`. Methods: Rose & van der Laan (2008, 2009,
> 2011 *Targeted Learning*, 2014 double-robust case-control). Implements Track 2 of
> `PHASE_8_CAUSAL_STRATEGY`.

## Scope

**In:** marginal causal effects (RD, RR, marginal OR) from binary-outcome case-control
samples via case-control weighting — CCW-g-formula, CCW-IPW, CCW-AIPW (all reusing
`causatr`), and CCW-TMLE (new targeting step). Known and cohort-estimated prevalence q₀.
Independent, matched, and nested CC samples.

**Out:** time-to-event marginal estimands under sampling (Phase 10), calibration of
weights (Phase 12), continuous/secondary outcomes (Phase 18).

## Key design decisions

- **Case-control weights map the CC sample to the cohort.** For a sample with cases and
  controls, weight cases by q₀/(sample case fraction) and controls by (1−q₀)/(sample
  control fraction) (equivalently the Rose–van der Laan CCW form). The weighted
  empirical distribution mimics the source population, so a cohort estimator on the
  weighted sample targets the marginal estimand.
- **Reuse causatr for g-formula / IPW / AIPW.** The weights enter as causatr's
  observation `weights`. CCW-g-formula = `causat(estimator="gcomp", weights=w_cc)` +
  `contrast()`; CCW-IPW and CCW-AIPW likewise. matchatr does NOT reimplement these.
- **CCW-TMLE is the one new engine.** Initial outcome model (optionally weighted) +
  a fluctuation step that solves the case-control-weighted efficient influence equation
  (a weighted logistic tilt with the clever covariate H(A,W) = A/g(W) − (1−A)/(1−g(W))).
  This lives in matchatr because the etverse has no targeted-learning code.
- **Weights are observation weights, never a data column** (`fit$details$weights`).
- **q₀ fixed vs estimated** changes the variance (extra IF term when estimated). The fit
  records `prevalence_known = TRUE/FALSE`.
- **Matched CC + CCW** must respect the matching in the weights/standardization
  (standardize over the matching-variable distribution); document the constraint that the
  effect modifier / matching variable handling follows causatr's baseline-covariate rule.

## API design

```r
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(prevalence = 0.02),
              confounders = ~ age + smoke, estimator = "ccw_gformula")
contrast(fit, type = "difference", ci_method = "sandwich")   # marginal RD
contrast(fit, type = "ratio")                                # marginal RR

matcha(..., estimator = "ccw_aipw")   # doubly-robust (causatr AIPW)
matcha(..., estimator = "ccw_tmle")   # targeted (new fluctuation step)
```

## Support matrix

| Design | Estimator | Estimand | Variance | Status |
|---|---|---|---|---|
| unmatched CC | ccw_gformula | RD/RR/mOR | sandwich (causatr) | ✅ done (Chunk 1) |
| unmatched CC | ccw_ipw | RD/RR/mOR | sandwich (causatr) | ✅ done (Chunk 2) |
| unmatched CC | ccw_aipw | RD/RR/mOR | sandwich (DR, causatr) | ✅ done (Chunk 2) |
| unmatched CC | ccw_tmle | RD/RR/mOR | EIF (DR, new) | ✅ done (Chunk 3) |
| matched CC | ccw_gformula/ipw/aipw/tmle | RD/RR/mOR | sandwich/EIF + boot | ✅ done (Chunk 4c) |
| nested CC | ccw_* | — | — | ⛔ `matchatr_bad_estimator` → `ipw_cox` (Chunk 4c) |
| any ccw_* | — (no prevalence) | — | — | ⛔ `matchatr_missing_prevalence` |
| q₀ estimated | ccw_* | RD/RR/mOR | IF with extra term | ✅ done (Chunk 4b) |

## Implementation plan

- `R/weights_cc.R` — `cc_weights(prevalence, outcome, design)` → the q₀ weight vector;
  fixed-vs-estimated q₀ bookkeeping.
- `R/ccw.R` — `fit_ccw()` dispatch: build weights, call `causatr::causat()` with
  `estimator ∈ {gcomp, ipw, aipw}` + weights, wrap the result; record variance kind.
- `R/tmle_ccw.R` — the NEW targeting step: initial Q̄ fit, clever covariate, weighted
  logistic fluctuation, update, marginalize; efficient-influence-function variance.
- `R/variance_ccw.R` — IF assembly including the estimated-q₀ correction; bootstrap
  refit path (resample within case/control strata, recompute weights).

## Variance / inference notes

- CCW-g-formula / IPW / AIPW: causatr's sandwich on the weighted fit gives the
  point-channel variance; add the sampling/q₀ correction. When q₀ is estimated from the
  cohort, the IF gains a term for q̂₀'s variance.
- CCW-TMLE: variance from the efficient influence function evaluated at the targeted fit
  (standard TMLE plug-in EIF variance), weighted by the CC weights.
- Bootstrap: resample within case and control strata separately (preserve the design),
  recompute q₀ weights each replicate.

## Oracle testing strategy

- **causatr on the explicitly reweighted pseudo-cohort** is the primary oracle: build the
  q₀-weighted dataset by hand, run `causatr::causat()`/`contrast()`, and assert
  `fit_ccw()` matches.
- **Truth-based**: simulate a cohort with known marginal RD/RR; draw an unmatched CC
  sample; confirm CCW-g-formula/AIPW/TMLE recover the marginal truth (which the
  conditional OR does NOT equal — a contrast test pinning marginal vs conditional).
- **CCW-TMLE** vs R `tmle::tmle()` run with case-control weights / on the reweighted
  sample.
- **Double robustness**: misspecify either the outcome or the propensity model for
  CCW-AIPW/TMLE and confirm consistency persists.

## Chunk plan

The 2026-06-11 causatr-reuse audit confirmed Chunks 2 and 4 are **delegation-first**
(causatr already provides the engines and the weight-aware machinery); only Chunk 3 is
genuinely new code.

1. ✅ `cc_weights()` + CCW-g-formula via causatr + pseudo-cohort oracle + missing-q₀
   rejection. (`R/weights_cc.R`, `R/ccw.R`; `test-weights_cc.R`, `test-ccw.R`)
2. ✅ CCW-IPW + CCW-AIPW + double-robustness tests. **Delegation-first:** `fit_ccw()`
   is parameterized over `fit$estimator` → `causatr::causat(estimator = "ipw" |
   "aipw", weights = cc_weights, …)` (both accept external `weights`); `contrast_ccw()`
   is reused unchanged. AIPW gives the doubly-robust marginal estimator with no new
   variance engine. The double-robustness test uses a functional-form misspecification
   (a `~ w`-linear working model that omits a quadratic term) so exactly one of the
   outcome / propensity models is wrong through matchatr's single `confounders`
   argument; CCW-AIPW recovers the marginal truth either way. (`R/ccw.R`; `test-ccw.R`,
   `helper-dgp.R::make_dr_cohort_ccw()`.)
3. ✅ CCW-TMLE targeting step (**new code** — causatr has no targeted learning) + EIF
   variance + `tmle` oracle. The shared `ccw_prepare()` (factored out of `fit_ccw()`)
   builds the weighted sample; `fit_ccw_tmle()` runs the clever-covariate logistic
   fluctuation and the EIF variance; `contrast_ccw_tmle()` reports RD / RR / OR.
   (`R/tmle_ccw.R`; `test-tmle_ccw.R`, `helper-tmle-oracle.R`.)
4. Estimated-q₀ variance correction + matched/nested CC support + bootstrap, split
   into sub-chunks:
   - **4a ✅ within-stratum bootstrap.** `ci_method = "bootstrap"` for all four CCW
     engines: resample cases / controls separately (design-preserving, so the q₀
     weights stay fixed), refit, percentile interval; drops the bootstrap rejection.
     matchatr owns the stratified loop (`ccw_bootstrap_ci()`, `R/variance_ccw.R`) —
     causatr's plain bootstrap mixes the strata and cannot preserve n1 / n0.
     (`test-variance_ccw.R`.)
   - **4b ✅ estimated-q₀ variance.** `unmatched_cc(prevalence = q0, prevalence_n =
     N)` declares q₀ estimated from N cohort members; the analytic interval adds the
     delta-method term (∂ψ/∂q₀)²·q₀(1−q₀)/N (`ccw_estimated_q0_term()` /
     `ccw_apply_estimated_q0()`, `R/variance_ccw.R`) and the bootstrap redraws q₀*
     per replicate; the fit records `prevalence_known`. (`test-variance_ccw.R`.)
   - **4c ✅ matched / nested CC support.** `matched_cc()` gains `prevalence`
     (`prevalence_n`), so the CCW estimators run on a matched case-control sample;
     the matching variable is a baseline covariate (it must be in `confounders`, so
     the effect is standardized over its distribution rather than conditioned on the
     matched sets), with the documented Rose & van der Laan (2009) efficiency
     caveat. A **nested** CC is risk-set sampled, so binary q₀ reweighting does not
     identify a marginal effect: `matcha(design = nested_cc(...), estimator =
     "ccw_*")` is rejected (`matchatr_bad_estimator`) toward `ipw_cox`. (`R/cc_design.R`,
     `R/matcha.R`; `test-ccw.R`, `helper-dgp.R::make_matched_cohort_ccw()`.)

**Cross-phase note (PHASE_13):** `causatr::causat_mice()` is estimator-agnostic, so a CCW
fit (any of the above) is automatically poolable over a `mice` mids object for
missing-by-design covariates — PHASE_13 should delegate to it rather than build MI pooling.
`causatr::diagnose()` (positivity / balance / weight ESS on the reweighted pseudo-cohort)
is available on any CCW fit at zero cost.

## Deferred items

Time-to-event marginal estimands (Phase 10), calibrated weights (Phase 12), secondary /
continuous outcomes (Phase 18), transportability (would compose with causatr transport).

**Missing data.** The CCW family complete-cases (listwise deletion in `ccw_prepare()`,
with a `matchatr_dropped_rows` warning) as the shipped interim policy. The principled
alternatives — multiple imputation with interactions for missing confounders, and an
outcome-missingness / IPCW extended-TMLE for missing outcomes (Dashti et al. 2024) —
are owned by **PHASE_13** (see its §(3); the CCW marginal estimators reuse
`causatr::causat_mice` directly, no congeniality construction needed).
