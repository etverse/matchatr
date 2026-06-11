# Phase 9 — Case-Control-Weighted Marginal Causal Contrasts

> **Status: IN PROGRESS — Chunk 1 (CCW-g-formula) complete.**
> `matcha(estimator = "ccw_gformula")` reports the marginal RD / RR / marginal OR
> from an unmatched case-control sample with a known q0, via `cc_weights()` +
> `causatr` g-computation (`R/weights_cc.R`, `R/ccw.R`). Chunks 2–4 (CCW-IPW /
> AIPW, CCW-TMLE, estimated-q0 variance + matched/nested CC + bootstrap) pending.
> Methods: Rose & van der Laan (2008, 2009, 2011 *Targeted Learning*, 2014 double-robust
> case-control). Implements Track 2 of `PHASE_8_CAUSAL_STRATEGY`.

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
| unmatched CC | ccw_ipw | RD/RR | sandwich / boot | needs-test |
| unmatched CC | ccw_aipw | RD/RR | sandwich (DR) / boot | needs-test |
| unmatched CC | ccw_tmle | RD/RR | EIF (new) / boot | needs-test |
| matched CC | ccw_gformula/aipw | RD/RR | boot (+IF) | needs-test |
| nested CC | ccw_* | RD/RR | boot | needs-test |
| any ccw_* | — (no prevalence) | — | — | ⛔ `matchatr_missing_prevalence` |
| q₀ estimated | ccw_* | RD/RR | IF with extra term | needs-test |

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
2. CCW-IPW + CCW-AIPW + double-robustness tests. **Delegation-first:** parameterize
   `fit_ccw()` over `fit$estimator` → `causatr::causat(estimator = "ipw" | "aipw",
   weights = cc_weights, …)` (both accept external `weights`); reuse `contrast_ccw()`
   unchanged. AIPW gives the doubly-robust marginal estimator with no new variance engine;
   the double-robustness tests use causatr's per-component `confounders_outcome` /
   `confounders_treatment` to misspecify one model at a time.
3. CCW-TMLE targeting step (**new code** — causatr has no targeted learning) + EIF
   variance + `tmle` oracle.
4. Estimated-q₀ variance correction + matched/nested CC support + bootstrap.
   **Delegation-first for the bootstrap:** causatr's `refit_gcomp` / `refit_ipw` /
   `refit_aipw` already resample and re-apply external `weights`, so matchatr only adds the
   within-case/control-strata resample + per-replicate q₀ reweighting, then drops the
   current `ci_method = "bootstrap"` rejection in `contrast_ccw()`. Matched-CC marginal
   standardization maps to causatr's `contrast(by = )` path, not new machinery.

**Cross-phase note (PHASE_13):** `causatr::causat_mice()` is estimator-agnostic, so a CCW
fit (any of the above) is automatically poolable over a `mice` mids object for
missing-by-design covariates — PHASE_13 should delegate to it rather than build MI pooling.
`causatr::diagnose()` (positivity / balance / weight ESS on the reweighted pseudo-cohort)
is available on any CCW fit at zero cost.

## Deferred items

Time-to-event marginal estimands (Phase 10), calibrated weights (Phase 12), secondary /
continuous outcomes (Phase 18), transportability (would compose with causatr transport).
