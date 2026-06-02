# Phase 9 — Case-Control-Weighted Marginal Causal Contrasts

> **Status: DESIGN**
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
| unmatched CC | ccw_gformula | RD/RR/mOR | sandwich (causatr) / boot | needs-test |
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

1. `cc_weights()` + CCW-g-formula via causatr + pseudo-cohort oracle + missing-q₀
   rejection.
2. CCW-IPW + CCW-AIPW via causatr + double-robustness tests.
3. CCW-TMLE targeting step (new) + EIF variance + `tmle` oracle.
4. Estimated-q₀ variance correction + matched/nested CC support + bootstrap.

## Deferred items

Time-to-event marginal estimands (Phase 10), calibrated weights (Phase 12), secondary /
continuous outcomes (Phase 18), transportability (would compose with causatr transport).
