# Phase 18 — Secondary Analysis of Case-Control Data

> **Status: DESIGN**
> Book chapters: 14 (Secondary Analysis of Case-Control Data).

## Scope

**In:** valid inference for a *secondary* outcome Y₂ (not the disease Y₁ that defined the
case-control sampling) — inverse-probability-weighted estimation (using only the
*relative* case/control sampling probabilities) and semiparametric maximum likelihood
(SPML2 / SPML3) that jointly models (Y₁, Y₂ | X); the AIPW hybrid (Sofer et al. 2017)
for robustness + efficiency. Binary and continuous Y₂.

**Out:** the primary-outcome analyses (Phases 2–9); transportability.

## Key design decisions

- **Weighted (IPW) route** is robust and simple: weight the Y₂ regression score by the
  inverse relative case/control sampling probability (absolute values not needed). No
  nuisance models, but inefficient. Reuse `survey::twophase`/`svyglm` or causatr with the
  relative-probability weights.
- **SPML route** is efficient but nuisance-model-sensitive: jointly model (Y₁, Y₂|X) via
  SPML2 (f₂(Y₂|X)·f_{1|2}(Y₁|Y₂,X)) or SPML3 (Palmgren: marginals + log-OR association
  model). Include Y₂×X interactions to protect against misspecification.
- **AIPW hybrid (Sofer 2017)** combines IPW robustness with SPML efficiency — the
  recommended default when feasible (and the natural fit for causatr's AIPW once
  reweighted).
- **Reuse vs new**: weighted route = reuse (IPW). SPML2/3 have published author software
  (Lin & Zeng; Ghosh) — wrap where licensable, else implement focused versions; the AIPW
  hybrid composes with causatr's AIPW. Rubin's rules do NOT apply (not imputation-based).

## API design

```r
matcha(data, outcome = "y2",                       # secondary outcome
       exposure = "x", design = unmatched_cc(primary = "y1"),
       confounders = ~ z, estimator = "secondary_ipw")     # relative-prob weighted
matcha(..., estimator = "secondary_spml")          # SPML2/3
matcha(..., estimator = "secondary_aipw")          # Sofer hybrid (causatr AIPW)
```

## Support matrix

| Y₂ | Estimator | Engine | Variance | Status |
|---|---|---|---|---|
| binary | secondary_ipw | survey / causatr | sandwich | needs-test |
| binary | secondary_spml (SPML2/3) | matchatr / author code | info matrix | smoke |
| continuous | secondary_spml | matchatr | info matrix | smoke |
| binary/continuous | secondary_aipw | causatr AIPW | sandwich (DR) | needs-test |

## Implementation plan

- `R/secondary.R` — relative-probability weights from the primary-outcome sampling;
  `fit_secondary_ipw()` (weighted regression), `fit_secondary_spml()` (SPML2/3 joint
  likelihood), `fit_secondary_aipw()` (causatr AIPW on the reweighted sample).

## Variance / inference notes

Sandwich for the weighted/AIPW routes; information matrix for SPML. SPML is sensitive to
nuisance-model misspecification — the severity scales with the Y₁–Y₂ association strength;
include interactions and prefer the AIPW hybrid for robustness.

## Oracle testing strategy

- Author software (Lin & Zeng SPML2; Ghosh SPML3; Sofer RECSO/AIPW) as oracles where
  available; `survey`/causatr for the weighted/AIPW routes.
- Truth-based: simulate (Y₁, Y₂, X) with a known Y₂–X association, draw a CC sample on
  Y₁, confirm the naive (unweighted) Y₂ regression is biased while IPW/SPML/AIPW recover
  truth; show double robustness for AIPW.

## Chunk plan

1. Relative-probability weights + `secondary_ipw` + biased-naive-vs-corrected test.
2. `secondary_aipw` via causatr + double-robustness.
3. SPML2/3 (smoke / oracle where author code is available).

## Deferred items

Ma & Carroll (2016) nonparametric-robust SPML (high-dimensional); continuous-Y₂ SPML
beyond the basic case.
