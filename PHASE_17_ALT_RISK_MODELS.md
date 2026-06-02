# Phase 17 — Alternative Risk Models: Additive, Multiplicative, and Mixed

> **Status: DESIGN**
> Book chapters: 11 (Alternative Formulation of Models), with 17 (additive hazards).

## Scope

**In:** odds-ratio / risk models beyond the standard log-linear (multiplicative) form —
linear (additive) odds-ratio models, exponential-additive, and the Breslow-Storer /
Thomas / Barlow mixed families that nest additive and multiplicative via a mixing
parameter; identity-link binomial (additive risk) and the lexpit (linear+logistic) model;
score tests of departure from multiplicativity. Cross-reference the additive-hazards
survival model (`addhazard`) for the time-to-event designs.

**Out:** the default multiplicative logistic/Cox (Phases 2–7).

## Key design decisions

- **Conditional logistic with a non-log-linear odds-ratio form**: linear (1 + xᵀβ) or
  exponential-additive (1 + Σ(e^{β_k}−1)x_k). Fit by maximizing the conditional
  likelihood with the chosen OR form.
- **Breslow-Storer λ family** (unifying additive↔multiplicative via λ) and the Thomas α
  / Barlow covariate-specific α_k mixed models — offer the unified family + a score test
  of H₀ (multiplicative) vs additive.
- **Additive-risk (identity-link) binomial** via `blm`; **lexpit** (Kovalchik 2013) mixed
  linear+logistic via `blm`. With known sampling fractions these recover cohort-level
  absolute-risk parameters.
- **Inference cautions**: for the linear OR model Wald is unreliable — prefer score/LR;
  the exponential-additive form has better-behaved Wald SEs; the constrained identity-link
  model needs robust (influence-function) variance.
- **Reuse** `survival::clogit` (conditional, non-log-linear via custom likelihood), `blm`
  (identity/lexpit), and `addhazard` (additive hazards). Barlow's covariate-specific mixed
  model needs custom code.

## API design

```r
matcha(data, outcome = "case", exposure = "x", design = matched_cc(strata = "set"),
       estimator = "clogit", or_form = "additive")        # linear OR
matcha(..., or_form = "breslow_storer")                    # mixed family + test
matcha(data, ..., design = unmatched_cc(), estimator = "additive_risk")  # blm identity link
```

## Support matrix

| Model | Design | Engine | Inference | Status |
|---|---|---|---|---|
| linear (additive) OR | matched/unmatched | clogit (custom lik) | score/LR | needs-test |
| exponential-additive OR | matched | clogit (custom) | Wald/LR | needs-test |
| Breslow-Storer mixed (λ) | matched | custom | LR + multiplicativity test | smoke |
| identity-link additive risk | unmatched (known fractions) | blm | robust | needs-test |
| lexpit (linear+logistic) | unmatched | blm | robust | smoke |
| additive hazards | NCC/case-cohort | addhazard | sandwich | cross-ref Phase 7 |

## Implementation plan

- `R/alt_risk_models.R` — custom conditional-likelihood objective for linear /
  exponential-additive OR forms; Breslow-Storer λ family + score test; `blm` wrappers for
  identity-link / lexpit; `addhazard` cross-ref.

## Variance / inference notes

Score/LR for the linear OR (Wald unreliable); robust variance for the constrained
identity-link model; standard for exponential-additive. Document the model-selection
guidance: choose scale by prior biology, not data-dredging.

## Oracle testing strategy

- `blm` (identity-link / lexpit) and `addhazard` (additive hazards) as engines + oracles.
- Truth-based: simulate under an additive-risk DGP; confirm the additive model recovers
  the true risk difference where the multiplicative model is misspecified, and the
  multiplicativity score test rejects.

## Chunk plan

1. Linear / exponential-additive conditional OR forms + score-test inference.
2. Identity-link additive risk + lexpit via `blm` + robust variance.
3. Breslow-Storer mixed family + multiplicativity test (smoke); `addhazard` cross-ref.

## Deferred items

Barlow covariate-specific mixed models (custom, low priority); broader hazard-model
families.
