# Phase 13 — Multiple Imputation for the Missing-by-Design Covariate

> **Status: DESIGN**
> Book chapters: 20 (Multiple Imputation for Sampled Cohort Data).

## Scope

**In:** multiple imputation of the expensive covariate X (missing by design for unsampled
cohort members) followed by full-cohort Cox regression and Rubin pooling — the
White-Royston approximate method and the Bartlett rejection-sampling (`smcfcs`) method;
auxiliary variables; competing risks / left-truncation / stratified Cox extensions.

**Out:** MI for missing confounders in non-survival CC (composes with `causatr::causat_mice`
— note the cross-package reuse), calibration-via-MI (Phase 12).

## Key design decisions

- **MAR by design**: X is missing only for unsampled units; the imputation model is the
  conditional X | (Z, δ, Λ̂₀(t)) implied by the Cox likelihood. Two methods:
  - **Approximate (White & Royston 2009)** — regress X on (Z, δ, Nelson-Aalen Λ̂₀(t)),
    draw K parameter vectors, impute. Fast; minor downward bias for large effects. Via
    `mice`.
  - **Rejection sampling (Bartlett et al. 2015)** — sample X from a proposal, accept via
    the Cox-likelihood rule; unbiased, handles nonlinear terms. Via `smcfcs`.
- **Include auxiliary Phase-1 predictors of X** (e.g. an FFQ for diet) for large
  efficiency gains.
- **Delegate imputation to `mice` / `smcfcs`** and pooling to Rubin's rules; matchatr
  supplies the Nelson-Aalen term and the design-aware imputation setup, then fits the
  full-cohort Cox per imputation. **Reuse causatr's pooling** (`pool_rubin`) where it
  applies, for consistency across the etverse.
- **Do not impute the outcome / event** — use Y, δ as predictors in the imputation model
  (mirrors causatr's MI invariant).

## API design

```r
fit <- matcha(cohort, outcome = "case", exposure = "x",   # x measured on sample only
              design = case_cohort(subcohort = "sub", time = "t"),
              confounders = ~ z, estimator = "mi_cox",
              impute = mi_spec(method = "smcfcs", m = 20, aux = ~ ffq))
contrast(fit)    # pooled HR (Rubin)
```

## Support matrix

| Method | Design | Engine | Pooling | Status |
|---|---|---|---|---|
| approximate (W&R) | NCC / case-cohort | mice + survival::coxph | Rubin | needs-test |
| rejection (Bartlett) | NCC / case-cohort | smcfcs | Rubin | needs-test |
| + auxiliary | both | mice/smcfcs | Rubin | needs-test |
| competing risks / truncation | both | smcfcs | Rubin | smoke |
| impute the outcome | — | — | ⛔ `matchatr_impute_outcome` |

## Implementation plan

- `R/mi_cox.R` — `mi_spec()`; build the Nelson-Aalen term; run `mice`/`smcfcs`; loop
  full-cohort `coxph` over completed datasets; pool via Rubin (reuse `causatr::pool_rubin`
  if exported, else local). Reject outcome imputation.

## Variance / inference notes

Rubin's rules: total variance W̄ + (1+1/K)B̄. The approximate method has a small downward
bias for large effects (document); rejection sampling is unbiased. CIs use the
Barnard-Rubin df.

## Oracle testing strategy

- `mice` / `smcfcs` engines + oracle. Truth-based: simulate a cohort with known Cox β,
  set X missing for unsampled units, impute + pool, confirm β recovery and that an
  informative auxiliary improves efficiency (handbook §20.6.2).
- Compare to the IPW (Phase 7) and full-cohort estimates on the same DGP — MI should be
  at least as efficient when the imputation model is good.

## Chunk plan

1. Approximate MI (`mice` + Nelson-Aalen) + Rubin pooling + truth oracle + outcome-impute
   rejection.
2. Rejection sampling via `smcfcs` + auxiliary variables.
3. Competing-risks / left-truncation / stratified extensions (smoke).

## Deferred items

MI for non-survival CC confounders (defer to `causatr::causat_mice`), MI-for-calibration
(Phase 12), full semiparametric MLE alternative (Phase 14).
