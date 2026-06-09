# Phase 6 — Case-Cohort: Prentice, Self-Prentice, and Borgan Estimators

> **Status: Chunks 1–2 IMPLEMENTED (Prentice + Self-Prentice + LinYing + Borgan
> I/II IPW + stratified subcohort + design rejections). Chunk 3 DESIGN.**
> Book chapters: 16 (Overview), 17 (Survival Analysis / Sample Survey).

## Scope

**In:** case-cohort analysis — a fixed subcohort sampled once from the cohort serves as
controls at all failure times. Cox hazard-ratio estimation via the Prentice and
Self-Prentice pseudo-likelihoods and the Borgan I/II IPW weighted estimators; simple
and stratified subcohort sampling; IPW Breslow cumulative baseline hazard and absolute
risk F_x(t).

**Out:** nested case-control (Phase 5/7), calibration of weights (Phase 12), marginal
causal contrasts under sampling (Phase 10), additive-hazards model (defer / Phase 17
cross-ref).

## Key design decisions

- **Delegate to `survival::cch`**, which implements Prentice, Self-Prentice (`"LinYing"`
  / `"SelfPrentice"`), and Borgan I/II (`"I"`, `"II"`) with their correct asymptotic
  variances. matchatr supplies the design object (subcohort indicator, time, strata) and
  the contrast/absolute-risk layer; it does NOT hand-roll the pseudo-likelihood.
- **The pseudo-likelihood is NOT a true likelihood.** Controls (subcohort) are reused
  across all failure times → dependent score factors. SEs do NOT come from the naive
  information matrix and LR statistics are not χ² (Self & Prentice 1988). Variance =
  cohort information term + sampling-variation term. (Invariant in `hard-rules.md`.)
- **Borgan IPW (II) for stratified subcohorts uses plug-in asymptotic variance, not the
  robust sandwich** — the chapter warns the sandwich substantially overestimates it
  there. Default to the estimator-appropriate variance returned by `cch`.
- **Time scale is flexible at analysis** (a case-cohort advantage over NCC) — document.
- **Reject** a `case_cohort` design lacking a subcohort indicator or a time variable
  (`matchatr_bad_design`).

## API design

```r
fit <- matcha(nwtco, outcome = "rel", exposure = "histol",
              design = case_cohort(subcohort = "subcohort", time = "edrel",
                                   method = "LinYing"),
              confounders = ~ stage + age, estimator = "cch")
contrast(fit)                      # HR + CI (Self-Prentice variance)
absolute_risk(fit, x = newdata, times = c(2, 5))   # IPW Breslow F_x(t)
```

## Support matrix

| Subcohort | Method | Estimand | Variance | Status |
|---|---|---|---|---|
| simple | Prentice | HR | asymptotic (Self-Prentice) | needs-test |
| simple | Self-Prentice / LinYing | HR | asymptotic | needs-test |
| simple | Borgan I (IPW) | HR | asymptotic | ✅ nwtco oracle + truth DGP |
| stratified | Borgan II (IPW) | HR | plug-in asymptotic | ✅ nwtco oracle + truth DGP |
| simple | — | absolute risk F_x(t) | Breslow + delta | needs-test |
| case_cohort, no subcohort/time | — | — | ⛔ `matchatr_bad_design` |

## Implementation plan

- `R/case_cohort.R` — `fit_cch()` (build the `Surv(time, status) ~ x + cov` call +
  subcohort/id wiring for `survival::cch`, select method), HR contrast, design checks.
- `R/absolute_risk.R` — IPW Breslow cumulative baseline hazard Λ̂₀(t) and
  F̂_x(t) = 1 − exp(−exp(β̂·x)Λ̂₀(t)) with pointwise CIs (`survival`/`survey::svykm`
  cross-ref).

## Variance / inference notes

Use the variance `cch` returns per method (Self-Prentice asymptotic; Borgan plug-in).
Absolute risk: Breslow estimator variance is a Gaussian-process covariance (Ch19 §19.5.2
form) — pointwise CIs via the delta method on (β̂, Λ̂₀). Do NOT use the naive
information-matrix SE for the pseudo-likelihood β̂.

## Oracle testing strategy

- `survival::cch` on the `survival::nwtco` (Wilms tumour) dataset — the canonical
  case-cohort example; assert each method's β̂/SE matches.
- Truth-based: simulate a cohort with known Cox β, draw a subcohort, confirm Prentice /
  Self-Prentice / Borgan recover β and that full-cohort `coxph` agrees within sampling
  error.
- Absolute risk: compare F̂_x(t) to the full-cohort `survfit`/`basehaz` curve.

## Chunk plan

1. `fit_cch()` (Self-Prentice + Prentice) + nwtco oracle + design rejections.
2. Borgan I/II IPW + stratified subcohort + variance-method selection.
3. IPW Breslow absolute risk F_x(t) + pointwise CIs.

## Deferred items

Weight calibration (Phase 12), additive-hazards (`addhazard`) cross-ref (Phase 17),
marginal causal contrasts under case-cohort sampling (Phase 10), MI for unmeasured
covariates (Phase 13).
