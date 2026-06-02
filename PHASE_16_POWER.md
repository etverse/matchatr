# Phase 16 — Power and Sample Size

> **Status: DESIGN**
> Book chapters: 9 (Power and Sample Size).

## Scope

**In:** analytical sample-size / power calculation for unmatched and matched
case-control designs (Wald/score tests on the log OR), plus a simulation-based power
engine that captures confounding, misclassification, and missing data. Optional
extensions for two-phase and G×E designs.

**Out:** power for the survival sampling designs' HR (could be a later extension), full
Bayesian design.

## Key design decisions

- **Two routes**: (1) closed-form formulas — unmatched logistic (information-matrix
  based) and matched discordant-pair (McNemar) sample sizes; (2) a simulation engine
  (generate a large population, repeatedly draw the study sample, fit, test, count
  rejections) that handles confounding / misclassification / missingness the formulas
  cannot.
- **Simulation engine reuses the package's own fitters** (Phases 2–3) and the design
  samplers (Phase 5) — power is "draw a design, run `matcha`, test" in a loop, so it
  composes with everything matchatr can fit.
- **Delegate closed forms where good packages exist** (`epiR`, `samplesizeCaseControl`,
  `powerAnalysis` for G×E, `osDesign` for two-phase) and wrap them in a consistent
  `power_cc()` / `n_cc()` API; provide the simulation engine natively.

## API design

```r
n_cc(or = 2, p_exposed = 0.3, ratio = 2, power = 0.8, alpha = 0.05,
     design = "unmatched")                       # closed form
power_cc(n_cases = 200, or = 2, ..., design = "matched")
power_sim(dgp, design = nested_cc(...), estimator = "clogit",
          reps = 1000, target = "x")             # simulation
```

## Support matrix

| Design | Route | Engine | Status |
|---|---|---|---|
| unmatched | closed form | epiR / samplesizeCaseControl | needs-test |
| matched (discordant pairs) | closed form | native / epiR | needs-test |
| any matchatr-fittable | simulation | native + matcha loop | needs-test |
| G×E, two-phase | closed form | powerAnalysis / osDesign | smoke |

## Implementation plan

- `R/power.R` — `n_cc()`, `power_cc()` (closed forms + wrappers), `power_sim()`
  (simulation loop over a user DGP + design + estimator, returning empirical power and a
  curve).

## Variance / inference notes

Closed forms use the information-matrix / discordant-pair variance under the alternative.
Simulation reports empirical power with a Monte-Carlo CI; document reps needed for a given
precision.

## Oracle testing strategy

- `epiR::epi.sscc` and `samplesizeCaseControl` as oracles for the closed forms.
- Self-consistency: `power_sim` empirical power should match the closed-form power for the
  simple no-confounding logistic case (validation of the simulation engine against the
  analytic formula).

## Chunk plan

1. Closed-form `n_cc()` / `power_cc()` (unmatched + matched) + `epiR` oracle.
2. `power_sim()` simulation engine + analytic-vs-empirical consistency test.
3. G×E / two-phase wrappers (smoke).

## Deferred items

Power for survival-sampling HR; sequential/adaptive designs; FDR-based multiple-testing
power.
