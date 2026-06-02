# Phase 11 — Two-Phase Sampling Framework

> **Status: DESIGN**
> Book chapters: 12 (Multi-Phase Sampling).

## Scope

**In:** the unifying two-phase view of case-control designs — Phase-1 variables on all
units, Phase-2 (expensive) variables on a subsample drawn by a known/estimable scheme π.
Ascertainment-corrected likelihood and survey-weighted estimation of regression
parameters; the case-cohort / NCC / stratified CC as special cases. K≥2 phases noted.

**Out:** calibration of the Phase-2 weights (Phase 12), survival-specific two-phase Cox
(Phase 6/7 cover the design-weighted Cox; this phase is the general two-phase regression
machinery), full semiparametric MLE for missing covariates (Phase 14).

## Key design decisions

- **Three estimation routes** (Ch12): (1) known-π conditional score (efficient
  baseline); (2) Breslow-Cain two-step plug-in (estimate π̂ from Phase-1, plug into the
  conditional likelihood; = semiparametric MLE under a saturated π model); (3) joint
  semiparametric-efficient estimation with discrete Phase-1 variables. Plus the
  survey-weighted estimating equation Σ R_i/π_i ∂log f/∂β = 0.
- **"Estimated better than known"** — estimating π recovers Phase-1 information from
  unsampled units and *lowers* variance vs using known π (Ch12 §12.4.2). Default to the
  estimated-π two-step route; offer known-π.
- **Delegate to `survey::twophase` + `svyglm`** for the survey-weighted route; the
  conditional/joint routes need a focused likelihood implementation (or `osDesign`).
- **This phase generalizes Phases 6–7**: case-cohort and NCC are two-phase designs; keep
  the design objects compatible so a two-phase analysis can consume them.

## API design

```r
fit <- matcha(data, outcome = "case", exposure = "x",
              design = two_phase(phase1 = ~ z, phase2 = ~ x,
                                 strata = "samp_stratum", weights = "auto"),
              confounders = ~ z, estimator = "twophase")   # survey-weighted
```

## Support matrix

| Route | π | Engine | Variance | Status |
|---|---|---|---|---|
| survey-weighted | known/estimated | survey::twophase/svyglm | sandwich (Phase I+II) | needs-test |
| two-step plug-in | estimated | matchatr/osDesign | sandwich | needs-test |
| known-π conditional | known | matchatr | info matrix | needs-test |
| joint efficient | estimated (discrete Z) | matchatr | observed info | smoke |

## Implementation plan

- `R/two_phase.R` — `two_phase()` design constructor; `fit_twophase()` routing to
  `survey::twophase()`/`svyglm()` (survey route) or a conditional-likelihood fitter
  (known-π / two-step).
- `R/variance_twophase.R` — Phase-I + Phase-II variance decomposition (or read it off
  the `survey` fit).

## Variance / inference notes

Variance = Phase-I + design-dependent Phase-II component; sandwich for the
weighted/two-step routes, observed information for the joint route. The `survey` package
returns the correct two-phase variance directly.

## Oracle testing strategy

- `survey::twophase` + `svyglm` as engine + oracle; `osDesign` for two-phase
  case-control. `survival::cch` as the case-cohort special case (consistency check:
  the two-phase Cox should reproduce `cch`).
- Truth-based: simulate a cohort, draw a stratified two-phase sample, confirm β recovery
  and the "estimated > known" variance ordering.

## Chunk plan

1. `two_phase()` design + survey-weighted `fit_twophase()` + oracle.
2. Two-step plug-in (estimated π) + variance + "estimated-better-than-known" test.
3. Known-π conditional / joint route (smoke + special-case consistency with `cch`).

## Deferred items

Calibration (Phase 12), continuous Phase-1 variables in the joint route (curse of
dimensionality), full semiparametric MLE (Phase 14).
