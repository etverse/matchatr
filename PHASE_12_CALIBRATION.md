# Phase 12 — Calibration of Sampling Weights

> **Status: DESIGN**
> Book chapters: 13 (Calibration), with 17, 19 (survival calibration).

## Scope

**In:** survey calibration of the case-control / design / inclusion-probability weights
(Phases 7, 9, 11) to known Phase-1 totals, for efficiency gains — raking / regression
calibration on auxiliary variables built from influence-function surrogates; the
imputation-based (Han 2015) doubly-robust calibration variant.

**Out:** the base weights themselves (Phases 7, 9), multiple imputation as a standalone
analysis (Phase 13 — though MI-for-calibration is referenced here).

## Key design decisions

- **Calibrate weights so weighted auxiliary totals match cohort totals**, minimizing
  distance to the base weights (Deville-Särndal). Efficiency gain comes from projecting
  the estimator onto the auxiliary space.
- **Auxiliary variables should be influence-function surrogates / Phase-1 predictors of
  the expensive covariate or outcome** — NOT variables already in the model (those give
  no gain; the model already controls them). This is the key practical rule (Ch13).
- **Delegate to `survey::calibrate()`** on the `twophase`/`svydesign` object; matchatr
  builds the auxiliary variables (from a Phase-1 model or imputed X) and wires them in.
- **MI-for-calibration (Han 2015)** is the most efficient IPW variant and is
  doubly-robust; implement via `mice`/`smcfcs` to construct the auxiliary, then calibrate.
- This phase is **largely reuse** (survey-sampling standard + MI infrastructure).

## API design

```r
fit <- matcha(data, outcome = "case", exposure = "x",
              design = nested_cc(strata = "set", time = "t",
                                 weights = ncc_weights("km"),
                                 calibrate = ~ surrogate_dose + ageRx),
              estimator = "ipw_cox")
```

## Support matrix

| Base weight | Auxiliary | Engine | Variance | Status |
|---|---|---|---|---|
| NCC KM | influence-fn surrogate | survey::calibrate | sandwich | needs-test |
| case-cohort Borgan | Phase-1 predictor | survey::calibrate | sandwich | needs-test |
| CCW q₀ | imputed-X (Han MI) | survey::calibrate + mice | sandwich / Rubin | needs-test |
| calibrate on in-model var | — | — | (warn: no gain) `matchatr_calibration_no_gain` |

## Implementation plan

- `R/calibrate_weights.R` — `calibrate_spec()` parsing; build auxiliary variables
  (Phase-1 model influence surrogates or imputed X); call `survey::calibrate()`;
  return calibrated weights into the existing estimator path.

## Variance / inference notes

Sandwich variance from the calibrated `survey` design (calibration reduces the Phase-II
variance via projection). MI-for-calibration combines with Rubin's rules. Cautions from
Ch13/17/19: calibrating on model variables yields no gain (warn).

## Oracle testing strategy

- `survey::calibrate` engine + oracle. Truth-based: simulate a cohort with a strong
  Phase-1 surrogate of the expensive covariate; confirm calibration shrinks the SE
  relative to the uncalibrated Phase-7 fit while keeping β unbiased (handbook §19.6,
  §17.4.3 breast-cancer dose example).

## Chunk plan

1. `calibrate_spec()` + survey::calibrate wiring on NCC KM weights + SE-reduction test.
2. Case-cohort calibration + no-gain warning for in-model auxiliaries.
3. MI-for-calibration (Han 2015) doubly-robust variant.

## Deferred items

Continuous Phase-1 auxiliaries; calibration for counter-matched NCC (Rivera & Lumley);
calibration of the CCW q₀ weights for marginal contrasts.
