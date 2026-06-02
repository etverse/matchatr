# Phase 1 — Design Taxonomy, Data Model, and Two-Step API

> **Status: IMPLEMENTED** (2026-06-02)
> Book chapters: 2 (Design Issues), with the sampling structures of 4, 16.
> Foundation phase — fixes the data model and API every later phase depends on.
> Code: `R/cc_design.R`, `R/matcha.R`, `R/dispatch.R`, `R/contrast.R`,
> `R/constructors.R`, `R/checks.R`, `R/print.R`. Tests: `test-cc_design.R`,
> `test-matcha.R`, `test-dispatch.R`, `test-contrast.R`, `test-constructors.R`,
> `test-print.R`, `test-rejections.R`.

## Scope

**In:** the unified `matchatr_design` S3 object and its constructors; the two-step
`matcha()` + `contrast()` API skeleton; the orthogonal `design =` / `estimator =`
axes; where weights live; the dispatch table from (design, estimator) → engine.

**Out:** any actual estimator (Phases 2+). No numeric results in this phase — it is
plumbing, validation, and dispatch only.

## Key design decisions

- **Two orthogonal axes: `design =` and `estimator =`.** The *design* object encodes
  the sampling structure (how cases/controls/subcohort were drawn): strata, matching
  ratio, time scale, prevalence q₀, inclusion probabilities/weights. The *estimator*
  decides the analysis (conditional vs marginal; OR vs HR vs RD/RR). This mirrors
  causatr's `estimator =` choice and avoids overloading one argument.
- **`matcha()` is the fit verb** (mirrors `causat()`/`surv_fit()`); `contrast()` is
  reused from the etverse convention. `matcha()` → `matchatr_fit`; `contrast()` →
  `matchatr_result`.
- **Weights are NEVER a data column.** Case-control weights (q₀-based) and design
  weights (inclusion-probability) are computed by the design layer and stored on the
  fit (`fit$details$weights`), then passed to the engine as observation weights.
  Mirrors causatr's invariant. (See `hard-rules.md`.)
- **Case-control weights and design weights are distinct slots**, not a single
  `weights` field — they have different variance consequences and a fit may, in
  principle, need to know which kind it holds.
- **Rejected:** a single `cc_data()` wrapper that mutates the data frame to carry
  design metadata. Instead the design object is passed alongside the untouched data,
  so the same data frame can be analysed under different designs.

## Design taxonomy (constructors)

| Constructor | Sampling structure | Primary estimators | Chapter |
|---|---|---|---|
| `unmatched_cc(prevalence=)` | independent CC sample | logistic OR, MH; CCW marginal | 3 |
| `matched_cc(strata=, ratio=)` | individually/freq matched | conditional logistic (clogit) | 4 |
| `nested_cc(strata=, time=, ...)` | risk-set sampling, m:1 | clogit/coxph strata; IPW Cox | 16,18,19 |
| `case_cohort(subcohort=, time=)` | fixed subcohort | Prentice/Self-Prentice/Borgan (cch) | 16,17 |
| `two_phase(phase1=, phase2=, ...)` | two-phase sample | survey / calibration | 12,13 |
| `counter_matched(strata=, ...)` | stratified risk-set | weighted (counter-match) Cox | 16,19 |

Each constructor returns a `matchatr_design` carrying: `type`, strata/time columns,
matching ratio, `prevalence` (q₀, optional), and a `weight_spec` describing how
inclusion/CC weights are to be computed (deferred to the weight phases).

## API design (concrete)

```r
# Matched CC -> conditional OR (Phase 3)
matcha(data, outcome = "case", exposure = "x",
       design = matched_cc(strata = "set", ratio = 2),
       confounders = ~ age + smoke, estimator = "clogit")

# NCC -> risk-set HR via conditional partial likelihood (Phase 5)
matcha(data, outcome = "case", exposure = "x",
       design = nested_cc(strata = "set", time = "t"), estimator = "clogit")

# Marginal causal RD from an unmatched CC sample (Phase 9, Rose & van der Laan)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(prevalence = 0.02),
              confounders = ~ age + smoke, estimator = "ccw_gformula")
contrast(fit, type = "difference", ci_method = "sandwich")
```

## Support matrix (this phase)

| (design, estimator) | Resolves to | Status |
|---|---|---|
| any design | unknown estimator string | ⛔ classed error `matchatr_bad_estimator` |
| `matched_cc`/`nested_cc` | `clogit` | dispatch stub → Phase 3/5 |
| `unmatched_cc` | `logistic`/`mh` | dispatch stub → Phase 2 |
| `case_cohort` | `cch` | dispatch stub → Phase 6 |
| any | `ccw_*` without `prevalence` | ⛔ `matchatr_missing_prevalence` |
| binary-only designs | non-binary `outcome` | ⛔ `matchatr_bad_outcome` |

## Implementation plan (per file)

- `R/cc_design.R` — the six constructors + `new_matchatr_design()` + `print` method +
  validators (strata exist, ratio sane, q₀ ∈ (0,1)).
- `R/matcha.R` — `matcha()`: validate inputs, build the design, resolve
  (design, estimator) via a dispatch table, call the (stubbed) engine, return
  `new_matchatr_fit()`.
- `R/constructors.R` — `new_matchatr_fit()`, `new_matchatr_result()`.
- `R/checks.R` — shared validators + the classed-error helpers.
- `R/print.R` — `print.matchatr_design`, `print.matchatr_fit`.

## Variance / inference notes

None in this phase. The fit object reserves a `details$variance_kind` slot so later
phases record which correction applies (model-info / sandwich / Self-Prentice /
Samuelsen / CCW-IF / bootstrap).

## Oracle testing strategy

No numeric oracle (no estimation yet). Tests assert: constructors build correct
objects; dispatch routes (design, estimator) pairs correctly; every rejection path
fires its classed error (`expect_snapshot(error = TRUE)`).

## Chunk plan

1. ✅ Design constructors + validators + print.
2. ✅ `matcha()` skeleton + dispatch table + fit constructor.
3. ✅ Rejection paths + tests.
4. ✅ `contrast()` skeleton (signature + `matchatr_result` contract; aborts
   `matchatr_not_estimated` until estimation lands).

## Deferred items (each owned by a planned phase)

- **All estimation** → Phases 2–6 (logistic/MH → P2; conditional logistic → P3;
  polytomous → P4; risk-set / weighted Cox → P5; case-cohort `cch` → P6).
- **Weight computation** (the `weight_spec` is declared here, realised later):
  inclusion-probability (Samuelsen/Borgan) → P7; case-control q₀ weights → P8/P9;
  counter-matching weights → P5/P7; survey / calibration design weights → P11/P12.
- **`contrast()` estimation body** (the skeleton lands here) → the causal-contrast
  phases P9 (CCW g-formula/IPW/AIPW/TMLE) and P10 (design-weighted survival).
- **Two-phase / counter-matched engines** (`survey_twophase`, `weighted_cox`) →
  P11/P12 and P5/P7 respectively.
- **Not planned (intentional non-feature):** supplying a prevalence q₀ on
  non-`unmatched_cc` designs. CCW (Rose & van der Laan) is an unmatched-CC
  technique; matched/nested/case-cohort designs use design-based inclusion
  weighting instead. Only `unmatched_cc(prevalence=)` carries q₀, so CCW on other
  designs is correctly gated by `matchatr_missing_prevalence`.
