# Feature Coverage Matrix

Single source of truth for what works, what's tested, and at what fidelity.
**Every PR that changes a feature MUST update this file.**

> **Status: first estimator landed (PHASE_2 Chunk 1).** On top of the PHASE_1
> design objects, `matcha()` fit verb, `(design, estimator)` dispatch, and
> input-validation / rejection paths, the unmatched case-control **logistic
> conditional OR** now runs end to end (`stats::glm` + the `contrast()` /
> `tidy()` / `summary()` OR layer). The remaining estimator cells stay pending
> until their phases land.

## Legend

| Symbol | Meaning |
|---|---|
| ✅ | Truth-based: estimate/SE/CI checked against analytical truth or external reference |
| 🟡 | Smoke: runs without error, finite output, target not pinned |
| ❌ | No test |
| ⛔ | Rejected by design (rejection path tested) |

References (planned): `survival::clogit` / `survival::cch`, `multipleNCC`, Hernán &
Robins / handbook book values, closed-form analytical truth, `causatr` / `survatr` on
the explicitly reweighted pseudo-cohort, R `tmle` (CCW-TMLE), `delicatessen`.

Columns: Design × Estimator × Weight × Outcome × Estimand × Contrast × Variance × Status × Test.

---

## Design & API layer (PHASE_1)

Plumbing, not estimation — fidelity symbols above (which grade estimate/SE/CI)
do not apply, so this section reports structural coverage instead.

| Feature | Status | Test |
|---|---|---|
| `unmatched_cc` / `matched_cc` / `nested_cc` / `case_cohort` / `two_phase` / `counter_matched` build a valid `matchatr_design` | ✅ built + asserted | `test-cc_design.R` |
| Constructor validation: q0 ∈ (0,1), ratio whole ≥ 1, strata non-empty character | ⛔ rejection tested | `test-cc_design.R`, `test-rejections.R` |
| `matcha()` returns a `matchatr_fit` (data.table copy, no mutation); runs the resolved engine (logistic populates `model`, unwired engines leave it `NULL`) | ✅ asserted | `test-matcha.R`, `test-unconditional.R` |
| `(design, estimator)` → engine dispatch; CCW family valid on any design | ✅ routing pinned | `test-dispatch.R` |
| Binary-outcome resolution (logical / 2-level factor / numeric 0/1) | ✅ + ⛔ | `test-matcha.R` |
| Reject unknown / design-incompatible estimator (`matchatr_bad_estimator`) | ⛔ | `test-dispatch.R`, `test-rejections.R` |
| Reject non-binary outcome (`matchatr_bad_outcome`) | ⛔ | `test-matcha.R`, `test-rejections.R` |
| Reject CCW without q0 (`matchatr_missing_prevalence`) | ⛔ | `test-matcha.R`, `test-rejections.R` |
| Reject missing columns / wrong design object (`matchatr_bad_design`) | ⛔ | `test-matcha.R`, `test-rejections.R` |
| Warn on uninformative conditional-likelihood strata (`matchatr_uninformative_stratum`) | ⚠️ warn tested | `test-matcha.R`, `test-rejections.R` |
| Reject column with two roles (outcome/exposure vs covariate/design; `matchatr_bad_input`) | ⛔ | `test-matcha.R` |
| Reject duplicated `data` column names (`matchatr_bad_input`) | ⛔ | `test-matcha.R` |
| `contrast()` verb skeleton: signature + `matchatr_result` contract; `matchatr_not_estimated` until estimation | ⛔ | `test-contrast.R` |
| `print.matchatr_design`, `print.matchatr_fit` | ✅ snapshot | `test-print.R` |

No estimator engine runs yet; no numeric oracle applies (per PHASE_1 design).

## Unmatched case-control (PHASE_2)

**Chunks 1–2 (logistic conditional OR, all exposure types) implemented.**
`matcha(estimator = "logistic")` fits `stats::glm(family = binomial)` (or a
pluggable `model_fn`, e.g. `mgcv::gam`); `contrast(type = "or")` reports the
exposure conditional odds ratio(s) with a Wald interval, `tidy()` / `summary()`
render the OR table, and RD / RR are rejected as unidentified without q0.

| Exposure | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|
| binary | logistic | cond. OR | OR | model | ✅ truth DGP + `glm` + 2×2 Woolf | `test-unconditional.R` |
| binary | logistic | cond. OR | OR | sandwich | ✅ vs `sandwich::sandwich` | `test-unconditional.R` |
| two-level factor | logistic | cond. OR | OR | model | ✅ == 0/1 coding | `test-unconditional.R` |
| continuous | logistic | cond. OR (per unit) | OR | model | ✅ vs `glm` | `test-unconditional.R` |
| categorical k>2 | logistic | cond. OR per level | OR | model | ✅ vs `glm`; `esoph` book oracle | `test-unconditional.R` |
| ordinal (numeric score) | logistic | cond. OR / trend | OR | model | ✅ vs `glm` | `test-unconditional.R` |
| continuous / smooth confounder | logistic (GAM via `model_fn`) | cond. OR | OR | model/sandwich | ✅ == `glm` (linear) + 🟡 smooth | `test-unconditional.R` |
| logistic | — | RD / RR | — | — | ⛔ `matchatr_unidentified_estimand` | `test-unconditional.R` |
| logistic OR | — | OR | — | bootstrap | ⛔ `matchatr_unsupported_variance` | `test-unconditional.R` |
| constant / collinear exposure | logistic | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-unconditional.R` |
| ordered-factor exposure | logistic | — | — | — | ⛔ `matchatr_bad_input` (polynomial contrasts) | `test-unconditional.R` |
| binary, stratified | mh | cond. OR | OR | RBG | ❌ Chunk 3 | — |

S3 surface: `tidy.matchatr_fit` (broom-style coefficient / OR table, model or
`robust` SE), `tidy.matchatr_result`, `summary.matchatr_fit`,
`print.matchatr_result` — all tested in `test-unconditional.R`. Smooth-of-exposure
(spline OR-curve) is deferred (the OR is then a value-vs-value contrast).

## Matched case-control (PHASE_3)

_Pending implementation._

## Multiple case / control groups (PHASE_4)

_Pending implementation._

## Nested case-control (PHASE_5)

_Pending implementation._

## Case-cohort (PHASE_6)

_Pending implementation._

## IPW for nested case-control (PHASE_7)

_Pending implementation._

## Case-control-weighted causal contrasts (PHASE_9)

_Pending implementation._

## Design-weighted causal survival (PHASE_10)

_Pending implementation._

## Two-phase / calibration (PHASE_11, PHASE_12)

_Pending implementation._

## Multiple imputation / semiparametric MLE (PHASE_13, PHASE_14)

_Pending implementation._

## Small-sample / power / alternative models / secondary analysis (PHASE_15–18)

_Pending implementation._

## Extensions: SCCS, response-selective (PHASE_19, PHASE_20)

_Pending implementation._
