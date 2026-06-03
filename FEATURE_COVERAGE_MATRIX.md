# Feature Coverage Matrix

Single source of truth for what works, what's tested, and at what fidelity.
**Every PR that changes a feature MUST update this file.**

> **Status: classical odds-ratio engines landing.** On top of the PHASE_1
> design objects, `matcha()` fit verb, `(design, estimator)` dispatch, and
> input-validation / rejection paths, the unmatched case-control **logistic /
> Mantel-Haenszel** ORs (PHASE_2) and the matched case-control **conditional
> logistic** OR (PHASE_3 Chunk 1, `survival::clogit`) run end to end through the
> shared `contrast()` / `tidy()` OR layer. The remaining estimator cells stay
> pending until their phases land.

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

**Chunks 1–3 implemented — the unmatched case-control layer is complete.**
`matcha(estimator = "logistic")` fits `stats::glm(family = binomial)` (or a
pluggable `model_fn`, e.g. `mgcv::gam`) and `estimator = "mh"` computes the
Mantel-Haenszel stratified OR; `contrast(type = "or")` reports the exposure
conditional / summary odds ratio(s) with a Wald interval, `tidy()` / `summary()`
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
| binary, stratified | mh | summary OR | OR | RBG | ✅ vs `stats::mantelhaen.test` (OR + CI) | `test-mantel_haenszel.R` |
| binary, crude (no strata) | mh | OR | OR | RBG | ✅ vs closed-form 2×2 | `test-mantel_haenszel.R` |
| non-binary exposure | mh | — | — | — | ⛔ `matchatr_bad_input` | `test-mantel_haenszel.R` |
| zero-margin / sandwich·bootstrap CI | mh | — | — | — | ⛔ `matchatr_unestimable_exposure` / `matchatr_unsupported_variance` | `test-mantel_haenszel.R` |

S3 surface: `tidy.matchatr_fit` (broom-style coefficient / OR table, model or
`robust` SE), `tidy.matchatr_result`, `summary.matchatr_fit`,
`print.matchatr_result` — all tested in `test-unconditional.R`. Smooth-of-exposure
(spline OR-curve) is deferred (the OR is then a value-vs-value contrast).

## Matched case-control (PHASE_3)

**Chunk 1 implemented — conditional logistic regression.**
`matcha(design = matched_cc(strata = ...), estimator = "clogit")` fits the
matched case-control conditional maximum likelihood via `survival::clogit`
(`outcome ~ exposure + confounders + strata(set)`, each matched set a stratum),
and `contrast(type = "or")` reports the exposure's conditional odds ratio with a
partial-likelihood-information Wald interval. The matching variables are
conditioned away (no estimable coefficient); only the exposure / adjustment ORs
are reported. McNemar 1:1 closed form (Chunk 2) and effect modification /
variable-ratio handling (Chunk 3) remain pending.

| Matching | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|
| 1:1, binary x | clogit | cond. OR | OR | partial-lik info | ✅ closed-form McNemar: OR = n10/n01 **and** Var(log OR) = 1/n10+1/n01 (independent of clogit) | `test-clogit.R` |
| M:1, binary x | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP (CMLE recovers β within 3.5 SE) + `survival::clogit` pass-through | `test-clogit.R` |
| variable ratio (mixed 1:1/1:2/1:3) | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP | `test-clogit.R` |
| continuous exposure (per unit) | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP | `test-clogit.R` |
| infert induced/spontaneous | clogit | cond. OR | OR | partial-lik info | 🟡 regression pin vs canonical `survival::clogit` example (OR ≈ 4.09, 7.29) | `test-clogit.R` |
| M:1 + non-matching covariate | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP (recovers adjusted β) + `survival::clogit` pass-through | `test-clogit.R` |
| factor exposure (per-level) | clogit | cond. OR per level | OR | partial-lik info | ✅ vs `survival::clogit` (+ reference) | `test-clogit.R` |
| multi-column strata (frequency matching) | clogit | cond. OR | OR | partial-lik info | ✅ vs `survival::clogit` (crossed strata) | `test-clogit.R` |
| clogit | — | RD / RR | — | — | ⛔ `matchatr_unidentified_estimand` | `test-clogit.R` |
| clogit OR | — | OR | — | sandwich / bootstrap | ⛔ `matchatr_unsupported_variance` | `test-clogit.R` |
| exposure constant within strata | clogit | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-clogit.R` |
| missing exposure / confounder | clogit | cond. OR | OR | partial-lik info | ⚠️ `matchatr_dropped_rows` (complete-case) | `test-clogit.R` |

The conditional OR assembly is shared with the unmatched logistic engine via
`conditional_or_result()` (exposure coefficient by term position, Wald interval
on the log scale, exponentiated). `tidy()` renders the per-term OR table (no
intercept row). McNemar (1:1 closed form), effect modification across strata, and
explicit variable-ratio handling are pending Chunks 2–3.

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
