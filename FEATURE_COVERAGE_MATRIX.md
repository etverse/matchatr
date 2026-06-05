# Feature Coverage Matrix

Single source of truth for what works, what's tested, and at what fidelity.
**Every PR that changes a feature MUST update this file.**

> **Status: classical odds-ratio engines landing.** On top of the PHASE_1
> design objects, `matcha()` fit verb, `(design, estimator)` dispatch, and
> input-validation / rejection paths, the unmatched case-control **logistic /
> Mantel-Haenszel** ORs (PHASE_2) and the matched case-control **conditional
> logistic** OR (PHASE_3 Chunk 1, `survival::clogit`) and the **polytomous**
> subtype ORs for multi-group outcomes (PHASE_4 Chunk 1, `nnet::multinom`) run
> end to end through the shared `contrast()` / `tidy()` OR layer. The remaining
> estimator cells stay pending until their phases land.

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

**Chunks 1–3 implemented — the matched case-control layer is complete.**
`matcha(design = matched_cc(strata = ...), estimator = "clogit")` fits the
matched case-control conditional maximum likelihood via `survival::clogit`
(`outcome ~ exposure + confounders + strata(set)`, each matched set a stratum),
and `contrast(type = "or")` reports the exposure's conditional odds ratio with a
partial-likelihood-information Wald interval. `estimator = "mcnemar"` computes the
1:1 matched-pair OR = n10/n01 with Var(log OR) = 1/n10 + 1/n01 in closed form
(no `clogit`), rejecting M:1 / richer matching toward `clogit`. Adding
`effect_modifier = "m"` fits `outcome ~ exposure * m + ... + strata(set)` and
reports the exposure's **stratum-specific** conditional OR within each modifier
level (`beta_x` at the reference level, `beta_x + beta_{x:level}` elsewhere,
Wald intervals from the joint partial-likelihood variance). The matching
variables are conditioned away (no estimable coefficient); only the exposure /
adjustment / interaction ORs are reported.

| Matching | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|
| 1:1, binary x | clogit | cond. OR | OR | partial-lik info | ✅ closed-form McNemar: OR = n10/n01 **and** Var(log OR) = 1/n10+1/n01 (independent of clogit) | `test-clogit.R` |
| 1:1, binary x | mcnemar | cond. OR | OR | McNemar 1/n10+1/n01 | ✅ closed form (hand-counted) **and** `survival::clogit` exact equality **and** truth DGP | `test-mcnemar.R` |
| M:1, binary x | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP (CMLE recovers β within 3.5 SE) + `survival::clogit` pass-through | `test-clogit.R` |
| variable ratio (mixed 1:1/1:2/1:3) | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP | `test-clogit.R` |
| continuous exposure (per unit) | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP | `test-clogit.R` |
| infert induced/spontaneous | clogit | cond. OR | OR | partial-lik info | 🟡 regression pin vs canonical `survival::clogit` example (OR ≈ 4.09, 7.29) | `test-clogit.R` |
| M:1 + non-matching covariate | clogit | cond. OR | OR | partial-lik info | ✅ truth DGP (recovers adjusted β) + `survival::clogit` pass-through | `test-clogit.R` |
| factor exposure (per-level) | clogit | cond. OR per level | OR | partial-lik info | ✅ vs `survival::clogit` (+ reference) | `test-clogit.R` |
| multi-column strata (frequency matching) | clogit | cond. OR | OR | partial-lik info | ✅ vs `survival::clogit` (crossed strata) | `test-clogit.R` |
| effect modifier (`x:m`), modifier = matching var (1:1) | clogit | stratum-specific OR | OR | partial-lik info | ✅ within-level McNemar: OR **and** Var(log OR) per level (independent of clogit) | `test-effect_modification.R` |
| effect modifier (`x:m`), multi-level + adjustment (M:1) | clogit | stratum-specific OR | OR | partial-lik info | ✅ hand-built `survival::clogit` linear combos (point/SE/CI) + truth DGP | `test-effect_modification.R` |
| effect modifier, character / logical coercion | clogit | stratum-specific OR | OR | partial-lik info | ✅ == factor coding | `test-effect_modification.R` |
| effect modifier on non-clogit engine | — | — | — | — | ⛔ `matchatr_bad_input` | `test-effect_modification.R` |
| continuous / numeric effect modifier | clogit | — | — | — | ⛔ `matchatr_bad_input` (bin / `factor()`) | `test-effect_modification.R` |
| effect modifier = outcome / exposure | clogit | — | — | — | ⛔ `matchatr_bad_input` | `test-effect_modification.R` |
| effect modifier + 3+-level factor exposure | clogit | — | — | — | ⛔ `matchatr_unsupported_combination` | `test-effect_modification.R` |
| effect modifier, aliased per-level interaction | clogit | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-effect_modification.R` |
| unconditional 1:1 MLE (OR² bias) | — | — | — | — | ✅ pins invariant: unconditional β = 2 × conditional (McNemar) β | `test-mcnemar.R` |
| clogit | — | RD / RR | — | — | ⛔ `matchatr_unidentified_estimand` | `test-clogit.R` |
| clogit OR | — | OR | — | sandwich / bootstrap | ⛔ `matchatr_unsupported_variance` | `test-clogit.R` |
| exposure constant within strata | clogit | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-clogit.R` |
| missing exposure / confounder | clogit | cond. OR | OR | partial-lik info | ⚠️ `matchatr_dropped_rows` (complete-case) | `test-clogit.R` |
| mcnemar, exposure coding (logical / 2-level factor / 0-1) | mcnemar | cond. OR | OR | McNemar | ✅ identical OR across codings | `test-mcnemar.R` |
| M:1 / variable-ratio / mixed | mcnemar | — | — | — | ⛔ `matchatr_not_one_to_one` (→ clogit) | `test-mcnemar.R` |
| non-binary exposure | mcnemar | — | — | — | ⛔ `matchatr_bad_input` (→ clogit) | `test-mcnemar.R` |
| one-sided / no discordant pairs | mcnemar | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-mcnemar.R` |
| RD / RR · sandwich / bootstrap CI | mcnemar | — | — | — | ⛔ `matchatr_unidentified_estimand` / `matchatr_unsupported_variance` | `test-mcnemar.R` |
| missing pair member | mcnemar | cond. OR | OR | McNemar | ⚠️ `matchatr_dropped_rows` (complete pairs) | `test-mcnemar.R` |

The conditional OR assembly is shared with the unmatched logistic engine via
`conditional_or_result()` (exposure coefficient by term position, Wald interval
on the log scale, exponentiated). `tidy()` renders the per-term OR table (no
intercept row). The McNemar 1:1 closed form (`estimator = "mcnemar"`) lives in
its own engine (`R/mcnemar.R`, mirroring the Mantel-Haenszel closed form) and
reports through the shared `matchatr_result` contract. Effect modification
(`R/effect_modification.R`, `stratum_specific_or_result()`) builds the per-level
exposure log OR as a linear combination of the `exposure * modifier`
coefficients and reads its variance from the joint partial-likelihood vcov; M:1
and variable-ratio matching need no special handling because the conditional
likelihood treats any matched-set composition uniformly (truth-tested in
`test-clogit.R`).

## Multiple case / control groups (PHASE_4)

**Chunk 1 implemented — the unconstrained polytomous odds ratios run end to
end.** `matcha(design = unmatched_cc(), estimator = "polytomous")` fits a
baseline-category multinomial logistic via `nnet::multinom`
(`outcome ~ exposure + confounders`) for an outcome with three or more groups
(`reference =` selects the baseline group, releveled to the front). Each
non-reference equation's exposure coefficient is that subtype's log odds ratio
versus the reference; `contrast(type = "or")` reports one OR row per
(subtype, exposure-coefficient) with an information-matrix Wald interval, and
`tidy()` renders the full per-equation table with a `y.level` column. RD / RR
are rejected as unidentified without the source-population prevalences; the
robust-sandwich and bootstrap variances do not apply. The constrained
(common-OR) fit and the homogeneity LRT remain Chunk 2 (pending).

| Exposure | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|
| binary, 3-group outcome (saturated) | polytomous | subtype OR | OR | info matrix | ✅ closed-form 2×2 Woolf: OR **and** Var(log OR) per subtype (independent of multinom) | `test-polytomous.R` |
| binary + confounder | polytomous | subtype OR | OR | info matrix | ✅ truth DGP (recovers per-subtype β within 3.5 SE) + `nnet::multinom` coef/vcov equality | `test-polytomous.R` |
| continuous exposure (per unit) | polytomous | subtype OR | OR | info matrix | ✅ vs `nnet::multinom` | `test-polytomous.R` |
| factor exposure (per level) | polytomous | subtype OR per level | OR | info matrix | ✅ vs `nnet::multinom` (+ exposure reference) | `test-polytomous.R` |
| reference choice / default / character outcome | polytomous | subtype OR | OR | info matrix | ✅ baseline releveled; == explicitly-releveled `multinom` | `test-polytomous.R` |
| unused outcome level | polytomous | subtype OR | OR | info matrix | ✅ dropped before fitting | `test-polytomous.R` |
| missing exposure / confounder | polytomous | subtype OR | OR | info matrix | ⚠️ `matchatr_dropped_rows` (complete-case) | `test-polytomous.R` |
| two-group outcome | polytomous | — | — | — | ⛔ `matchatr_bad_outcome` (→ binary estimators) | `test-polytomous.R` |
| numeric / logical outcome | polytomous | — | — | — | ⛔ `matchatr_bad_outcome` | `test-polytomous.R` |
| out-of-range `reference` | polytomous | — | — | — | ⛔ `matchatr_bad_input` | `test-polytomous.R` |
| `reference` on a non-polytomous estimator | — | — | — | — | ⛔ `matchatr_bad_input` | `test-polytomous.R` |
| constant exposure | polytomous | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-polytomous.R` |
| ordered-factor exposure | polytomous | — | — | — | ⛔ `matchatr_bad_input` (polynomial contrasts) | `test-polytomous.R` |
| effect modifier on polytomous | — | — | — | — | ⛔ `matchatr_bad_input` (clogit-only) | `test-polytomous.R` |
| polytomous on a matched design | — | — | — | — | ⛔ `matchatr_bad_estimator` | `test-polytomous.R` |
| polytomous | — | RD / RR | — | — | ⛔ `matchatr_unidentified_estimand` | `test-polytomous.R` |
| polytomous OR | — | OR | — | sandwich / bootstrap | ⛔ `matchatr_unsupported_variance` | `test-polytomous.R` |

The per-subtype OR assembly lives in `R/polytomous.R`
(`multinom_exposure_or()`): it locates the exposure columns by term position
(`term_assign()`, shared with the logistic / clogit engines) and reads each
log-OR's variance by the `level:predictor` name `nnet::multinom` gives its
`vcov()`, so the lookup never assumes a coefficient ordering. `matcha()`
resolves the multi-group outcome through `resolve_polytomous_outcome()` (≥3
groups, reference releveled to the baseline) and routes it via the new
`outcome_kind` field on the dispatch; `tidy.matchatr_fit` branches to
`tidy_multinom()` for the matrix-coefficient fit. `nnet` is a recommended
(zero-dependency) R package, so no skip guards are needed.

### Constrained common-OR fit + homogeneity LRT (PHASE_4 Chunk 2)

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

## Documentation (website articles)

Every implemented estimator is demonstrated in a rendered Quarto article on the
[package website](https://etverse.github.io/matchatr/) (sources in
`vignettes/`). The site style mirrors the other etverse packages (altdoc +
Quarto, `lumen` theme).

| Article | Covers |
|---|---|
| `introduction.qmd` | Design taxonomy, the two orthogonal axes, the two-step `matcha()` / `contrast()` / `tidy()` API, what-is-identified, what-works-today |
| `unmatched-cc.qmd` | `logistic` (binary / continuous / categorical / trend / GAM-adjusted) and `mh` (Mantel-Haenszel) ORs |
| `matched-cc.qmd` | `clogit` conditional OR, `mcnemar` 1:1 matched-pair OR, `effect_modifier` stratum-specific ORs |
| `multiple-groups.qmd` | `polytomous` per-subtype ORs vs reference, the `y.level` tidy table, collinearity guard |

Articles document only implemented features; the pending phases above are not
yet covered.
