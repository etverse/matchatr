# Feature Coverage Matrix

Single source of truth for what works, what's tested, and at what fidelity.
**Every PR that changes a feature MUST update this file.**

> **Status: classical odds-ratio engines landing.** On top of the PHASE_1
> design objects, `matcha()` fit verb, `(design, estimator)` dispatch, and
> input-validation / rejection paths, the unmatched case-control **logistic /
> Mantel-Haenszel** ORs (PHASE_2) and the matched case-control **conditional
> logistic** OR (PHASE_3 Chunk 1, `survival::clogit`) and the **polytomous**
> subtype ORs for multi-group outcomes (PHASE_4, `nnet::multinom`) — including
> `test_homogeneity()`, the Wald test of whether the exposure OR is constant
> across subtypes plus the efficient pooled common OR — run end to end through
> the shared `contrast()` / `tidy()` OR layer. The **nested case-control** risk-set
> conditional partial likelihood (PHASE_5 Chunk 1) reuses that same `clogit`
> engine and reports the **hazard ratio** (`type = "hr"`; OR = HR exactly under
> risk-set sampling). The remaining estimator cells stay pending until their
> phases land.

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
| continuous / smooth confounder | logistic (GAM via `model_fn`) | cond. OR | OR | model/sandwich | ✅ == `glm` (linear) + 🟡 smooth (intentionally smoke-only: no closed-form truth for a penalized spline OR) | `test-unconditional.R` |
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

**Chunks 1 & 2 implemented — the polytomous layer is complete.**
`matcha(design = unmatched_cc(), estimator = "polytomous")` fits a
baseline-category multinomial logistic via `nnet::multinom`
(`outcome ~ exposure + confounders`) for an outcome with three or more groups
(`reference =` selects the baseline group, releveled to the front). Each
non-reference equation's exposure coefficient is that subtype's log odds ratio
versus the reference; `contrast(type = "or")` reports one OR row per
(subtype, exposure-coefficient) with an information-matrix Wald interval, and
`tidy()` renders the full per-equation table with a `y.level` column. RD / RR
are rejected as unidentified without the source-population prevalences; the
robust-sandwich and bootstrap variances do not apply. `test_homogeneity()`
(Chunk 2, below) tests whether the exposure OR is constant across subtypes and
reports the efficient pooled common OR.

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

### Common-OR pooling + homogeneity test (PHASE_4 Chunk 2)

**Implemented — the polytomous layer is complete.** `test_homogeneity(fit)`
takes an unconstrained polytomous fit and, for each exposure term, tests whether
the exposure odds ratio is constant across the disease subtypes
(H0: beta_1 = ... = beta_M) and reports the efficient pooled ("common") odds
ratio that holds under homogeneity. The test is the canonical **Wald** test of
etiologic heterogeneity (Begg & Gray 1984; `riskclustr::eh_test_subtype`):
W = (C b)' (C V C')^-1 (C b) on M − 1 df, where `b` / `V` are the stacked subtype
log-ORs and their multinomial-information covariance (reused from
`multinom_exposure_or()`). The common OR is the minimum-variance (GLS /
inverse-variance) restricted estimator, asymptotically equal to the constrained
MLE — so no constrained refit is needed and continuous confounders are handled
directly. A binary or continuous exposure gives one test; an unordered factor
exposure one test per level. `print()` / `tidy()` render the per-term common OR +
homogeneity statistic.

| Exposure | Estimand | Statistic | Variance | Status | Test |
|---|---|---|---|---|---|
| binary, saturated 3-group | pooled OR + homogeneity | Wald χ² (df = M−1) | info matrix | ✅ closed-form 2×2 Woolf: χ² **and** pooled OR (independent of multinom vcov) | `test-homogeneity.R` |
| binary + continuous confounder | pooled OR + homogeneity | Wald χ² | info matrix | ✅ GLS / `C V C'` functional of `contrast()` (exact) + size/power DGP | `test-homogeneity.R` |
| any | homogeneity p-value | Wald χ² | info matrix | ✅ vs `riskclustr::eh_test_subtype` (mlogit engine) | `test-homogeneity.R` |
| factor exposure (per level) | pooled OR + homogeneity | Wald χ² per level | info matrix | ✅ vs hand-built `nnet::multinom` contrast | `test-homogeneity.R` |
| binary, true common OR | pooled OR | — | info matrix | ✅ efficiency: pooled SE < each subtype SE (Begg & Gray) | `test-homogeneity.R` |
| non-polytomous fit (logistic / mh / clogit) | — | — | — | ⛔ `matchatr_bad_input` | `test-homogeneity.R` |
| non-estimated fit (no model) | — | — | — | ⛔ `matchatr_not_estimated` | `test-homogeneity.R` |
| bad `conf_level` / non-`matchatr_fit` | — | — | — | ⛔ `matchatr_bad_input` | `test-homogeneity.R` |

`test_homogeneity()` and `print` / `tidy.matchatr_homogeneity` live in
`R/homogeneity.R`; the per-column Wald + GLS kernel is `homogeneity_one_term()`,
reusing `multinom_exposure_or()`'s stacked log-ORs + cross-subtype covariance and
the `C V C'` contrast pattern shared with `R/effect_modification.R`.

## Nested case-control (PHASE_5)

**Chunk 1 implemented — the risk-set conditional partial likelihood (hazard
ratio).** `matcha(design = nested_cc(strata = ..., time = ...), estimator =
"clogit")` fits the same conditional partial likelihood as the matched design via
`survival::clogit` (`outcome ~ exposure + confounders + strata(set)`, each
sampled risk set a stratum), and `contrast()` reports the exposure's **hazard
ratio** (`type = "hr"`) with a partial-likelihood-information Wald interval —
OR = HR exactly under risk-set (incidence-density) sampling, with no rare-disease
caveat (Prentice & Breslow 1978). The design's `time` column records how controls
were sampled; the risk-set membership is read from `strata`, so the conditional
likelihood does not enter `time` (it feeds the later inclusion-weight / weighted-
Cox phases). Each conditional design identifies exactly one scale: requesting an
odds ratio from a risk-set design (or a hazard ratio from a matched design) is
`matchatr_unidentified_estimand`. **PHASE_5 Chunk 2** adds the exported
`sample_ncc()`, which generates an analysis-ready NCC dataset from a cohort by
risk-set (incidence-density) control sampling (with optional population-stratum
matching and delayed entry); counter-matching is the remaining chunk.

| Sampling | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|
| simple m:1 risk-set | clogit | cond. HR | HR | partial-lik info | ✅ truth DGP (CMLE recovers cohort Cox β within 3.5 SE) + `survival::clogit` pass-through | `test-nested_cc.R` |
| m:1 risk-set + confounder | clogit | cond. HR | HR | partial-lik info | ✅ truth DGP (recovers adjusted β) + `survival::clogit` pass-through | `test-nested_cc.R` |
| risk-set vs full cohort (OR = HR) | clogit | cond. HR | HR | partial-lik info | ✅ NCC β agrees with full-cohort `survival::coxph` β within combined SE | `test-nested_cc.R` |
| efficiency at β = 0 | clogit | cond. HR | HR | partial-lik info | ✅ Var ratio ≈ (m+1)/m (Goldstein & Langholz 1992), Monte-Carlo pin | `test-nested_cc.R` |
| factor exposure (per-level) | clogit | cond. HR per level | HR | partial-lik info | ✅ vs `survival::clogit` (+ reference) | `test-nested_cc.R` |
| effect modifier (`x:m`) | clogit | stratum-specific HR | HR | partial-lik info | ✅ hand-built `survival::clogit` linear combos (per-level point) | `test-nested_cc.R` |
| nested_cc, non-binary outcome | — | — | — | — | ⛔ `matchatr_bad_outcome` (continuous / 3-level) | `test-nested_cc.R` |
| odds ratio from a risk-set design | clogit | — | — | — | ⛔ `matchatr_unidentified_estimand` | `test-nested_cc.R` |
| hazard ratio from a matched design | clogit | — | — | — | ⛔ `matchatr_unidentified_estimand` | `test-nested_cc.R` |
| nested_cc | — | RD / RR | — | — | ⛔ `matchatr_unidentified_estimand` | `test-nested_cc.R` |
| nested_cc HR | — | HR | — | sandwich / bootstrap | ⛔ `matchatr_unsupported_variance` | `test-nested_cc.R` |
| risk set with no control | clogit | — | — | — | ⚠️ `matchatr_uninformative_stratum` (clogit drops it) | `test-nested_cc.R` |
| exposure constant within risk set | clogit | — | — | — | ⛔ `matchatr_unestimable_exposure` | `test-nested_cc.R` |

The NCC analysis shares the `clogit` engine (`R/clogit.R`, `fit_clogit()`) and the
`conditional_or_result()` / `stratum_specific_or_result()` assemblies with the
matched design; the only differences are the design→scale mapping (`"hr"` vs
`"or"`, resolved in `default_contrast_type()` and `contrast_clogit()`) and the
estimand label. The truth oracle is a cohort with a known Cox log-HR
(`make_ncc_cohort()`) from which an incidence-density NCC sample is drawn
(`sample_ncc_riskset()`); the full-cohort `survival::coxph` β is the
design-faithful HR the subsample recovers. Deferred: `multipleNCC` IPW
cross-checks (Phase 7).

### Risk-set control sampling — `sample_ncc()` (PHASE_5 Chunk 2)

**Implemented.** `sample_ncc(cohort, time, event, m, match, entry)` generates a
nested case-control dataset from a cohort by risk-set (incidence-density)
sampling: each event anchors a matched set holding the case and `m` controls
drawn without replacement from the subjects at risk at that failure time. The
output is an analysis-ready `data.table` (the cohort columns plus `set`, the
per-set `case` indicator, and `risk_time`) that feeds straight into
`matcha(design = nested_cc(strata = "set", time = "risk_time"))`. The sampler is
native base-R/data.table (always available, deterministically seedable via the
ambient RNG); `Epi::ccwc` is an external cross-check in the tests, not a runtime
dependency. The test-only `sample_ncc_riskset()` DGP fixture now delegates to it,
so the Chunk 1 analysis tests also exercise the exported sampler.

| Feature | Estimand | Status | Test |
|---|---|---|---|
| m:1 risk-set draw: one case + ≤m at-risk controls per set, no within-set reuse | NCC sample | ✅ structural invariants | `test-risk_set_sampling.R` |
| risk-set definition (the at-risk pool) | — | ✅ vs `Epi::ccwc` (every sampled control ∈ our eligible pool) | `test-risk_set_sampling.R` |
| sampled NCC analysed → log-HR | cond. HR | ✅ recovers cohort Cox β within 3.5 SE + agrees with full-cohort `coxph` | `test-risk_set_sampling.R` |
| efficiency at β = 0 | — | ✅ Var ratio ≈ (m+1)/m (Goldstein & Langholz 1992) | `test-risk_set_sampling.R` |
| additional matching (`match = ~ s`) | NCC sample | ✅ every control shares the case's stratum + analysis-ready | `test-risk_set_sampling.R` |
| delayed entry (`entry`) | NCC sample | ✅ not-yet-entered subjects excluded from the risk set | `test-risk_set_sampling.R` |
| fewer than m eligible at late times | NCC sample | ✅ smaller set, no error | `test-risk_set_sampling.R` |
| a control may serve before its own later event | NCC sample | ✅ sampled controls include future cohort cases | `test-risk_set_sampling.R` |
| case with no eligible control | — | ⛔ `matchatr_empty_risk_set` (snapshot) | `test-risk_set_sampling.R` |
| missing `time` / `event` column | — | ⛔ `matchatr_bad_design` | `test-risk_set_sampling.R` |
| non-0/1 or event-free `event` | — | ⛔ `matchatr_bad_outcome` | `test-risk_set_sampling.R` |
| non-whole / sub-1 / `NULL` `m` | — | ⛔ `matchatr_bad_ratio` | `test-risk_set_sampling.R` |
| `match` not a formula / names absent column | — | ⛔ `matchatr_bad_input` / `matchatr_bad_design` | `test-risk_set_sampling.R` |
| output-name collision (`set`/`case`/`risk_time` present) | — | ⛔ `matchatr_bad_input` (snapshot) | `test-risk_set_sampling.R` |
| non-data.frame cohort / non-numeric `time` | — | ⛔ `matchatr_bad_input` | `test-risk_set_sampling.R` |

`sample_ncc()` lives in `R/risk_set_sampling.R` (with the `@noRd`
`eligible_controls()` / `resolve_event_indicator()` / `check_sample_m()` /
`reject_empty_risk_set()` helpers); the `Epi::ccwc` oracle wrapper is
`tests/testthat/helper-ncc-oracle.R`.

### Counter-matched control sampling and weighted partial likelihood (PHASE_5 Chunk 3)

**Implemented.** `sample_ncc_counter_matched(cohort, time, event, surrogate, m,
match, entry)` generates a counter-matched NCC dataset: at each event time the
case is matched to `m` controls drawn from the *opposite* surrogate stratum.
The output appends `set`, `case`, `risk_time`, and `log_w` (the Langholz-Borgan
1995 log-sampling-weight). `matcha(design = counter_matched(strata = "set",
time = "risk_time", weights = "log_w"), estimator = "weighted_cox")` fits
`survival::coxph(outcome ~ exposure + confounders + strata(set) + offset(log_w))`
and `contrast()` reports the hazard ratio (`type = "hr"`). The weights encode
that the case represents its entire same-stratum risk set (log_w = log(n_same + 1))
and each control represents the opposite stratum divided by the controls drawn
(log_w = log(n_opp / m_take)); the unweighted clogit is biased for this design.

| Feature | Estimand | Status | Test |
|---|---|---|---|
| counter-matched draw: one case + ≤m opp-stratum controls per set | CM sample | ✅ structural invariants (one case, controls from opp. stratum) | `test-weighted_cox.R` |
| log_w formula (case = log(n_same+1), ctrl = log(n_opp/m)) | — | ✅ exact check on controlled micro-cohort with known risk-set counts | `test-weighted_cox.R` |
| truth recovery: weighted CM HR recovers cohort Cox β | cond. HR | ✅ within 3.5 reported SEs | `test-weighted_cox.R` |
| full-cohort coxph oracle agreement | cond. HR | ✅ CM β agrees with full-cohort `survival::coxph` β within combined SE | `test-weighted_cox.R` |
| fewer than m opp-stratum controls → smaller set, no error | CM sample | ✅ structural assertion | `test-weighted_cox.R` |
| a control may be a later cohort case | CM sample | ✅ asserted | `test-weighted_cox.R` |
| `type = "or"` or `"difference"` from counter-matched | — | ⛔ `matchatr_unidentified_estimand` | `test-weighted_cox.R` |
| `ci_method = "sandwich"` / `"bootstrap"` | — | ⛔ `matchatr_unsupported_variance` | `test-weighted_cox.R` |
| no weights column in design | — | ⛔ `matchatr_bad_design` (at fit time) | `test-weighted_cox.R` |
| continuous / multi-level / NA surrogate | — | ⛔ `matchatr_bad_input` | `test-weighted_cox.R` |
| no opposite-stratum controls at risk | — | ⛔ `matchatr_empty_risk_set` (snapshot) | `test-weighted_cox.R` |
| `log_w` column clash | — | ⛔ `matchatr_bad_input` | `test-weighted_cox.R` |
| missing surrogate column | — | ⛔ `matchatr_bad_design` | `test-weighted_cox.R` |

`sample_ncc_counter_matched()` and `resolve_surrogate()` live in
`R/risk_set_sampling.R`; the engine is `R/weighted_cox.R`
(`fit_weighted_cox()` / `contrast_weighted_cox()`). The test fixture
`sample_ncc_counter_matched_fixture()` is in `tests/testthat/helper-dgp.R`.

## Case-cohort (PHASE_6)

**Chunk 1 implemented — Prentice, Self-Prentice, and Lin-Ying pseudo-likelihood
hazard ratio.** `matcha(design = case_cohort(subcohort, time, method, id),
estimator = "cch")` subsets to cases + subcohort members, builds the
`Surv(time, status) ~ exposure + confounders` formula, and delegates to
`survival::cch`. `contrast()` reports the exposure's **hazard ratio** using the
variance `survival::cch` returns for the chosen method (Self-Prentice asymptotic;
Lin-Ying robust). The naive information-matrix SE is never used (dependent score
factors due to subcohort reuse across failure times). Borgan I/II (stratified
subcohort) and absolute risk are deferred to Chunks 2–3.

| Method | Estimand | Variance | Status | Test |
|---|---|---|---|---|
| Prentice | HR | Self-Prentice asymptotic | ✅ nwtco oracle (vs `survival::cch` exact equality) + truth DGP (within 3.5 SE) + full-cohort `coxph` agreement | `test-case_cohort.R` |
| SelfPrentice | HR | Self-Prentice asymptotic | ✅ nwtco oracle (vs `survival::cch` exact equality) | `test-case_cohort.R` |
| LinYing | HR | Lin-Ying robust | ✅ nwtco oracle (vs `survival::cch` exact equality) + truth DGP (within 3.5 SE) | `test-case_cohort.R` |
| any, `type = "or"` | — | — | ⛔ `matchatr_unidentified_estimand` | `test-case_cohort.R` |
| any, `type = "difference"` / `"ratio"` | — | — | ⛔ `matchatr_unidentified_estimand` | `test-case_cohort.R` |
| any, `ci_method = "sandwich"` / `"bootstrap"` | — | — | ⛔ `matchatr_unsupported_variance` | `test-case_cohort.R` |
| invalid method string | — | — | ⛔ `matchatr_bad_input` (snapshot) | `test-case_cohort.R` |
| missing subcohort / time column | — | — | ⛔ `matchatr_bad_design` | `test-case_cohort.R` |
| I.Borgan (stratified) | HR | Borgan asymptotic | ✅ nwtco oracle (vs `survival::cch` exact equality) + truth DGP (within 3.5 SE) | `test-case_cohort.R` |
| II.Borgan (stratified) | HR | Borgan asymptotic | ✅ nwtco oracle (vs `survival::cch` exact equality) + truth DGP (within 3.5 SE) | `test-case_cohort.R` |
| Borgan I/II without `stratum` | — | — | ⛔ `matchatr_bad_design` (snapshot) | `test-case_cohort.R` |
| missing `stratum` column in data | — | — | ⛔ `matchatr_bad_design` | `test-case_cohort.R` |
| absolute risk F_x(t), Prentice / SelfPrentice / LinYing | F̂_x(t) | IPW Breslow + delta-method log-log CI | ✅ nwtco oracle (vs full-cohort survfit, tolerance 0.06) + truth DGP (exponential, CI covers truth) | `test-absolute_risk.R` |
| absolute risk F_x(t), Borgan I/II stratified | F̂_x(t) | per-stratum IPW Breslow + delta-method CI | ✅ structural + CI-ordering check | `test-absolute_risk.R` |
| `absolute_risk()` on a non-cch/ipw_cox/ipw_aft engine (e.g. clogit) | — | — | ⛔ `matchatr_not_implemented` | `test-absolute_risk.R` |
| mismatched `newdata` columns | — | — | ⛔ `matchatr_bad_input` | `test-absolute_risk.R` |

`fit_cch()` / `contrast_cch()` / `cch_exposure_coef_names()` live in
`R/case_cohort.R`. The coefficient-name lookup builds a fresh `model.matrix()`
from the user-facing formula (matching the `cch` coefficient names exactly, which
use standard R contrasts) rather than using `term_assign()`, because
`survival::cch` rewrites its internal formula in a non-standard form.

## IPW for nested case-control (PHASE_7)

**Phase complete (Chunks 1–5).**
`sample_ncc(incl_prob = TRUE)` computes Samuelsen KM inclusion probabilities and
appends `ipw_weight` (1/π_j) and `.cohort_row` to the NCC data.
`matcha(estimator = "ipw_cox")` deduplicates controls by `.cohort_row`, fits
`coxph(weights = ipw_weight, robust = TRUE)`, and reports the hazard ratio with
the Lin-Wei robust sandwich variance. Oracle: `multipleNCC::wpl(weight.method =
"KM")` agrees exactly on log-HR and SE.
**Chunk 2** adds `compute_ncc_weights(ncc, cohort, method, selection_formula, time, entry)`:
replaces the `ipw_weight` column with GLM or GAM working-model inclusion
probabilities fitted via logistic regression (`method = "glm"`) or a generalised
additive model (`method = "gam"`) of the binary control-selection indicator across
all (eligible subject, event-time) pairs. Requires the full Phase-1 cohort; omitting
`cohort` or `time` aborts with `matchatr_missing_phase1`.
**Chunk 3** extends `absolute_risk()` to the `ipw_cox` engine: a native
inverse-probability-weighted Breslow cumulative baseline hazard over the deduplicated,
Samuelsen-weighted NCC analysis sample gives `F̂_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t))`
with delta-method complementary-log-log CIs.
**Chunk 4** reuses one control set across **multiple endpoints**. Two modes feed
the same `ipw_cox` weighted Cox: (A) sampling on the union "any-failure" event
ascertains every endpoint's cases at once, so each cause-specific endpoint is
analysed directly with `matcha(outcome = "<cause>", estimator = "ipw_cox")`; (B)
`reuse_ncc_endpoint(ncc, cohort, time, event)` augments a primary-endpoint NCC
with the secondary endpoint's unsampled cohort cases (weight 1), keeping the
controls' primary inclusion weights 1/π_j. Competing-endpoint cases (ascertained
by the sampling) keep weight 1: `ncc_ipw_analysis_data()` forces weight 1 for any
subject that is a case of the analysed endpoint **or** the failing subject of some
sampled risk set.
**Chunk 5** adds two non-Cox alternative models (Ch19 §19.5) on the same
deduplicated Samuelsen-weighted sample. `estimator = "ipw_aft"` fits a weighted
Weibull accelerated failure time model (`survival::survreg`, robust sandwich) and
`contrast(type = "af")` reports the time ratio exp(β) (acceleration factor; Kang,
Lu & Liu 2017). `estimator = "ipw_aalen"` fits the weighted constant
additive-hazards model (`lin_ying_additive()`, the Lin & Ying 1994 closed form
with a robust sandwich) and `contrast(type = "excess")` reports the excess hazard
γ (additive rate difference; Borgan & Langholz 1997) — a linear-scale, possibly
negative estimand, so its Wald interval is symmetric, not exponentiated. Both
identify one scale and reject the others; `timereg::aalen` is the additive
oracle, `survival::survreg` / `multipleNCC::KMprob` the AFT oracle.
**AFT absolute risk** extends `absolute_risk()` to the `ipw_aft` engine: the
fitted weighted Weibull is a parametric survival curve, so
F̂_x(t) = 1 − exp(−exp((log t − η̂)/σ̂)) is read directly off (β̂, σ̂) with a
delta-method complementary-log-log CI over θ = (β, log σ) using the robust
survreg sandwich — no Breslow step. It shares the cloglog inversion and result
assembly with the Cox-type engines (`R/absolute_risk_aft.R` +
`cloglog_risk_ci()` / `new_matchatr_absolute_risk()` in `R/absolute_risk.R`).

| Weight | Estimator | Estimand | Variance | Status | Test |
|---|---|---|---|---|---|
| KM (Samuelsen 1997) | `ipw_cox` | HR | Lin-Wei robust sandwich | ✅ `multipleNCC::wpl` oracle (exact); truth-DGP HR recovery (3.5-SE band) | `test-ipw_ncc.R` |
| `ipw_cox` without `incl_prob = TRUE` data | — | — | — | ⛔ `matchatr_missing_ipw_weights` | `test-ipw_ncc.R` |
| `ipw_cox` with `type = "or"` or `"difference"` | — | — | — | ⛔ `matchatr_unidentified_estimand` | `test-ipw_ncc.R` |
| `ipw_cox` with `ci_method = "bootstrap"` | — | — | — | ⛔ `matchatr_unsupported_variance` | `test-ipw_ncc.R` |
| GLM working-model weights (`method = "glm"`) | `ipw_cox` | HR | Lin-Wei robust sandwich | ✅ `multipleNCC::wpl(glm)` oracle (2e-2); truth-DGP HR recovery (3.5-SE band); structural (cases = 1, controls ≥ 1, finite) | `test-ipw_ncc.R` |
| GAM working-model weights (`method = "gam"`) | `ipw_cox` | HR | Lin-Wei robust sandwich | ✅ structural (cases = 1, controls ≥ 1, finite) | `test-ipw_ncc.R` |
| `compute_ncc_weights` without `cohort` | — | — | — | ⛔ `matchatr_missing_phase1` | `test-ipw_ncc.R` |
| `compute_ncc_weights`, `time` col absent from `cohort` | — | — | — | ⛔ `matchatr_missing_phase1` | `test-ipw_ncc.R` |
| `compute_ncc_weights` without `.cohort_row` in `ncc` | — | — | — | ⛔ `matchatr_bad_input` | `test-ipw_ncc.R` |
| IPW absolute risk F_x(t) (Chunk 3) | `ipw_cox` | F̂_x(t) | IPW Breslow + delta-method log-log CI | ✅ exact vs `survival::survfit` (weighted Breslow, 1e-8) + full-cohort survfit (sampling tol) + truth DGP (exponential, CI covers truth) | `test-absolute_risk_ncc.R` |
| `absolute_risk()` on a non-cch/ipw_cox/ipw_aft engine (clogit, ipw_aalen) | — | — | — | ⛔ `matchatr_not_implemented` | `test-absolute_risk.R`, `test-absolute_risk_aft.R` |
| multiple endpoints, combined-event reuse (mode A) | `ipw_cox` | HR per endpoint | Lin-Wei robust sandwich | ✅ `multipleNCC::wpl` oracle (exact, per endpoint) + truth-DGP recovery (both causes, 3.5-SE band) | `test-multi_endpoint.R` |
| multiple endpoints, cohort-augmented reuse (mode B) | `ipw_cox` | HR (secondary) | Lin-Wei robust sandwich | ✅ independent `KMprob` + `survival::coxph` reconstruction (1e-6) + truth-DGP recovery (3.5-SE band) | `test-multi_endpoint.R` |
| competing-endpoint case reused as a control | `ipw_cox` | HR | — | ✅ ascertained → weight 1 in the analysis sample | `test-multi_endpoint.R` |
| `reuse_ncc_endpoint` structural (augments unsampled cases, weight 1; no-op on combined-event NCC) | — | — | — | ✅ structural invariants | `test-multi_endpoint.R` |
| `reuse_ncc_endpoint` without `cohort` / `time` col absent | — | — | — | ⛔ `matchatr_missing_phase1` | `test-multi_endpoint.R` |
| `reuse_ncc_endpoint` without `.cohort_row` / missing bookkeeping col / event absent | — | — | — | ⛔ `matchatr_bad_input` | `test-multi_endpoint.R` |
| `reuse_ncc_endpoint`, secondary endpoint with no cases | — | — | — | ⛔ `matchatr_bad_outcome` | `test-multi_endpoint.R` |
| accelerated failure time, Weibull (`estimator = "ipw_aft"`) | `ipw_aft` | time ratio exp(β) | survreg robust sandwich | ✅ full-cohort `survreg` recovery (3.5-SE) + independent `KMprob` + `survreg` reconstruction (1e-6, incl. complex continuous-exposure / factor-confounder set) | `test-aft_ncc.R` |
| AFT, non-Weibull baseline (`dist = "exponential"` / `"lognormal"` / `"loglogistic"`) | `ipw_aft` | time ratio exp(β) | survreg robust sandwich | ✅ `KMprob` + `survreg` reconstruction per distribution (1e-6); right `survreg` dist threaded; exponential fixes σ=1 | `test-aft_ncc.R` |
| invalid / off-estimator `dist` | — | — | — | ⛔ `matchatr_bad_input` | `test-aft_ncc.R` |
| AFT absolute risk F_x(t), S(t\|x) for all four baselines | `ipw_aft` | F̂_x(t) | survreg robust + delta-method CI (log-log for extreme-value baselines, standardised-residual otherwise) | ✅ `predict.survreg` quantile round-trip (1e-7) + numDeriv gradient reconstruction of estimate + CI via each error CDF (1e-7, incl. factor contrasts) + full-cohort `survreg` recovery (sampling tol) + truth DGP (Weibull, CI covers); weight-agnostic (KM + GLM) | `test-absolute_risk_aft.R` |
| additive hazards, constant effects (`estimator = "ipw_aalen"`) | `ipw_aalen` | excess hazard γ (rate difference) | Lin-Ying robust sandwich | ✅ truth-DGP recovery (3.5-SE, binary + 3-level factor exposure) + `timereg::aalen` oracle (point exact 1e-6 incl. full coef vector on a complex covariate set; robust SE within 5%) | `test-additive_ncc.R` |
| time-varying additive cumulative regression B_j(t) (`excess_risk()`) | `ipw_aalen` | cumulative excess hazard B_j(t) = ∫β_j(s)ds | Aalen martingale | ✅ `timereg::aalen` (no `const()`) `cum` + `var.cum` exact (1e-8, incl. 3-level factor exposure) + truth DGP (B_x(t) = β_x·t, SE band, CI covers) | `test-excess_risk.R` |
| `excess_risk()` on a non-`ipw_aalen` engine / bad `times` | — | — | — | ⛔ `matchatr_not_implemented` / `matchatr_bad_input` | `test-excess_risk.R` |
| `ipw_aft` / `ipw_aalen` off-scale `type` / `ci_method` / non-`incl_prob` / non-nested | — | — | — | ⛔ `matchatr_unidentified_estimand` / `matchatr_unsupported_variance` / `matchatr_missing_ipw_weights` / `matchatr_bad_estimator` | `test-aft_ncc.R`, `test-additive_ncc.R` |

## Case-control-weighted causal contrasts (PHASE_9)

**Chunk 1 implemented — CCW-g-formula.** `matcha(estimator = "ccw_gformula")`
on an unmatched case-control sample carrying a known prevalence `q0`
(`unmatched_cc(prevalence = q0)`) computes Rose & van der Laan case-control
weights (`cc_weights()`, `R/weights_cc.R`) that reweight the sample to the source
population, fits a weighted g-computation via `causatr::causat(estimator =
"gcomp")` (`fit_ccw()`, `R/ccw.R`), and `contrast()` reports the **marginal**
effect — risk difference (`type = "difference"`, the default), risk ratio
(`type = "ratio"`), or marginal odds ratio (`type = "or"`) — by forwarding to
`causatr::contrast()` over the treat-all / treat-none static interventions.
Point estimate and variance are delegated to causatr; matchatr owns only the
weighting layer. The remaining `ccw_*` estimators (IPW / AIPW / TMLE) and
matched / nested CC support are pending later chunks.

| Design | Exposure | Estimator | Estimand | Contrast | Variance | Status | Test |
|---|---|---|---|---|---|---|---|
| unmatched CC + q0 | binary | ccw_gformula | marginal RD | difference | sandwich (causatr) | ✅ truth DGP + exact pseudo-cohort `causatr` oracle | `test-ccw.R` |
| unmatched CC + q0 | binary | ccw_gformula | marginal RR | ratio | sandwich (causatr) | ✅ truth DGP + exact pseudo-cohort `causatr` oracle | `test-ccw.R` |
| unmatched CC + q0 | binary | ccw_gformula | marginal OR | or | sandwich (causatr) | ✅ truth DGP (≠ conditional OR) + exact pseudo-cohort oracle | `test-ccw.R` |
| Rose & van der Laan weights `cc_weights()` | — | — | — | — | — | ✅ weighted case-fraction == q0 closed form | `test-weights_cc.R` |
| any | non-binary | ccw_gformula | — | — | — | ⛔ `matchatr_bad_input` (binary ATE only) | `test-ccw.R` |
| any (no q0) | — | ccw_gformula | — | — | — | ⛔ `matchatr_missing_prevalence` | `test-ccw.R`, `test-matcha.R` |
| unmatched CC + q0 | binary | ccw_gformula | — | hr / af / excess | — | ⛔ `matchatr_unidentified_estimand` | `test-ccw.R` |
| unmatched CC + q0 | binary | ccw_gformula | RD/RR/OR | — | bootstrap | ⛔ `matchatr_unsupported_variance` | `test-ccw.R` |

CCW-IPW / CCW-AIPW (Chunk 2), CCW-TMLE (Chunk 3, the one new targeting engine),
and estimated-q0 variance + matched/nested CC + within-stratum bootstrap
(Chunk 4) stay pending.

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
| `multiple-groups.qmd` | `polytomous` per-subtype ORs vs reference, the `y.level` tidy table, `test_homogeneity()` (Wald test + pooled common OR), collinearity guard |
| `nested-cc.qmd` | `clogit` risk-set hazard ratio (`type = "hr"`), OR = HR equivalence, `survival::clogit` / full-cohort `coxph` agreement |
| `case-cohort.qmd` | Prentice / SelfPrentice / LinYing / Borgan I/II HRs, stratified subcohort, `absolute_risk()` IPW Breslow F_x(t), design rejections |
| `ipw-ncc.qmd` | `ipw_cox` IPW weighted Cox HR, `sample_ncc(incl_prob = TRUE)` Samuelsen KM weights, GLM/GAM working-model weights, `ipw_aft` time ratio (Weibull / exponential / lognormal / loglogistic baselines) + `ipw_aalen` excess hazard, `absolute_risk()` IPW Breslow (Cox) and parametric AFT F_x(t) + cumulative-incidence plot, `excess_risk()` time-varying cumulative excess hazard B_j(t), classical vs IPW comparison, rejection paths |
| `ccw-marginal.qmd` | `ccw_gformula` case-control-weighted marginal RD / RR / marginal OR, why the marginal effect differs from the conditional OR (non-collapsibility), `cc_weights()` / `unmatched_cc(prevalence = q0)`, the `causatr` g-computation delegation, rejection paths |

Articles document only implemented features; the pending phases above are not
yet covered.

## Cross-language oracle coverage (Python / statsmodels)

Every implemented classical estimator is additionally cross-checked against an
independent **Python** (`statsmodels`) fit, so a bug shared between matchatr and
its R engine cannot hide behind a same-package comparison. Each oracle reads the
SAME committed dataset both languages share, and the Python output is a committed
CSV fixture — tests never invoke Python, so CI needs no Python toolchain (each is
guarded with `skip_if(!file.exists(...))`). Fixtures and the regeneration recipe
live in `tests/testthat/fixtures/python/` (see its `README.md`); the comparisons
are in `test-python-oracle.R`. `statsmodels` anchors the classical MLEs;
`delicatessen` is reserved for the causal / sandwich estimands of later phases.

| Estimator | matchatr | statsmodels oracle | Agreement |
|---|---|---|---|
| unmatched logistic OR | `logistic` | `Logit` | ✅ 1e-4 |
| matched conditional OR | `clogit` | `ConditionalLogit` | ✅ 1e-3 (independent partial-likelihood optimisers) |
| nested case-control HR | `clogit` | `ConditionalLogit` | ✅ 1e-3 |
| Mantel–Haenszel summary OR | `mh` | `StratifiedTable` (RBG) | ✅ 1e-4 |
| polytomous subtype ORs | `polytomous` | `MNLogit` | ✅ 1e-4 |
| homogeneity Wald χ² + pooled OR | `test_homogeneity` | `MNLogit` + GLS (hand-built) | ✅ 1e-3 |
| GLM working-model ipw_weight | `compute_ncc_weights(glm)` | `statsmodels.Logit` + product formula | ✅ 1e-6 (same algorithm, double-precision rounding only) |
