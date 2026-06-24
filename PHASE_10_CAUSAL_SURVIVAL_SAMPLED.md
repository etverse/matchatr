# Phase 10 — Design-Weighted Causal Survival for Sampled Cohorts

> **Status: Chunk 1 (case-cohort) COMPLETE.** **Approach revised 2026-06-24: the
> survatr-delegation premise was disproved by a full-cohort truth oracle and
> dropped.** The marginal survival contrast is now g-computed on matchatr's own
> design-weighted Cox + IPW-Breslow `absolute_risk()` machinery (see **Why not
> survatr** below). Book chapters: 16, 17, 19 (sampling + absolute risk);
> implements Track 3 of `PHASE_8_CAUSAL_STRATEGY`.

## Why not survatr (the disproved premise)

The original plan fed the Samuelsen / Borgan inclusion weights into
`survatr::surv_fit(estimator = "gcomp", weights = )` as observation weights. A
full-cohort truth oracle (run survatr / the estimator on the full cohort, draw a
case-cohort sample, check recovery) showed this is **biased**, for two reasons:

1. survatr's `compute_survival_curve()` marginalizes the counterfactual survival
   with an **unweighted** `mean()` — the design weights reach the hazard fit but
   not the g-formula standardization (its own comment: "IPW weighting enters in
   chunk 5").
2. More fundamentally, a **constant per-subject** inclusion weight cannot express
   the **time-varying** case-cohort risk-set weighting (a non-subcohort case is
   only ascertained at its failure, so its pre-event person-time must not enter
   the risk set at weight 1). Fed into a pooled-logistic discrete-time hazard,
   constant weights attenuate the coefficients (~20–30 Monte-Carlo SE on β in the
   reference DGP).

Verified: the subcohort *alone* (a genuine random cohort sample) is unbiased; the
subcohort + cases with constant weights is biased. The Cox partial likelihood
(`survival::cch` / the Samuelsen weighted Cox) handles the risk sets correctly,
so the fix is to **g-compute on matchatr's own weighted Cox**, not survatr's
pooled logistic.

## Revised approach — g-computation on the design's weighted Cox

The marginal risk under do(X = a) is the inclusion-weighted (Horvitz-Thompson)
average of the fitted absolute risk over the cohort covariate distribution:

  F^a(t) = Σ_j w_j · F(t | X = a, W_j) / Σ_j w_j,
  F(t | x, W) = 1 − exp(−Λ̂₀(t) · exp(β̂ᵀ x)).

- **β̂** from the design's weighted Cox — `fit_cch()` (case-cohort) or
  `fit_ipw_cox()` (nested CC) — which weight the *partial likelihood* correctly.
- **Λ̂₀(t)** from the inverse-probability Breslow baseline (`ipw_breslow_cch()` /
  `ipw_breslow_ncc()`), reused from `absolute_risk()`.
- **Standardization sample**: the subcohort (a random cohort sample, weighted by
  the inverse stratum sampling fraction `N_s/n_s`) for case-cohort; the
  deduplicated Samuelsen-weighted risk-set sample for nested CC. The constant
  inclusion weights are valid *here* because the marginalization is a survey
  (Horvitz-Thompson) population mean of a fixed risk function — unbiased — whereas
  in a hazard model they misrepresented the time-varying risk sets.
- **Contrast**: RD(t) = F^1 − F^0, RR(t) = F^1/F^0, RMST difference =
  ∫₀^τ (S^1 − S^0) over the failure-time grid. Variance is a design-preserving
  bootstrap (the cohort resample for case-cohort).
- **Scales**: new `contrast(type = "rmst")`; `"difference"` / `"ratio"` are the
  risk difference / ratio at the requested `times`. The conditional OR / HR and
  the IPW NCC scales are rejected (`matchatr_unidentified_estimand`).

## Current status & roadmap (2026-06-24)

**Shipped — Chunk 1 (case-cohort):** `matcha(design = case_cohort(...), estimator =
"surv_gcomp")` → `contrast(type = "difference" | "ratio" | "rmst", times = )` with a
design-preserving cohort-resample bootstrap. Files: `R/causal_survival_sampled.R`,
`R/contrast_surv_sampled.R`, `R/variance_surv_sampled.R`. Validated by the full-cohort
g-computation truth oracle + per-subject `absolute_risk()` agreement + an independent
dense-grid RMST oracle.

**Pending — Chunk 2 (nested case-control):** wire the `nested_cc` dispatch row (the
`fit_ipw_cox()` / `ipw_breslow_ncc()` branches already exist in the engine helpers) and a
**design-aware Samuelsen bootstrap** (resample the cohort, redraw the risk sets, recompute
the KM inclusion weights, refit). Not started.

**Deferred — Chunk 3 (doubly-robust `surv_aipw`):** see the Chunk 3 references below.

**Open questions (the current rejection paths, under literature review):** whether to
support a non-binary exposure (modified-treatment-policy / dose-response standardization),
which additional marginal survival scales to add (RMTL, survival-difference, CIF, quantile
survival), and whether an analytic influence-function / sandwich SE can replace or augment
the bootstrap. The unmatched/matched-CC rejection (no time-to-event structure) and the
`times` input validation are settled.

## Scope

**In:** marginal causal-survival estimands (absolute risk F_x(t), risk difference, risk
ratio, RMST difference) under case-cohort (Chunk 1) and nested case-control (Chunk 2)
sampling, g-computed on the design's inclusion-weighted Cox.

**Out:** classical HR estimation (Phases 5–7), the non-survival CCW family (Phase 9),
weight calibration (Phase 12), competing risks under sampling.

## API

```r
# Case-cohort -> marginal risk difference / ratio at follow-up times (Chunk 1).
fit <- matcha(cohort, outcome = "event", exposure = "x",
              design = case_cohort(subcohort = "sub", time = "t"),
              confounders = ~ age, estimator = "surv_gcomp")  # cch under the hood
contrast(fit, type = "difference", times = c(2, 5))           # marginal RD(t)
contrast(fit, type = "ratio", times = c(2, 5))                # marginal RR(t)
contrast(fit, type = "rmst", times = 5)                       # marginal RMST diff up to 5
# n_boot controls the design-preserving bootstrap (default 500).
```

`matcha()` is called with the **full cohort** (as for `cch`); the contrast standardizes
over the subcohort. RMST is a `type =` of `contrast()`, not a separate verb. The exposure
is recoded to 0/1.

## Support matrix

| Design | Engine | Estimand | Contrast | Variance | Status |
|---|---|---|---|---|---|
| case-cohort | `cch` | marginal F_x(t) | difference / ratio | bootstrap | ✅ Chunk 1 |
| case-cohort | `cch` | marginal RMST | rmst | bootstrap | ✅ Chunk 1 |
| nested CC | `ipw_cox` | marginal F_x(t) / RMST | difference / ratio / rmst | design-aware bootstrap | ⏳ Chunk 2 |
| case-cohort / NCC | DR | marginal RD (DR) | difference | — | ⏳ Chunk 3 (deferred) |
| unmatched / matched CC | — | — | — | — | ⛔ `matchatr_bad_estimator` (no time-to-event) |
| non-binary exposure | — | — | — | — | ⛔ `matchatr_bad_input` (under review) |
| off-scale `type` (or/hr/af/excess) | — | — | — | — | ⛔ `matchatr_unidentified_estimand` |
| `ci_method = "sandwich"` | — | — | — | — | ⛔ `matchatr_unsupported_variance` (under review) |
| missing / non-positive `times` | — | — | — | — | ⛔ `matchatr_bad_input` |

## Implementation

- `R/causal_survival_sampled.R` — `fit_surv_gcomp()` (delegates to `fit_cch()` /
  `fit_ipw_cox()` after `surv_gcomp_recode_exposure()`), `surv_gcomp_std_sample()` +
  `subcohort_std_weights()` (the Horvitz-Thompson standardization sample and weights).
- `R/contrast_surv_sampled.R` — `contrast_surv_gcomp()`, `surv_gcomp_marginal_risk()`
  (closed-form F^a(t) from `surv_gcomp_breslow()` + `ar_lp_from_newdata()`),
  `surv_gcomp_rmst_diff()` (left-Riemann integral of the step survival).
- `R/variance_surv_sampled.R` — `surv_gcomp_boot_ci()` + `boot_percentile_ci()` (the
  design-preserving bootstrap; Chunk 2 adds the Samuelsen design-aware variant).

## Variance / inference notes

The shipped variance is the design-preserving bootstrap (cohort resample for case-cohort;
a Samuelsen risk-set redraw for NCC in Chunk 2). An analytic influence-function / sandwich
SE for the design-weighted standardized survival is an open question (under review) — it
must account for both the marginalization (coupling subjects) and the sampling variability
of the inclusion weights, so the bootstrap is the correct default for now; an explicit
`ci_method = "sandwich"` request is rejected rather than silently substituted.

## Chunk 3 references — doubly-robust treatment-specific survival

A DR estimator is consistent for the treatment-specific survival if **either** the outcome
(hazard) model **or** the propensity (+ censoring) models is correct. References (verified):

- Robins & Rotnitzky (1992) — AIPW / locally efficient estimating-equation theory.
- Hubbard, van der Laan & Robins (2000) — locally efficient survival estimation.
- **Zhang & Schaubel (2012)**, *Statistics in Medicine* 31(30): 4255–4268 — DR for
  treatment-specific survival **and RMST** (logistic treatment or Cox hazard model).
- **Bai, Tsiatis & O'Brien (2013)**, *Biometrics* 69(4) — locally efficient AIPWCC for
  treatment-specific survival **under stratified sampling** — directly relevant to the
  sampled designs here.

Since the pivot, the DR estimator would be matchatr's own (g-computation + augmentation on
the weighted Cox), **not** a survatr delegation — so Chunk 3 is deferred on effort, not
blocked on a survatr addition.

## Oracle testing strategy

- **Full-cohort g-computation is the truth**: simulate a cohort with a confounded
  binary exposure and a survival outcome, g-compute the marginal RD(t)/RR(t)/RMST on
  the full cohort (everyone in the subcohort), then draw a case-cohort sample and
  confirm the design-weighted estimate recovers it. Because a case-cohort sample
  retains *every* case, the failure-time grid is identical between the full cohort
  and the sample, so the comparison is exact (no discretization confound). Tested as
  Monte-Carlo **unbiasedness** (the mean over subcohort draws recovers the truth; the
  per-draw band uses `expect_lt(abs(mean − truth), band)` because the marginal RD is
  small-magnitude).
- **Per-subject agreement with `absolute_risk()`**: the marginalized closed form is a
  weighted average of the validated F_x(t), so it must equal `absolute_risk()`
  averaged over the standardization sample (machine precision).
- No Python (`delicatessen`) oracle covers this estimand.

## Chunk plan

1. **(COMPLETE)** Case-cohort `surv_gcomp`: `fit_surv_gcomp()` (delegates to
   `fit_cch()`), `contrast_surv_gcomp()` (marginal RD/RR/RMST via the closed-form
   `F^a(t)` over the subcohort), `surv_gcomp_boot_ci()` (cohort-resample bootstrap),
   `contrast(type = "rmst")` + `times =`. Full-cohort truth oracle + per-subject
   `absolute_risk()` agreement + the rejection paths.
2. NCC Samuelsen-weighted `surv_gcomp` (delegates to `fit_ipw_cox()` +
   `ipw_breslow_ncc()`, marginalize over the deduplicated risk-set sample) + a
   design-aware bootstrap (resample the cohort, redraw the risk sets, recompute the
   Samuelsen weights, refit).
3. `surv_aipw` (doubly-robust) — deferred: a DR treatment-specific survival estimator
   (Zhang & Schaubel 2012; Bai, Tsiatis & O'Brien 2013) is well-established but not
   yet implemented anywhere in the etverse. The singly-robust path from Chunks 1–2
   ships regardless.

## Deferred items

Weight calibration for efficiency (Phase 12), competing risks under sampling,
transportability, and the doubly-robust `surv_aipw` (Chunk 3).
