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

**Near-term (active) — Chunks 2–6 (see the Chunk plan):** the nested case-control engine
(2), two more marginal scales (3: RMST ratio + quantile survival difference), categorical /
ordinal exposures (4), continuous / shift exposures via causatr's intervention DSL (5), and
an opt-in analytic sandwich SE (6).

**Deferred-on-effort — Chunk 7:** the doubly-robust `surv_aipw` (the DR version of this
phase's own estimand). The larger extensions that surfaced in the rejection review live in
**other phase docs**, not here — efficient stochastic-MTP survival in `PHASE_8` (scope),
competing-risks / CIF in `PHASE_7` (its classical foundation), weight calibration in
`PHASE_12` — see *Where the larger extensions live*.

**Settled rejections (a 2026-06-24 literature review confirmed these stay):** the
conditional/marginal **hazard ratio** (`type = "hr"`) — a g-standardized marginal HR is
time-varying and non-collapsible, properly a separate marginal-structural-Cox engine
(Hernán 2010; Aalen, Cook & Røysland 2015; Martinussen, Vansteelandt & Andersen 2020); the
**survival difference / ratio at t** — algebraically `−RD(t)`, redundant; the
**unmatched / matched-CC designs** — no time-to-event structure to standardize; and the
`times` input validation (NA / non-positive). RMTL and years-of-life-lost reduce to a sign
flip of the RMST in a single-endpoint analysis and become distinct only with a
competing-risks engine (Deferred).

## Scope

**In:** marginal causal-survival estimands (absolute risk F_x(t), risk difference, risk
ratio, RMST difference and ratio, quantile survival difference) under case-cohort and
nested case-control sampling, g-computed on the design's inclusion-weighted Cox; binary,
categorical/ordinal, and continuous/shift exposures (the last via causatr's intervention
DSL); a design-preserving bootstrap with an opt-in analytic sandwich.

**Out:** classical HR estimation (Phases 5–7), the non-survival CCW family (Phase 9),
weight calibration (Phase 12); the efficient stochastic-MTP survival estimator and the
competing-risks / CIF engine (both Deferred — future engines, not this one).

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

| Design / input | Engine | Estimand | Contrast | Variance | Status |
|---|---|---|---|---|---|
| case-cohort | `cch` | marginal F_x(t) | difference / ratio | bootstrap | ✅ Chunk 1 |
| case-cohort | `cch` | marginal RMST | rmst | bootstrap | ✅ Chunk 1 |
| nested CC | `ipw_cox` | marginal F_x(t) / RMST | difference / ratio / rmst | design-aware bootstrap | ⏳ Chunk 2 |
| case-cohort / NCC | `cch` / `ipw_cox` | marginal RMST ratio, quantile survival diff | rmst_ratio / quantile | bootstrap | ⏳ Chunk 3 |
| categorical / ordinal exposure | `cch` / `ipw_cox` | marginal F_x(t) per level | difference / ratio | bootstrap | ⏳ Chunk 4 (each-vs-reference) |
| continuous / shift exposure | `cch` / `ipw_cox` | dose-response / shift F^a(t) | difference / ratio / rmst | bootstrap | ⏳ Chunk 5 (causatr DSL; positivity warning) |
| any | `cch` / `ipw_cox` | marginal contrast | — | analytic IF / sandwich | ⏳ Chunk 6 (opt-in, approximate) |
| case-cohort / NCC | DR | marginal RD (DR) | difference | bootstrap | ⏳ Deferred (`surv_aipw`) |
| unmatched / matched CC | — | — | — | — | ⛔ `matchatr_bad_estimator` (no time-to-event) |
| off-scale `type` (or / hr / af / excess) | — | — | — | — | ⛔ `matchatr_unidentified_estimand` (settled — HR is a separate MSM-Cox estimand) |
| `ci_method = "sandwich"` (until Chunk 6) | — | — | — | — | ⛔ `matchatr_unsupported_variance` |
| missing / non-positive `times` | — | — | — | — | ⛔ `matchatr_bad_input` (settled — input validation) |

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

The **bootstrap stays the principled default** — it is the only method that, without
caveat, captures the full sampling variability (the subcohort redraw for case-cohort; the
Samuelsen weight estimation for NCC). An analytic influence-function / sandwich SE is the
**Chunk 6** opt-in, documented as an approximate fast alternative:

- **NCC (`ipw_cox`): moderate.** Reuse the robust Lin-Wei `vcov(β̂)` the engine already
  returns + the `survfit.coxph` infinitesimal-jackknife baseline SE + a standardization
  influence function (the `riskRegression::ate` mechanics; Ozenne et al. 2017, 2020). The
  Samuelsen weights are estimated, but treating them as known is the established
  conservative practice (Samuelsen 1997; Støer & Samuelsen 2016, `multipleNCC`).
- **case-cohort (`cch`): substantial.** The correct variance carries a two-phase
  design/sampling term beyond the model term; it is derived and packaged externally
  (Etievant & Gail 2024/25, `CaseCohortCoxSurvival`, *IJE* 54(2):dyaf016; Rebora et al.
  2016; Breslow & Wellner 2007). Native is real work — the alternative is to delegate the
  pure-risk SE. Until Chunk 6, an explicit `ci_method = "sandwich"` is rejected rather than
  silently substituting the bootstrap.

The standardized-risk gradient the sandwich needs (∂F/∂β = (1−F)Λ_x x, ∂F/∂Λ₀ =
(1−F)e^{βᵀx}) is already computed in `assemble_absolute_risk()`; the missing piece is the
Cov(Λ̂₀, β̂) cross-term (the `survfit.coxph` formula) and the weighted standardization
average of the per-subject gradient.

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
2. **NCC Samuelsen-weighted `surv_gcomp`** (delegates to `fit_ipw_cox()` +
   `ipw_breslow_ncc()`, marginalize over the deduplicated risk-set sample) + a
   design-aware bootstrap (resample the cohort, redraw the risk sets, recompute the
   Samuelsen weights, refit). The engine helpers already carry the `nested_cc` branches;
   this wires the dispatch row and the Samuelsen bootstrap.
3. **Additional marginal scales**: `contrast(type = "rmst_ratio")` (RMST¹/RMST⁰ from the
   two areas already integrated, log-scale delta CI; Royston & Parmar 2013) and
   `type = "quantile"` (the difference of marginal survival-time quantiles, default the
   median, read off the two standardized curves; needs an out-of-horizon NA guard —
   Heinzl & Mittlböck 2018).
4. **Categorical / ordinal exposure**: drop the binary-only restriction. A k-level
   exposure reports each-vs-reference marginal contrasts `F̂^aₘ(t) − F̂^a₀(t)` (and the
   ratio), one curve per non-reference level — the same g-formula identification applied
   level-by-level (Keil et al. 2014; Westreich et al. 2012). Positivity is checked per
   level. An ordinal exposure may additionally render the ordered (step) dose-response.
5. **Continuous / shift exposure via causatr's intervention DSL**: `contrast()` accepts a
   `causatr_intervention` (`shift` / `scale_by` / `threshold` / `dynamic` / `static`),
   applies it to the exposure column in the standardization, and reports the dose-response
   / shift `F^a(t)` contrast with the bootstrap CI. **Prerequisite (cross-package):** causatr
   exports its internal `apply_intervention()` so matchatr reuses the application logic
   rather than `:::` or re-implementing it. Ships with a **positivity / out-of-support
   warning** and a documented Cox-linearity caveat (the dose-response shape inherits the
   `β̂ᵀx` form). This is the deterministic-plug-in estimand; the **efficient**
   stochastic-MTP estimator is a scope decision recorded in `PHASE_8` (it needs the
   sibling packages' MTP machinery, not matchatr's design layer).
6. **Opt-in analytic sandwich SE** (`ci_method = "sandwich"`): NCC first (moderate — reuse
   robust `vcov(β̂)` + survfit IJ baseline + standardization IF, conservative re weight
   estimation), case-cohort second (substantial — the two-phase term, or delegate to
   `CaseCohortCoxSurvival` / `riskRegression::ate`). Documented as approximate; the
   bootstrap stays the default. See **Variance / inference notes**. Independent of Chunks
   2–5; can land any time after Chunk 1.
7. **(Deferred) Doubly-robust `surv_aipw`.** The DR version of *this phase's* estimand, so
   it lives here: matchatr's own g-computation + augmentation on the weighted Cox (no longer
   a survatr delegation), consistent if **either** the hazard model **or** the propensity
   (+ censoring) models is correct, with the augmentation IF as the variance. Builds on
   Chunks 1–2 (the singly-robust path) + Chunk 6 (the same standardization IF machinery).
   Substantial. Refs: Robins & Rotnitzky 1992; Zhang & Schaubel 2012; Bai, Tsiatis & O'Brien
   2013 (stratified sampling).

## Where the larger extensions live (not PHASE_10 chunks)

These came up in the 2026-06-24 rejection review but belong to other phase docs, not here:

- **Efficient stochastic-MTP survival** (the positivity-respecting modified-treatment-policy
  estimand for continuous exposures, with density-ratio / TMLE inference) → a scope decision
  in `PHASE_8` ("Rejected / deferred alternatives"): matchatr ships the deterministic
  plug-in (Chunk 5) by reusing causatr's intervention DSL; the efficient estimator is the
  sibling packages' (causatr / survatr) concern, with matchatr supplying the design weights.
- **Competing-risks / CIF under sampling** → `PHASE_7`, which already owns the classical
  cause-specific / multiple-endpoint NCC analysis. A causal marginal-CIF engine (CIF
  difference / ratio, and the RMTL / years-of-life-lost that follow) extends *that*
  infrastructure; it is a candidate future phase, not a contrast scale on this
  single-endpoint engine.
- **Weight calibration** → `PHASE_12`; **transportability** → a separate future direction.
  Both compose with this engine's contrasts but are built elsewhere.

## References (verified 2026-06-24 — author / title / venue confirmed; some page spans not digit-checked)

- Hernán (2010), *Epidemiology* 21(1):13–15 — the hazards of hazard ratios.
- Aalen, Cook & Røysland (2015), *Lifetime Data Analysis* 21(4):579–593 — Cox ≠ causal HR.
- Martinussen, Vansteelandt & Andersen (2020), *Lifetime Data Analysis* 26(4):833–855 —
  hazard-contrast interpretation.
- Royston & Parmar (2013), *BMC Med Res Methodol* 13:152 — RMST (difference and ratio).
- Heinzl & Mittlböck (2018), *J Eval Clin Pract* 24(4):708–712 — quantile survival difference.
- Keil et al. (2014), *Epidemiology* 25(6):889–897; Westreich et al. (2012), *Stat Med*
  31(18):2000–2009 — parametric g-formula survival under interventions.
- Muñoz & van der Laan (2012), *Biometrics* 68(2):541–549; Kennedy (2019), *JASA*
  114(526):645–656; Díaz, Williams, Hoffman & Schenck (2023), *JASA* 118(542):846–857;
  Hejazi et al. (2021), *Biometrics* 77(4):1241–1253; Young, Hernán & Robins (2014),
  *Epidemiol Methods* 3(1):1–19 — stochastic / modified-treatment-policy interventions.
- Lee, Hudgens, Cai & Cole (2016), *Statistica Sinica* 26(2):509–526 — MSM-Cox under
  case-cohort sampling.
- Ozenne et al. (2017), *R Journal* 9(2):440–460; Ozenne et al. (2020), *Biometrical
  Journal* 62(3):751–763 — `riskRegression::ate` standardized-risk influence function.
- Etievant & Gail (2024/25), `CaseCohortCoxSurvival`, *IJE* 54(2):dyaf016; Rebora et al.
  (2016), *BMC Med Res Methodol* 16:5; Breslow & Wellner (2007), *Scand J Stat*
  34(1):86–102 — two-phase / case-cohort standardized-risk variance.
- Robins & Rotnitzky (1992); Hubbard, van der Laan & Robins (2000); Zhang & Schaubel
  (2012), *Stat Med* 31(30):4255–4268; Bai, Tsiatis & O'Brien (2013), *Biometrics*
  69(4):830–839 — locally efficient / doubly-robust treatment-specific survival
  (`surv_aipw`): consistent if either the outcome (hazard) model or the propensity
  (+ censoring) models is correct.

(The MTP references above support Chunk 5's deterministic plug-in, which is this phase's;
the *efficient* MTP scope decision is recorded in `PHASE_8` and the competing-risks / CIF
extension in `PHASE_7`, where those larger items are owned.)
