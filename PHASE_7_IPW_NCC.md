# Phase 7 — Inverse Probability Weighting for Nested Case-Control

> **Status: Chunk 1 complete (KM weights + `ipw_cox`). Chunks 2-3 remain.**
> Book chapters: 19 (IPW in NCC), with 16, 18 background.

## Scope

**In:** the IPW reformulation of nested case-control data that *breaks the matching* and
reuses controls as a case-cohort-like sample — Samuelsen design-based (KM) weights,
working-model (GLM/GAM, Chen local-averaging) weights, the weighted Cox partial
likelihood, robust variance, multiple endpoints from one control set, IPW absolute
risk, alternative time scales, and additional/counter-matching weight adjustments.

**Out:** the classical conditional-likelihood NCC (Phase 5), calibration of the IPW
weights (Phase 12), marginal causal contrasts (Phase 10).

## Key design decisions

- **The whole point is control reuse.** Classical NCC uses each control only at its
  case's failure time (Phase 5). IPW treats the union of cases+controls as a biased
  cohort sample and weights by inverse inclusion probability π_j, enabling: reuse across
  failure times, multiple endpoints with one control set, and non-PH / alternative-time-
  scale models. This is `multipleNCC`'s territory.
- **Design-based (KM) weight** (Samuelsen 1997):
  π_{j0} = 1 − ∏_{i: j∈R(t_i)} (1 − m_i/(n(t_i)−1)), with w_j = 1/π_j, w_case = 1.
  **Working-model weights** fit a logistic/GAM for P(sampled | t_j, v_j). Both are
  supported; `weight_spec` on the `nested_cc()` design selects.
- **Design weights ≠ case-control weights.** These inclusion-probability weights answer
  a different question from the Rose & van der Laan q₀ weights (Phase 9) and have a
  different variance correction. (Invariant in `hard-rules.md`.)
- **Robust variance, not naive.** Naive inverse-Hessian SEs underestimate; use the
  Lin-Wei robust sandwich (the conservative form drops the O_p(N⁻¹) covariance term).
  Delegate to `multipleNCC::wpl` which returns the correct SEs.
- **Reject** working-model weights when Phase-1 entry/event times are unavailable
  (`matchatr_missing_phase1`).

## API design

```r
fit <- matcha(ncc, outcome = "case", exposure = "dose",
              design = nested_cc(strata = "set", time = "t",
                                 weights = ncc_weights("km")),   # or "glm"/"gam"
              confounders = ~ ageRx, estimator = "ipw_cox")
contrast(fit)                                  # HR + robust CI

# multiple endpoints reuse one control set
fit2 <- matcha(ncc, outcome = "case_alav", exposure = "sbp",
               design = nested_cc(strata = "set", time = "t",
                                  weights = ncc_weights("km")),
               estimator = "ipw_cox")
```

## Support matrix

| Weight | Estimator | Estimand | Variance | Status |
|---|---|---|---|---|
| KM (design) | ipw_cox | HR | robust sandwich | needs-test |
| GLM working-model | ipw_cox | HR | robust sandwich | needs-test |
| GAM / Chen | ipw_cox | HR | robust sandwich | needs-test |
| KM | ipw_cox (multi-endpoint) | HR per endpoint | robust | needs-test |
| KM | ipw_aft / additive | param/excess | robust | needs-test (Ch19 §19.5) |
| KM | — | absolute risk F_x(t) | IPW Breslow | needs-test |
| working-model, no Phase-1 times | — | — | ⛔ `matchatr_missing_phase1` |

## Implementation plan

- `R/weights_design.R` — `ncc_weights()` spec + Samuelsen KM weights, working-model
  (GLM/GAM) and Chen local-averaging weights.
- `R/weighted_cox.R` — `fit_ipw_cox()` (delegate to `multipleNCC::wpl`; fall back to
  `survival::coxph(weights=, robust=TRUE)` for the weighted partial likelihood),
  multiple-endpoint dispatch.
- `R/absolute_risk.R` (shared with Phase 6) — IPW Breslow Λ̂₀ and F̂_x(t) under reuse
  weights.

## Variance / inference notes

Weighted Cox score U_w(β) = Σ_{i∈D}[x_i − m̂(β; t_i)]; asymptotic covariance Σ₁ + Σ₂
where Σ₂ = sampling variance (first term Σ (1−π_i)/π_i W_iW_iᵀ; conservative form drops
the cross term). Use `multipleNCC`'s robust SE. KM/Chen weights have smaller variance
estimators than alternatives (Ch19 §19.3).

## Oracle testing strategy

- `multipleNCC::wpl` is the engine + oracle (radiation/breast-cancer & CVDNOR examples,
  handbook §19.3.1, §19.4.1).
- Truth-based: simulate a cohort with known Cox β, draw an NCC sample, compare IPW-Cox
  (KM weights) β̂ to full-cohort `coxph` and to the classical conditional NCC fit
  (Phase 5) — IPW should be at least as efficient and unbiased.
- Multiple endpoints: confirm one control set yields valid HRs for two outcomes.

## Chunk plan

1. ✅ Samuelsen KM weights + `fit_ipw_cox()` + robust variance + `multipleNCC` oracle.
   - `sample_ncc(incl_prob = TRUE)` appends `ipw_weight` (1/π_j) and `.cohort_row`
   - `samuelsen_km_weights()` implements the KM product formula (O(n × K))
   - `fit_ipw_cox()` deduplicates by `.cohort_row`, fits `coxph(robust = TRUE)`
   - `contrast_ipw_cox()` reports HR with Lin-Wei robust variance
   - Oracle: `multipleNCC::wpl(weight.method = "KM")` — exact agreement on log-HR and SE
2. Working-model (GLM/GAM/Chen) weights + Phase-1-missing rejection.
3. Multiple endpoints + IPW absolute risk + (optional) additive/AFT models.

## Deferred items

Weight calibration (Phase 12), marginal causal contrasts under NCC sampling (Phase 10),
quota-matching weights, counter-matching weighted analysis (cross-ref Phase 5).
