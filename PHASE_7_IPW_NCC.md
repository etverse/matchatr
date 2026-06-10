# Phase 7 — Inverse Probability Weighting for Nested Case-Control

> **Status: Chunks 1–3 complete. Chunk 4 remains.**
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
| KM | — | absolute risk F_x(t) | IPW Breslow | ✅ done (Chunk 3) |
| KM | ipw_cox (multi-endpoint) | HR per endpoint | robust | deferred (Chunk 4) |
| KM | ipw_aft / additive | param/excess | robust | deferred (Chunk 4, Ch19 §19.5) |
| working-model, no Phase-1 times | — | — | ⛔ `matchatr_missing_phase1` |

## Implementation plan

- `R/weights_design.R` — `ncc_weights()` spec + Samuelsen KM weights, working-model
  (GLM/GAM) and Chen local-averaging weights.
- `R/weighted_cox.R` — `fit_ipw_cox()` (delegate to `multipleNCC::wpl`; fall back to
  `survival::coxph(weights=, robust=TRUE)` for the weighted partial likelihood),
  multiple-endpoint dispatch.
- `R/absolute_risk_ncc.R` (Chunk 3) — IPW Breslow Λ̂₀ and F̂_x(t) under reuse weights,
  via the shared `assemble_absolute_risk()` core in `R/absolute_risk.R` (split out of
  the Phase-6 `cch` path, now in `R/absolute_risk_cch.R`).

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
2. ✅ Working-model (GLM/GAM) weights + Phase-1-missing rejection.
   - `compute_ncc_weights(ncc, cohort, method, selection_formula, time, entry)` in `R/weights_design.R`
   - Builds augmented selection dataset (one row per eligible subject × event time)
   - Fits logistic GLM or GAM of selection indicator; applies product formula for π_j
   - `matchatr_missing_phase1` fires when `cohort = NULL` or the time column is absent
   - Oracle: `multipleNCC::wpl(weight.method = "glm")` agrees within 2e-2 in log-HR
3. ✅ IPW absolute risk F_x(t) under reuse weights.
   - `absolute_risk()` gains an `ipw_cox` engine path alongside the existing `cch` path
   - `ipw_breslow_ncc()` (native, `R/absolute_risk_ncc.R`) computes the inverse-
     probability-weighted Breslow cumulative baseline hazard over the deduplicated
     NCC analysis sample: dΛ̂₀(t_k) = (Σ events) / (Σ_{at risk} w_j exp(β̂ᵀ x_j)),
     with the case weight 1 and each unique control's weight 1/π_j
   - F̂_x(t) = 1 − exp(−exp(β̂ᵀ x) Λ̂₀(t)) with delta-method complementary-log-log CIs,
     reusing the shared `assemble_absolute_risk()` core split out of the `cch` path
   - Oracle: full-cohort `survival::survfit(coxph)` F_x(t) (hand-rolled Breslow agrees
     to machine precision); truth DGP (exponential, analytical F_x(t)) covered by CI
4. Multiple endpoints (HR per endpoint from one reused control set) + (optional)
   additive/AFT models (Ch19 §19.5). The multi-endpoint reuse weight is a genuine
   design decision: the theoretically-correct inverse-sampling-probability weight is
   the *primary* sampling's π_j (what `sample_ncc(incl_prob = TRUE)` already attaches),
   whereas `multipleNCC::wpl` recomputes π per endpoint from that endpoint's risk sets.
   This warrants its own chunk so the reuse semantics, the cohort-case augmentation,
   and the (approximate) oracle strategy are decided deliberately.

## Deferred items

Weight calibration (Phase 12), marginal causal contrasts under NCC sampling (Phase 10),
quota-matching weights, counter-matching weighted analysis (cross-ref Phase 5).
