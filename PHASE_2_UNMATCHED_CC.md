# Phase 2 — Unmatched Case-Control: Logistic OR and Mantel-Haenszel

> **Status: DESIGN**
> Book chapters: 3 (Basic Concepts and Analysis).

## Scope

**In:** the classical unmatched case-control analysis — odds ratios from 2×2 and 2×K
tables, the Mantel-Haenszel stratified OR, and unconditional logistic regression for
adjusted (conditional) ORs. Categorical, ordinal (grouped-linear), and continuous
exposures; confounder adjustment; effect modification (interaction terms).

**Out:** marginal causal effects (Phase 9 — the OR here is the *conditional* OR, which
case-control data identifies directly), matched analysis (Phase 3), small-sample/exact
(Phase 15).

## Key design decisions

- **The case-control intercept is not interpretable** — it is offset by
  log(case sampling fraction / control sampling fraction). Only the slope coefficients
  (log ORs) carry over from the cohort logistic model. The fit must document this and
  `contrast()` must refuse to report a baseline risk from an unmatched CC logistic fit
  (no q₀). Marginal risks require Phase 9 (CCW).
- **OR is the native estimand; RR/RD are not identified** from unmatched CC data alone
  without external prevalence. Enforce: `contrast(type = "difference"/"ratio")` on a
  plain `estimator = "logistic"` fit → classed error pointing to `ccw_*` + `prevalence`.
- **Delegate the regression to `stats::glm(family = binomial)`**; matchatr adds the OR
  contrast layer + the MH closed form. The model_fn is pluggable (GAM via `mgcv::gam`).
- **Mantel-Haenszel** for stratified 2×2 as a closed-form alternative to logistic when
  a single stratified OR is wanted (Robins-Breslow-Greenland variance).
- **Collapsibility caution** baked into docs: do not use crude-vs-adjusted OR change to
  diagnose confounding (OR is non-collapsible even absent confounding).

## API design

```r
fit <- matcha(data, outcome = "case", exposure = "alcohol",
              design = unmatched_cc(),
              confounders = ~ age + tobacco, estimator = "logistic")
summary(fit)            # log OR table, Wald CIs
tidy(fit, exponentiate = TRUE)   # OR + 95% CI per term

# closed-form stratified OR
matcha(data, outcome = "case", exposure = "x",
       design = unmatched_cc(strata = "agegrp"), estimator = "mh")
```

## Support matrix

| Exposure | Estimator | Estimand | Contrast | Variance | Status |
|---|---|---|---|---|---|
| binary | logistic | cond. OR | OR | model/sandwich | needs-test |
| categorical k>2 | logistic | cond. OR | OR | model/sandwich | needs-test |
| ordinal (grouped-linear) | logistic | cond. OR/trend | OR | model | needs-test |
| continuous (linear/spline) | logistic (GLM/GAM) | cond. OR | OR | model | needs-test |
| binary, stratified | mh | cond. OR | OR | RBG | needs-test |
| logistic | — | RD/RR | — | — | ⛔ `matchatr_unidentified_estimand` |

## Implementation plan

- `R/unconditional.R` — `fit_logistic_cc()` (wraps `glm`), `fit_mh()` (Mantel-Haenszel
  + Robins-Breslow-Greenland variance), OR contrast assembly.
- `R/contrast.R` (shared) — OR contrast on the log scale with Wald CI; reject
  unidentified RD/RR for plain logistic.
- S3: `tidy`/`summary`/`print` for the OR table.

## Variance / inference notes

Logistic: model-based information matrix (default) or `sandwich::sandwich()` (robust).
For OR contrasts the SE is the coefficient SE on the log scale; exponentiate the CI
bounds. MH: Robins-Breslow-Greenland variance for log OR_MH. No causatr dependency
here — this is the classical layer.

## Oracle testing strategy

- `survival`/base `glm` agreement (trivial — we wrap it; assert pass-through correctness
  on coefficients/SE).
- Book values: Chapter 3 esophageal-cancer alcohol ORs; Framingham age-adjusted OR
  e^0.88 ≈ 2.42; UK oral-contraceptive ORs (Table 3.11).
- MH: cross-check against `stats::mantelhaen.test` and `epitools::epitab`.
- 2×2 closed form: hand-computed `n11 n00 / (n01 n10)` with the Woolf variance.

## Chunk plan

1. `fit_logistic_cc()` + OR contrast + unidentified-estimand rejection.
2. Categorical / ordinal / continuous exposure handling + book-value tests.
3. Mantel-Haenszel closed form + RBG variance + oracle.

## Deferred items

Marginal RD/RR/OR with prevalence (Phase 9). Exact / Firth small-sample (Phase 15).
Polytomous outcome / multiple case groups (Phase 4).
