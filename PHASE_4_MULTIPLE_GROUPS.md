# Phase 4 — Multiple Case or Control Groups: Polytomous Logistic

> **Status: Chunks 1 & 2 IMPLEMENTED — phase complete.**
> Book chapters: 5 (Multiple Case or Control Groups).
> Chunk 1 (`fit_polytomous()` via `nnet::multinom`, reference handling,
> ≥3-group rejection, the per-subtype OR contrast, and the `nnet::multinom` +
> closed-form 2×2 Woolf oracles) is in `R/polytomous.R` / `test-polytomous.R`.
> Chunk 2 (`test_homogeneity()` — the Wald test of homogeneity of the exposure
> OR across subtypes + the GLS-pooled common OR) is in `R/homogeneity.R` /
> `test-homogeneity.R`.

## Scope

**In:** case-control studies with more than two outcome groups (multiple disease
subtypes, or multiple control groups) analysed by polytomous (multinomial) logistic
regression; tests of homogeneity of exposure ORs across case subtypes; optional
cross-equation constraints (common OR / multiplicative structure) for efficiency.

**Out:** ordinal-outcome proportional-odds models (could be a later extension), matched
polytomous (conditional polytomous — note as deferred).

## Key design decisions

- **Polytomous logistic via `nnet::multinom`** with a common reference category; each
  non-reference equation's coefficients are the log ORs for that case group vs the
  shared reference. Choice of reference affects interpretation, not inference.
- **Pooling beats single subtype ORs** under homogeneity (Begg & Gray 1984): the
  common odds ratio is more efficient than any one subtype OR. Report it as the
  minimum-variance (GLS / inverse-variance) combination of the unconstrained
  subtype log-ORs — the restricted estimator, asymptotically equal to the
  constrained MLE — so no constrained refit is needed.
- **Homogeneity test** of the exposure OR across case groups as a **Wald** test of
  H0: β₁ = … = β_M on the stacked subtype coefficients (the canonical etiologic-
  heterogeneity test; Begg & Gray 1984; `riskclustr::eh_test_subtype`),
  asymptotically equivalent to the saturated-vs-constrained likelihood-ratio test.
- **Largely a wrapper** — no genuinely new statistical machinery; the value is the
  unified design API + the homogeneity test + the pooled-OR convenience, all
  reusing the unconstrained `nnet::multinom` fit.
- **Reject** a single binary outcome routed to `estimator = "polytomous"`
  (`matchatr_bad_outcome` — needs ≥3 groups).

## API design

```r
fit <- matcha(data, outcome = "subtype",      # factor: control, caseA, caseB
              exposure = "x", design = unmatched_cc(),
              confounders = ~ age, estimator = "polytomous",
              reference = "control")
tidy(fit, exponentiate = TRUE)          # OR per subtype vs reference
test_homogeneity(fit)                    # Wald test + pooled common OR
```

## Support matrix

| Groups | Estimator | Structure | Estimand | Variance | Status |
|---|---|---|---|---|---|
| k≥3 outcome groups | polytomous | unconstrained | subtype OR | info matrix | ✅ |
| k≥3 | polytomous | GLS pool | common OR | info matrix | ✅ |
| k≥3 | polytomous | — | homogeneity (Wald) | Wald χ² (df = M−1) | ✅ |
| 2 groups | polytomous | — | — | — | ⛔ `matchatr_bad_outcome` |

## Implementation plan

- `R/polytomous.R` — `fit_polytomous()` (wrap `nnet::multinom`) + the per-subtype
  OR contrast.
- `R/homogeneity.R` — `test_homogeneity()` (per-exposure Wald homogeneity test +
  GLS-pooled common OR via a contrast matrix on the stacked coefficients) +
  `print` / `tidy` methods.

## Variance / inference notes

Likelihood-based: information matrix from `multinom`; Wald CIs per equation,
exponentiated to ORs. The homogeneity test is the Wald statistic
W = (C b)' (C V C')⁻¹ (C b) ~ χ²₍M−1₎ for the difference contrast `C` on the
stacked subtype coefficients (eq. 5.4–5.9 in the chapter; a constraint matrix on
the stacked coefficients), and the common OR is the GLS / inverse-variance
restricted estimator `(1' V⁻¹ b)/(1' V⁻¹ 1)`.

## Oracle testing strategy

- `nnet::multinom` is the engine; assert wrapper fidelity.
- Truth-based: simulate a cohort with two case subtypes having known (equal or unequal)
  exposure ORs; confirm recovery and that the homogeneity Wald test has correct
  size/power.
- Cross-check the efficiency ordering (pooled common-OR SE < each subtype SE under a
  true common OR), and the homogeneity p-value against `riskclustr::eh_test_subtype`.

## Chunk plan

1. ✅ `fit_polytomous()` + reference handling + ≥3-group rejection + oracle.
2. ✅ Common-OR pooling + `test_homogeneity()` (Wald test).

### Chunk 2 decision: Wald / GLS on the unconstrained fit (no refit)

The constrained common-OR fit needs cross-equation constraints that
`nnet::multinom` does not expose. Three mechanisms were weighed: (a) the Poisson
/ log-linear surrogate via `stats::glm` (LRT by `anova()`, no new dependency —
but it needs *discrete* covariate patterns: one nuisance intercept per pattern,
which a continuous confounder like the DGP's `age` blows up); (b)
`VGAM::vglm(constraints = ...)` (native constrained MLE + a genuine LRT, but a
heavy new hard dependency, diverging from sibling `causatr`, which commits to
`nnet`); and (c) a **Wald / GLS contrast on the already-fitted unconstrained
multinomial**.

**Decision: (c).** It is the canonical applied method — the reference R package
for this feature, `riskclustr::eh_test_subtype` (Zabor, MSKCC), and the methods
literature test homogeneity of a risk factor's effect across subtypes with a
**Wald** test of H0: β₁ = … = β_M, df = M − 1; no applied standard uses a literal
LRT or a constrained refit here. With the stacked subtype log-ORs `b` and their
multinomial-information covariance `V` (reused from `multinom_exposure_or()`) and
a full-rank difference contrast `C`:

    W = (C b)' (C V C')⁻¹ (C b)  ~  χ²₍M−1₎.

The common OR is the minimum-variance (GLS / inverse-variance) restricted
estimator, `b_c = (1' V⁻¹ b)/(1' V⁻¹ 1)`, `Var = 1/(1' V⁻¹ 1)` — asymptotically
equal to the constrained MLE the Poisson / VGAM routes would compute. Because the
constraint is imposed on the existing fit, there is **no refit**, **no new
dependency**, and continuous confounders are handled directly; it mirrors the
`C V C'` construction already used for stratum-specific ORs in
`R/effect_modification.R`. (The Wald test is asymptotically equivalent to the
LRT; the doc's earlier "homogeneity LRT" wording is corrected to "Wald test".)

## Deferred items

Matched polytomous (conditional polytomous logistic). Proportional-odds / ordinal
outcomes. Marginal causal contrasts for subtype effects (would compose with Phase 9).

### Pluggable multinomial fitter (`multinom_fn`) — deferred, opt-in only

`nnet::multinom` stays the **sole engine and default**. The choice is settled on
three independent grounds: it is a *recommended* R package (zero install cost),
it is the handbook chapter's choice, and — decisively — it is the **same engine
the sibling `causatr` already commits to** (`causatr` Imports `nnet` and does its
categorical-treatment multinomial residualisation via `nnet::multinom` in the
internal `prepare_model_if_multinom()` / `mv_ht_closure`; those helpers are not
exported, so there is no API to reuse, but the engine choice must stay
consistent across the etverse). `causatr` pulls in **no** other multinomial
package (`VGAM` / `brglm2` / `mlogit` / `mclogit` appear nowhere), so adding one
here would diverge matchatr from its sibling.

A pluggable `multinom_fn` argument — mirroring the logistic engine's `model_fn`
(default `nnet::multinom`, `trace = FALSE`) — is the path to an *opt-in*
alternative for the **separation** edge case, where the ML fit diverges to ±∞.
The natural alternative is `brglm2::brmultinom` (Firth-type bias reduction, the
multinomial analog of the `logistf` option already in `Suggests`). It was
verified to be a near-drop-in: same `coef()` `[levels × predictors]` matrix, same
`level:predictor` `vcov()` names, compatible `model.matrix()` / `terms()` /
`$lev`. Two adapters are needed before it can be wired:

- detection must key off the matchatr engine (`fit$engine == "multinom"`), not
  `inherits(model, "multinom")` — `brmultinom` is class
  `c("brmultinom", "brglmFit", "glm")`;
- the analysis `n` must read `nrow(stats::model.matrix(model))` (fitter-agnostic)
  rather than `nrow(model$residuals)` (which `brmultinom` lacks).

This is left for a later chunk and, if pursued, should be coordinated with
`causatr` so the etverse adopts any bias-reduced multinomial option uniformly
(it overlaps the small-sample / Firth work in PHASE_15).
