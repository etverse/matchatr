# Phase 4 — Multiple Case or Control Groups: Polytomous Logistic

> **Status: Chunk 1 IMPLEMENTED · Chunk 2 DESIGN**
> Book chapters: 5 (Multiple Case or Control Groups).
> Chunk 1 (`fit_polytomous()` via `nnet::multinom`, reference handling,
> ≥3-group rejection, the per-subtype OR contrast, and the `nnet::multinom` +
> closed-form 2×2 Woolf oracles) is in `R/polytomous.R` / `test-polytomous.R`.
> Chunk 2 (constrained common-OR fit + `test_homogeneity()` LRT) remains DESIGN.

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
- **Joint (constrained) fitting beats pairwise** when constraints are defensible
  (Begg & Gray 1984): fitting all groups in one model is more efficient than separate
  binary logistic fits. Offer an unconstrained fit by default and a constrained
  ("common OR across subtypes") option via a constraint specification.
- **Homogeneity test** of the exposure OR across case groups as a likelihood-ratio test
  between the saturated (group-specific OR) and constrained (common OR) models.
- **Largely a wrapper** — no genuinely new statistical machinery; the value is the
  unified design API + the homogeneity test + constrained-fit convenience.
- **Reject** a single binary outcome routed to `estimator = "polytomous"`
  (`matchatr_bad_outcome` — needs ≥3 groups).

## API design

```r
fit <- matcha(data, outcome = "subtype",      # factor: control, caseA, caseB
              exposure = "x", design = unmatched_cc(),
              confounders = ~ age, estimator = "polytomous",
              reference = "control")
tidy(fit, exponentiate = TRUE)          # OR per subtype vs reference
test_homogeneity(fit)                   # LRT: common vs subtype-specific OR
```

## Support matrix

| Groups | Estimator | Structure | Estimand | Variance | Status |
|---|---|---|---|---|---|
| k≥3 outcome groups | polytomous | unconstrained | subtype OR | info matrix | needs-test |
| k≥3 | polytomous | common-OR constraint | shared OR | info matrix | needs-test |
| k≥3 | polytomous | — | homogeneity LRT | LRT | needs-test |
| 2 groups | polytomous | — | — | — | ⛔ `matchatr_bad_outcome` |

## Implementation plan

- `R/polytomous.R` — `fit_polytomous()` (wrap `nnet::multinom`), `test_homogeneity()`
  (LRT between saturated and common-OR fits via a constraint matrix), tidy/summary.

## Variance / inference notes

Likelihood-based: information matrix from `multinom`; Wald CIs per equation,
exponentiated to ORs; LRT for homogeneity and for nested constraints. Constrained
inference (eq. 5.4–5.9 in the chapter) via a constraint matrix on the stacked
coefficients.

## Oracle testing strategy

- `nnet::multinom` is the engine; assert wrapper fidelity.
- Truth-based: simulate a cohort with two case subtypes having known (equal or unequal)
  exposure ORs; confirm recovery and that the homogeneity LRT has correct size/power.
- Cross-check pairwise-vs-joint efficiency ordering (joint SE ≤ pairwise SE under a true
  common OR).

## Chunk plan

1. ✅ `fit_polytomous()` + reference handling + ≥3-group rejection + oracle.
2. Constrained (common-OR) fit + `test_homogeneity()` LRT.

### Chunk 2 carry-forward notes

The constrained common-OR fit needs cross-equation constraints that
`nnet::multinom` does not expose. The leading candidate is the Poisson /
log-linear surrogate via `stats::glm` (one Poisson cell per outcome group ×
covariate pattern), which represents the multinomial likelihood exactly and lets
the common-OR constraint be imposed by *omitting* the exposure × group
interaction — so the homogeneity LRT is a plain `anova()` of the constrained vs
saturated Poisson glm, with **no new dependency**. `VGAM::vglm(family =
multinomial, constraints = ...)` is the alternative (native constraint matrices,
but a new hard dependency). Decide in Chunk 2.

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
