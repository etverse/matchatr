# matchatr

Causal inference for **(matched) case-control, nested case-control (NCC), and
case-cohort** study designs. Provides design-faithful classical estimators
(conditional logistic regression, risk-set-sampled / weighted Cox, Prentice /
Self-Prentice / Borgan case-cohort) and marginal causal effects via **case-control
weighting** (the Rose & van der Laan g-formula / IPW / AIPW / TMLE family) and
**design-based inclusion weighting** (Samuelsen, Borgan). Part of the
[etverse](https://github.com/etverse) ecosystem; delegates estimation to `causatr`
(g-comp / IPW / AIPW + sandwich/bootstrap variance) wherever possible. The sibling
`survatr` covers causal survival on full-cohort person-period data; matchatr's
design-weighted marginal survival is g-computed on its own weighted Cox instead.

**Status.** PHASE_1–PHASE_9 are implemented and tested; PHASE_10 Chunk 1
(case-cohort `surv_gcomp`) has landed, with Chunks 2–3 pending. PHASE_11 onward
remain `Status: DESIGN`. `FEATURE_COVERAGE_MATRIX.md` is authoritative for
per-combination status — do not restate it here.

## Guide files

- `FEATURE_COVERAGE_MATRIX.md` — **single source of truth for "what works".** Every
  PR that changes a feature MUST update this file.
- `PHASE_*.md` — per-phase design docs in the project root (see roadmap below). They
  follow the `implement-feature` Step-1b 10-point structure. Design rationale and
  rejected alternatives live here.
- `NEWS.md` — per-chunk change record; the fullest narrative of why each piece
  landed the way it did.
- `.claude/hard-rules.md` — architecture invariants review must not re-flag as bugs.

## Project structure

This is an R package: `R/` (source), `tests/testthat/` (tests, `test-foo.R` mirrors
`R/foo.R`), `man/` (generated — do not edit), `NAMESPACE` (generated — do not edit),
`vignettes/` (long-form docs). The website is built with `altdoc` + Quarto
(`altdoc/quarto_website.yml`, `lumen` theme) to match the other etverse packages.

### R/ layout

Filenames name their contents; read the directory for detail. The groupings:

- **Design + API** — `cc_design.R`, `matcha.R`, `dispatch.R`, `contrast.R`,
  `constructors.R`; validators in `checks*.R` / `resolve.R` (all errors classed
  `matchatr_*`)
- **Sampling + weights** — `risk_set_sampling.R`, `weights_design.R`,
  `multi_endpoint.R`, `weights_cc.R`
- **Classical estimators** — `unconditional.R`, `mantel_haenszel.R`, `clogit.R`,
  `mcnemar.R`, `effect_modification.R`, `polytomous.R`, `homogeneity.R`,
  `weighted_cox.R`, `ipw_cox.R`, `aft_ncc.R`, `additive_ncc.R`, `lin_ying.R`,
  `excess_risk.R`, `aalen_cumulative.R`, `case_cohort.R`, `absolute_risk*.R`;
  shared extraction in `coef_extract.R`
- **Causal layer** — `ccw_prepare.R`, `ccw.R`, `tmle_ccw.R`, `ccw_tmle_contrast.R`,
  `variance_ccw.R`; design-weighted causal survival in `causal_survival_sampled.R`,
  `contrast_surv_sampled.R`, `variance_surv_sampled.R` (see
  `PHASE_10_CAUSAL_SURVIVAL_SAMPLED.md` for why this does not delegate to `survatr`)
- **Sampling-variance corrections** — `variance_self_prentice.R`, `variance_samuelsen.R`
- **S3 + support** — `print.R`, `summary.R`, `tidy.R`, `plot.R`, `coef.R`,
  `confint.R`, `data.R`, `matchatr-package.R`, `zzz.R`

New files: split at roughly 300 lines, following the existing per-engine split.

## Two-step API

The verb mirrors the siblings (`causatr::causat()`, `survatr::surv_fit()`):

```r
# Unmatched case-control -> conditional OR
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(),
              confounders = ~ age + smoke, estimator = "logistic")
summary(fit)                      # OR table, Wald CIs
contrast(fit, type = "or")        # exposure conditional OR + CI
# contrast(fit, type = "difference")  # -> matchatr_unidentified_estimand (need q0)

# Matched CC -> conditional OR; nested CC -> risk-set HR (swap the design object)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = matched_cc(strata = "set"),
              confounders = ~ age + smoke, estimator = "clogit")

# Nested case-control -> IPW weighted HR (Samuelsen KM weights; breaks matching)
ncc <- sample_ncc(cohort, time = "t", event = "d", m = 3, incl_prob = TRUE)
fit <- matcha(ncc, outcome = "d", exposure = "x",
              design = nested_cc(strata = "set", time = "t"), estimator = "ipw_cox")
contrast(fit)                         # HR with Lin-Wei robust sandwich variance

# Marginal causal effect from a case-control sample (Rose & van der Laan)
fit <- matcha(data, outcome = "case", exposure = "x",
              design = unmatched_cc(prevalence = 0.02),   # q0
              confounders = ~ age + smoke, estimator = "ccw_gformula")

result <- contrast(fit, type = "difference", ci_method = "sandwich")  # marginal RD
```

`matcha()` returns a `matchatr_fit`; `contrast()` (reused) returns a `matchatr_result`.

## Development commands

```r
devtools::load_all()     # load for dev
devtools::test()         # run tests
devtools::check()        # R CMD check
devtools::document()     # regenerate roxygen
```

Shell: `air format .` (format all R files).

## Code style

- `pkg::fun()` for external calls (no bare `library()`)
- `rlang::abort()` / `rlang::warn()` / `rlang::inform()` (not `stop()` / `warning()` /
  `message()`); rejection errors are classed `matchatr_*`
- `data.table` internally; return `data.table` from user-facing functions
- Roxygen on every function, including `@noRd` helpers
- Generous inline comments for math, design rationale, subtle invariants
- Do NOT remove existing comments unless the related code is also removed
- Exported functions at top of files, internal helpers below

## Testing rules

- `expect_snapshot(error = TRUE)` for error conditions
- NEVER delete/mock failing tests — fix the source
- Truth-based simulation tests mandatory for new features (cohort DGP with known truth,
  then sample a CC / NCC / case-cohort from it)
- External oracle cross-checks: `survival::clogit` / `survival::cch` (classical),
  `Epi::ccwc` (risk-set sampling), `multipleNCC` (NCC IPW), `causatr`/`survatr` on
  the explicitly reweighted pseudo-cohort (CCW), R `tmle`/`tmle3` (CCW-TMLE)
- Python cross-language oracles for every implemented classical estimator: committed
  data + result CSVs under `tests/testthat/fixtures/python/`, compared in
  `test-python-oracle.R`, `skip_if(!file.exists())`-guarded so CI needs no Python.
  `statsmodels` for the classical MLEs; `delicatessen` is reserved for the causal /
  sandwich estimands of the later CCW phases.
- Update `FEATURE_COVERAGE_MATRIX.md` in the same PR as test changes

## Cost discipline

- **Targeted tests** with `devtools::test(filter = "foo")` during development; full
  `devtools::test()` only before committing.
- **Foreground** test/check commands with `timeout: 600000`; never `run_in_background`
  for `devtools::test()` / `check()`.
- **Batch R scripts** — combine diagnostics into one `Rscript -e '...'` call.

## Constraints

- Run `devtools::test()` before committing
- Do not modify `man/` or `NAMESPACE` directly
- Run `devtools::document()` after changing roxygen comments

## Scope

matchatr owns the **sampling-design + weighting layer** for case-control-type designs
and the marginal causal contrasts they support. It **delegates** point estimation and
variance to `causatr` (g-comp / IPW / AIPW), to
`survival` (clogit / coxph / cch), `nnet` (polytomous), `multipleNCC` (NCC IPW), and
`survey` (two-phase / calibration). The ONE genuinely new estimator engine is
**CCW-TMLE** (targeting step), because the etverse has no targeted-learning code.
NOT in scope: genetics designs (handbook Ch23-28), measurement-error correction
(Ch10), case-crossover (Ch7, possible future module).

## R ecosystem integration

`DESCRIPTION` is authoritative for the dependency tiers. The non-obvious relationships:

| Package | Relationship |
|---|---|
| `causatr` | Imports — the delegated causal engine (g-comp / IPW / AIPW + variance) |
| `survatr` | **Suggests**, sibling only. matchatr's design-weighted survival deliberately does **not** delegate to it |
| `multipleNCC` | Suggests, **oracle only** — NCC IPW weights are hand-rolled, fit via `survival` |
| `timereg` (`aalen`) | Suggests, **oracle only** — the additive-hazards estimator is matchatr's |
| `tmle` | Suggests, test-only oracle for CCW-TMLE |
| `Epi` (`ccwc`) | Suggests, **test oracle** for risk-set sampling — not a runtime dependency |

## Phase roadmap (handbook chapter -> phase)

Per-phase status is in `FEATURE_COVERAGE_MATRIX.md`; each `PHASE_*.md` carries its own
design. Reference: *Handbook of Statistical Methods for Case-Control Studies* (Borgan,
Breslow, Chatterjee, Gail, Scott, Wild, 2018).

`PHASE_1_DESIGN_TAXONOMY` (Ch2) · `PHASE_2_UNMATCHED_CC` (Ch3) ·
`PHASE_3_MATCHED_CC` (Ch4) · `PHASE_4_MULTIPLE_GROUPS` (Ch5) ·
`PHASE_5_NESTED_CC` (Ch16,18) · `PHASE_6_CASE_COHORT` (Ch16,17) ·
`PHASE_7_IPW_NCC` (Ch19) · `PHASE_8_CAUSAL_STRATEGY` · `PHASE_9_CCW_CONTRASTS` ·
`PHASE_10_CAUSAL_SURVIVAL_SAMPLED` · `PHASE_11_TWO_PHASE` (Ch12) ·
`PHASE_12_CALIBRATION` (Ch13) · `PHASE_13_MULTIPLE_IMPUTATION` (Ch20) ·
`PHASE_14_SEMIPARAMETRIC_MLE` (Ch21) · `PHASE_15_SMALL_SAMPLE` (Ch8) ·
`PHASE_16_POWER` (Ch9) · `PHASE_17_ALT_RISK_MODELS` (Ch11) ·
`PHASE_18_SECONDARY_ANALYSIS` (Ch14) · `PHASE_19_SCCS` (Ch22) ·
`PHASE_20_RESPONSE_SELECTIVE` (Ch15).

## Key design decisions

- **`estimator =` selects the analysis; `design =` selects the sampling structure.**
  Two orthogonal axes. The design object carries strata, time, prevalence q₀, and
  inclusion weights; the estimator decides conditional vs marginal, OR vs HR vs RD/RR.
- **Case-control weights ≠ design weights** (see `hard-rules.md`). Both are observation
  weights into the engine, never data columns.
- **Marginal causal effects reuse causatr/survatr** by passing the q₀ / inclusion
  weights as `weights`; only the sampling-variance correction is matchatr's own.
- **CCW-TMLE is the single new estimator engine** (targeting/fluctuation) — the rest is
  delegation + a weighting/design layer.
- **Conditional likelihood for matched CC and NCC** (CMLE / partial likelihood), never
  unconditional MLE on matched-set indicators.
