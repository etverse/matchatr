# matchatr (development version)

## 2026-06-02 — PHASE_1 validation hardening (critical review)

Follow-up fixes from an adversarial review of the PHASE_1 foundation, closing
input-validation gaps before the estimation layers depend on them.

- `ratio = Inf` (or `-Inf`) now raises the classed `matchatr_bad_ratio` instead
  of an unclassed base-R error. `Inf %% 1` is `NaN` and `NaN != 0` is `NA`, so
  the old guard reached `if (NA)`; a finiteness check now runs before the modulo.
- A degenerate single-class outcome (all cases, all controls, or all `NA`) is
  now rejected with `matchatr_bad_outcome` for every encoding. The contrast
  check previously fired only for numeric columns, so an all-controls logical or
  two-level factor slipped through with `n_cases = 0`; `resolve_binary_outcome()`
  now coerces to 0/1 first and applies one uniform both-classes-present check.
- `matcha()` rejects `data` with duplicated column names (`matchatr_bad_input`).
  `[[` resolves a duplicated name to its first match, so a duplicated outcome /
  exposure / strata column would otherwise be silently chosen and the
  existence check could even blame the wrong column.
- `matcha()` rejects a column assigned to two incompatible roles
  (`matchatr_bad_input`): the outcome or exposure also appearing as a confounder
  or a design column (e.g. the exposure double-entered in `confounders`, or the
  outcome used as the matched-set id). Confounders and design columns may still
  overlap, so frequency-matching on a variable and adjusting for it remains valid.

## 2026-06-02 — Design taxonomy, data model, and two-step API

First implemented layer: the sampling-design objects and the `matcha()` fit
verb that the rest of the package builds on. No estimator runs yet — this phase
is plumbing, validation, and dispatch.

- Six design constructors returning a unified `matchatr_design` S3 object:
  `unmatched_cc()`, `matched_cc()`, `nested_cc()`, `case_cohort()`,
  `two_phase()`, and `counter_matched()`. Each carries its sampling structure
  (strata, time, matching ratio, prevalence q0, subcohort / phase columns) and a
  `weight_spec` declaring the intended weighting scheme. Structural arguments are
  validated at construction (q0 strictly in (0, 1); ratio a whole number >= 1;
  strata a non-empty character vector) with classed `matchatr_*` errors.
- `matcha(data, outcome, exposure, design, confounders, estimator)` validates
  the request, resolves the orthogonal `(design, estimator)` axes to an
  estimation engine via an internal dispatch table, and returns a
  `matchatr_fit`. The case-control-weighted family (`"ccw_gformula"`,
  `"ccw_ipw"`, `"ccw_aipw"`, `"ccw_tmle"`) routes on any design but requires a
  prevalence q0; the classical estimators are design-specific (`"logistic"` /
  `"mh"`, `"clogit"`, `"cch"`, ...). The `model` slot is `NULL` until an
  estimation engine is run.
- Weights are never a data column: the fit reserves distinct
  `details$cc_weights` and `details$design_weights` slots (case-control weights
  and inclusion-probability weights have different variance consequences) plus a
  `details$variance_kind` slot for the eventual sampling-variance correction.
- Classed rejection paths, all matching on a `matchatr_error` parent:
  `matchatr_bad_estimator` (unknown / design-incompatible estimator),
  `matchatr_bad_outcome` (non-binary case indicator),
  `matchatr_missing_prevalence` (CCW without q0),
  `matchatr_bad_prevalence` / `matchatr_bad_ratio` / `matchatr_bad_strata`
  (malformed construction), `matchatr_bad_design` (missing columns / wrong
  object), and a `matchatr_uninformative_stratum` warning when a conditional
  likelihood would drop a matched set with no case or no control.
- `print` methods for `matchatr_design` and `matchatr_fit`.

## 2026-06-02 — Package scaffold and design roadmap

Initial bootstrap of matchatr: causal inference for (matched) case-control,
nested case-control, and case-cohort study designs, as part of the etverse
ecosystem.

- Package scaffold: DESCRIPTION (Imports `causatr`, `survatr`, `survival`),
  MIT license, testthat 3 edition, Air formatting, `altdoc` + Quarto website
  (matching the other etverse packages), GitHub Actions (R-CMD-check,
  test-coverage / Codecov, format-check, altdoc), Makefile, and Claude Code
  configuration (`CLAUDE.md`, `.claude/hard-rules.md`, symlinked etverse skills,
  posit-dev skills bundle).
- Full design roadmap in `PHASE_1`–`PHASE_20` at the repository root, mapping the
  *Handbook of Statistical Methods for Case-Control Studies* (Borgan et al., 2018)
  to an implementation plan: design taxonomy and two-step API (PHASE_1); classical
  estimators (unmatched / matched / multiple-group, PHASE_2–4); time-to-event
  sampling designs (nested case-control, case-cohort, IPW-NCC, PHASE_5–7); the
  causal layer (strategy + case-control weighting g-formula / IPW / AIPW / TMLE +
  design-weighted causal survival, PHASE_8–10); and efficiency / advanced /
  extension phases (PHASE_11–20).

No estimator code is implemented yet — every phase is at `Status: DESIGN`.
