# matchatr (development version)

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
