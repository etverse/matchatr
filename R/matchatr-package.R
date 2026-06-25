#' matchatr: Causal Inference for Case-Control, Nested Case-Control, and
#' Case-Cohort Studies
#'
#' Classical and causal estimation for (matched) case-control, nested
#' case-control, and case-cohort study designs. matchatr supplies the
#' sampling-design and weighting layer — case-control weights (the Rose and
#' van der Laan g-formula / IPW / AIPW / TMLE family) and design-based
#' inclusion weights (Samuelsen, Borgan) — delegating point estimation and
#' variance to the etverse engine \pkg{causatr} (g-computation / IPW / AIPW)
#' and to \pkg{survival} (conditional logistic, weighted Cox, case-cohort).
#' Marginal causal-survival contrasts from a sampled design are g-computed on
#' the design's own weighted Cox model (the sibling \pkg{survatr} covers causal
#' survival on full-cohort person-period data).
#'
#' The package is under active development; see the `PHASE_*.md` design docs and
#' `CLAUDE.md` at the repository root for the implementation roadmap.
#'
#' @keywords internal
#' @importFrom survival clogit coxph strata Surv
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
