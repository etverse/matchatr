#' Case-control weights that map a sample back to the source population
#'
#' Computes the Rose & van der Laan (2008) case-control weights that reweight a
#' case-control sample so its weighted empirical distribution mimics the source
#' cohort. With these weights a cohort estimator (g-computation, IPW, AIPW)
#' applied to the weighted sample targets the **marginal** estimand instead of
#' the conditional odds ratio a case-control fit reports unweighted.
#'
#' @details
#' Let q0 = P(Y = 1) be the marginal outcome prevalence in the source
#' population, and let the sample hold n1 cases and n0 controls (n = n1 + n0).
#' Each subject is weighted by the ratio of its population to its sample outcome
#' frequency:
#'
#'   w_i = q0 / (n1 / n)        if subject i is a case (y_i = 1)
#'   w_i = (1 - q0) / (n0 / n)  if subject i is a control (y_i = 0)
#'
#' The weights sum to n, and the weighted case fraction is exactly q0, so the
#' weighted sample reproduces the population outcome margin. Cases are
#' oversampled in a case-control study, so they receive the smaller weight and
#' controls the larger one (for a rare outcome). This is the Rose & van der Laan
#' case-control-weighted (CCW) estimating-equation weight (Rose & van der Laan
#' 2008, *Int. J. Biostat.* 4(1); 2014, *Biometrics* 70(1)).
#'
#' The marginal prevalence q0 may be **known** (from a registry / the
#' literature) or **estimated** from the full cohort. When estimated, its
#' sampling uncertainty adds a term to the influence function; the
#' `prevalence_known` argument records which case applies so the variance layer
#' can branch on it. The weight values themselves are identical either way.
#'
#' @param prevalence Single numeric in (0, 1): the marginal outcome prevalence
#'   q0 in the source population.
#' @param y Integer/numeric 0/1 vector of outcome (case) indicators, one per
#'   sampled subject. `NA` entries are carried through as `NA` weights (the
#'   downstream estimation engine drops them) and are excluded from the n1 / n0
#'   sample-fraction counts.
#' @param prevalence_known Logical: `TRUE` (default) if q0 is supplied as a known
#'   constant, `FALSE` if it was estimated from the full cohort. Stored as the
#'   `prevalence_known` attribute on the returned vector for the variance layer.
#' @returns A numeric vector of case-control weights the same length as `y`,
#'   carrying a logical `prevalence_known` attribute. `NA` outcome entries map to
#'   `NA` weights.
#' @family weights
#' @seealso [matcha()], [contrast()], [unmatched_cc()]
#' @noRd
cc_weights <- function(prevalence, y, prevalence_known = TRUE) {
  # n1 / n0 are the realised sample case / control counts (NA excluded); the
  # sample outcome fractions n1/n and n0/n are what the population fractions
  # q0 / (1 - q0) are rescaled against.
  is_case <- !is.na(y) & y == 1L
  is_ctrl <- !is.na(y) & y == 0L
  n1 <- sum(is_case)
  n0 <- sum(is_ctrl)
  n <- n1 + n0

  # A degenerate sample (all cases or all controls) has no contrast to reweight
  # and would divide by a zero sample fraction; the binary-outcome resolver
  # rejects this upstream, so this guard is defensive.
  if (n1 == 0L || n0 == 0L) {
    rlang::abort(
      c(
        "Case-control weights need both cases and controls in the sample.",
        i = paste0(
          "Found ",
          n1,
          " case(s) and ",
          n0,
          " control(s)."
        )
      ),
      class = c("matchatr_bad_outcome", "matchatr_error")
    )
  }

  w <- rep(NA_real_, length(y))
  w[is_case] <- prevalence / (n1 / n)
  w[is_ctrl] <- (1 - prevalence) / (n0 / n)
  attr(w, "prevalence_known") <- prevalence_known
  w
}
