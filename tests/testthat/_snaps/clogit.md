# clogit rejection messages read clearly

    Code
      contrast(fit, type = "or", ci_method = "sandwich")
    Condition
      Error in `contrast()`:
      ! `ci_method = "sandwich"` is not available for the conditional logistic estimator.
      i It reports the partial-likelihood information-matrix interval; use `ci_method = "model"` (the default).

