# a two-group outcome is rejected (needs >=3 groups)

    Code
      matcha(d, outcome = "g", exposure = "x", design = unmatched_cc(), estimator = "polytomous")
    Condition
      Error in `matcha()`:
      ! Outcome `g` has 2 group(s); the polytomous estimator needs at least three.
      i A two-group outcome is a binary case-control analysis (`estimator = "logistic"`, `"mh"`, or `"clogit"`).

# an out-of-range reference is rejected

    Code
      matcha(d, outcome = "g", exposure = "x", design = unmatched_cc(), estimator = "polytomous",
      reference = "nope")
    Condition
      Error in `matcha()`:
      ! `reference` `nope` is not one of the outcome groups.
      i Observed groups: "control", "caseA", "caseB".

# reference on a non-polytomous estimator is rejected

    Code
      matcha(d, outcome = "case", exposure = "x", design = unmatched_cc(), estimator = "logistic",
      reference = "0")
    Condition
      Error in `matcha()`:
      ! `reference` is only used by the polytomous estimator.
      i Use `estimator = "polytomous"` with a multi-group outcome; got engine `glm_logistic`.

