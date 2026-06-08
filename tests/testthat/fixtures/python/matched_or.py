"""Independent oracle for matchatr's matched case-control conditional OR.

Reads ``matched_or_data.csv`` (a 1:3 matched case-control sample: case indicator,
binary exposure ``x``, binary covariate ``z``, matched-set id ``set``), fits the
conditional logistic (conditional maximum likelihood) model with each set as a
stratum, and writes ``matched_or_results.csv``: per term, the conditional
log-odds-ratio, its SE, and the OR with a 95% Wald interval (log scale,
exponentiated). The conditional likelihood has no intercept (it is conditioned
out with the matched-set nuisance parameters).

Run from this directory:  python3 matched_or.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "matched_or_data.csv"))
    # ConditionalLogit conditions on the per-set totals, so the exog carries NO
    # intercept (it is removed with the matched-set nuisance parameters), exactly
    # like survival::clogit.
    fit = sm.ConditionalLogit(
        data["case"], data[["x", "z"]], groups=data["set"]
    ).fit(disp=0)
    z = stats.norm.ppf(0.975)
    beta = fit.params
    se = fit.bse
    out = pd.DataFrame(
        {
            "term": beta.index,
            "estimate": beta.to_numpy(),  # conditional log odds ratio
            "std_error": se.to_numpy(),
            "odds_ratio": np.exp(beta.to_numpy()),
            "conf_low": np.exp(beta.to_numpy() - z * se.to_numpy()),
            "conf_high": np.exp(beta.to_numpy() + z * se.to_numpy()),
        }
    )
    out.to_csv(os.path.join(HERE, "matched_or_results.csv"), index=False)


if __name__ == "__main__":
    main()
