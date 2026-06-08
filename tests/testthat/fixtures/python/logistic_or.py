"""Independent oracle for matchatr's unmatched case-control logistic odds ratio.

Reads the shared dataset ``logistic_or_data.csv`` (a case-control sample with a
binary case indicator, a binary exposure ``x``, and a continuous confounder
``age``), fits the conditional logistic-regression model ``case ~ x + age`` with
statsmodels, and writes ``logistic_or_results.csv``: per term, the log-odds-ratio
estimate, its standard error, and the odds ratio with a 95% Wald interval taken
on the log scale and exponentiated (matchatr's interval convention).

Run from this directory:  python3 logistic_or.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "logistic_or_data.csv"))
    # Maximum-likelihood logistic regression: the exposure coefficient is the
    # conditional log odds ratio for x adjusting for age (Prentice & Pyke 1979:
    # case-control sampling shifts only the intercept).
    fit = smf.logit("case ~ x + age", data=data).fit(disp=0)
    z = stats.norm.ppf(0.975)
    beta = fit.params
    se = fit.bse
    out = pd.DataFrame(
        {
            "term": beta.index,
            "estimate": beta.to_numpy(),  # log odds ratio
            "std_error": se.to_numpy(),
            "odds_ratio": np.exp(beta.to_numpy()),
            "conf_low": np.exp(beta.to_numpy() - z * se.to_numpy()),
            "conf_high": np.exp(beta.to_numpy() + z * se.to_numpy()),
        }
    )
    out.to_csv(os.path.join(HERE, "logistic_or_results.csv"), index=False)


if __name__ == "__main__":
    main()
