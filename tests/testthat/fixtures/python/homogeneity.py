"""Independent oracle for matchatr's etiologic-homogeneity test + pooled OR.

Reads ``homogeneity_data.csv`` (a three-group outcome ``g`` with levels control /
caseA / caseB and a binary exposure ``x``), fits the baseline-category
multinomial logistic model with statsmodels, and reconstructs -- from the
subtype exposure log-ORs and their multinomial-information covariance -- the
Wald test that the exposure OR is constant across the non-reference subtypes
(H0: beta_caseA = beta_caseB) on M-1 df, together with the efficient GLS
(inverse-variance) pooled common OR that holds under homogeneity. Writes
``homogeneity_results.csv``: the chi-square, df, p-value, and the pooled OR with
its 95% Wald interval.

This mirrors matchatr's test_homogeneity(), which computes the same Wald
statistic and GLS-pooled OR on the unconstrained multinomial fit.

Run from this directory:  python3 homogeneity.py
"""

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    data = pd.read_csv(os.path.join(HERE, "homogeneity_data.csv"))
    g = pd.Categorical(data["g"], categories=["control", "caseA", "caseB"])
    exog = sm.add_constant(data[["x"]])
    fit = sm.MNLogit(g.codes, exog).fit(disp=0)

    names = list(exog.columns)  # ['const', 'x']
    n_par = len(names)
    x_pos = names.index("x")
    n_eq = fit.params.shape[1]  # non-baseline equations (caseA, caseB)
    # MNLogit ravels cov_params equation-major ([const_A, x_A, const_B, x_B]),
    # so the exposure coefficient of equation j sits at j * n_par + x_pos.
    idx = [j * n_par + x_pos for j in range(n_eq)]
    b = np.array([fit.params.iloc[x_pos, j] for j in range(n_eq)])
    cov_full = np.asarray(fit.cov_params())
    cov = cov_full[np.ix_(idx, idx)]  # covariance of the stacked subtype log-ORs

    # Wald test of homogeneity: contrast each subtype against the first
    # (H0: all subtype log-ORs equal), W = (C b)' (C V C')^-1 (C b), df = M - 1.
    n_sub = n_eq
    contrast = np.zeros((n_sub - 1, n_sub))
    for i in range(n_sub - 1):
        contrast[i, 0] = 1.0
        contrast[i, i + 1] = -1.0
    cb = contrast @ b
    chisq = float(cb @ np.linalg.solve(contrast @ cov @ contrast.T, cb))
    df = n_sub - 1
    p_value = float(stats.chi2.sf(chisq, df))

    # GLS (inverse-variance) pooled common log-OR and its variance.
    cov_inv = np.linalg.inv(cov)
    ones = np.ones(n_sub)
    denom = float(ones @ cov_inv @ ones)
    beta_pool = float(ones @ cov_inv @ b) / denom
    var_pool = 1.0 / denom
    z = stats.norm.ppf(0.975)

    out = pd.DataFrame(
        {
            "term": ["x"],
            "chisq": [chisq],
            "df": [df],
            "p_value": [p_value],
            "pooled_or": [np.exp(beta_pool)],
            "pooled_or_low": [np.exp(beta_pool - z * np.sqrt(var_pool))],
            "pooled_or_high": [np.exp(beta_pool + z * np.sqrt(var_pool))],
        }
    )
    out.to_csv(os.path.join(HERE, "homogeneity_results.csv"), index=False)


if __name__ == "__main__":
    main()
