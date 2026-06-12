"""
Cross-language oracle for matchatr's case-control-weighted MARGINAL contrasts
(the PHASE_9 CCW family: ccw_gformula / ccw_ipw / ccw_aipw).

matchatr maps a case-control sample to its source cohort with the Rose & van der
Laan case-control weights and then runs a cohort causal estimator on the weighted
sample (delegated to causatr). This oracle re-derives the same marginal risk
difference / risk ratio / marginal odds ratio by **M-estimation** with the
`delicatessen` sandwich, stacking the working-model scores and the
treatment-specific means into one estimating-equation system, so the point
estimate AND the sandwich standard error are cross-checked in a second language.

The case-control weights are FIXED observation weights (q0 is known here), so they
enter every estimating equation as a constant multiplier — exactly as causatr
treats them. The three estimators differ only in which working models are stacked:

  g-formula : weighted outcome logistic Q(A,W); mu_a = E_w[Q(a,W)].
  IPW       : weighted propensity g(W);         mu_a = E_w[1{A=a} Y / g_a(W)].
  AIPW      : both, combined in the doubly-robust EIF; mu_a = E_w[ 1{A=a}/g_a (Y-Q_a) + Q_a ].

Contrasts off the two means (mu1, mu0), with delta-method SEs from the sandwich
covariance of (mu1, mu0):

  difference : mu1 - mu0
  ratio      : mu1 / mu0                       (delta on log scale)
  or         : [mu1/(1-mu1)] / [mu0/(1-mu0)]   (delta on log scale)

Compared against matchatr's contrast() in test-python-oracle.R: point estimates
to ~1e-6 (same data, same canonical estimating equations), sandwich SEs to ~1e-2
(causatr's influence-function variance vs the delicatessen sandwich).

Usage:
    cd tests/testthat/fixtures/python
    python3 ccw_marginal.py

Requires: numpy, pandas, scipy, delicatessen (>= 4.0)
"""

import numpy as np
import pandas as pd
from scipy.stats import norm
from delicatessen import MEstimator

# ── load shared data ──────────────────────────────────────────────────────────
dat = pd.read_csv("ccw_marginal_data.csv")
y = dat["case"].to_numpy(dtype=float)
a = dat["x"].to_numpy(dtype=float)
w = dat["w"].to_numpy(dtype=float)
q0 = float(dat["q0"].iloc[0])
n = len(y)

# ── Rose & van der Laan case-control weights (fixed; q0 known) ────────────────
# case:    q0 / (n1 / n)            control: (1 - q0) / (n0 / n)
n1 = float(np.sum(y == 1))
n0 = float(np.sum(y == 0))
wt = np.where(y == 1, q0 / (n1 / n), (1.0 - q0) / (n0 / n))


def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


# Design matrices: outcome ~ 1 + A + W, propensity ~ 1 + W.
Xout = np.column_stack([np.ones(n), a, w])  # (n, 3)
Xout1 = np.column_stack([np.ones(n), np.ones(n), w])  # A := 1
Xout0 = np.column_stack([np.ones(n), np.zeros(n), w])  # A := 0
Zps = np.column_stack([np.ones(n), w])  # (n, 2)

z_crit = norm.ppf(0.975)


def contrasts_from_means(mu1, mu0, cov2):
    """RD / RR / mOR + delta-method SEs from the 2x2 cov of (mu1, mu0)."""
    out = {}
    # difference: g = mu1 - mu0, grad = [1, -1]
    g = np.array([1.0, -1.0])
    rd = mu1 - mu0
    rd_se = float(np.sqrt(g @ cov2 @ g))
    out["difference"] = (rd, rd_se, rd - z_crit * rd_se, rd + z_crit * rd_se)
    # ratio on log scale: log(mu1) - log(mu0), grad = [1/mu1, -1/mu0]
    gl = np.array([1.0 / mu1, -1.0 / mu0])
    log_rr = np.log(mu1) - np.log(mu0)
    rr_se = float(np.sqrt(gl @ cov2 @ gl))
    out["ratio"] = (
        np.exp(log_rr),
        np.exp(log_rr) * rr_se,  # SE on the RR scale (delta back-transform)
        np.exp(log_rr - z_crit * rr_se),
        np.exp(log_rr + z_crit * rr_se),
    )
    # marginal OR on log scale: logit(mu1) - logit(mu0),
    # grad d logit(mu)/d mu = 1 / (mu (1 - mu))
    go = np.array([1.0 / (mu1 * (1 - mu1)), -1.0 / (mu0 * (1 - mu0))])
    log_or = np.log(mu1 / (1 - mu1)) - np.log(mu0 / (1 - mu0))
    or_se = float(np.sqrt(go @ cov2 @ go))
    out["or"] = (
        np.exp(log_or),
        np.exp(log_or) * or_se,
        np.exp(log_or - z_crit * or_se),
        np.exp(log_or + z_crit * or_se),
    )
    return out


records = []


def record(estimator, res):
    for scale, (est, se, lo, hi) in res.items():
        records.append(
            {
                "estimator": estimator,
                "scale": scale,
                "estimate": est,
                "se": se,
                "ci_lower": lo,
                "ci_upper": hi,
            }
        )


# ── g-formula: theta = [b0, b1, b2, mu1, mu0] ────────────────────────────────
def psi_gformula(theta):
    b = np.asarray(theta[0:3])
    mu1, mu0 = theta[3], theta[4]
    p = expit(Xout @ b)
    score = wt * (y - p) * Xout.T  # (3, n)
    g1 = expit(Xout1 @ b)
    g0 = expit(Xout0 @ b)
    ee_mu1 = wt * (g1 - mu1)  # (n,)
    ee_mu0 = wt * (g0 - mu0)
    return np.vstack([score, ee_mu1, ee_mu0])


m = MEstimator(psi_gformula, init=[0.0, 0.0, 0.0, q0, q0])
m.estimate(solver="lm")
mu1, mu0 = m.theta[3], m.theta[4]
cov2 = m.variance[np.ix_([3, 4], [3, 4])]
record("ccw_gformula", contrasts_from_means(mu1, mu0, cov2))

# ── IPW: theta = [a0, a1, mu1, mu0] ──────────────────────────────────────────
def psi_ipw(theta):
    al = np.asarray(theta[0:2])
    mu1, mu0 = theta[2], theta[3]
    g = expit(Zps @ al)
    score = wt * (a - g) * Zps.T  # (2, n)
    ee_mu1 = wt * (a * y / g - mu1)
    ee_mu0 = wt * ((1 - a) * y / (1 - g) - mu0)
    return np.vstack([score, ee_mu1, ee_mu0])


m = MEstimator(psi_ipw, init=[0.0, 0.0, q0, q0])
m.estimate(solver="lm")
mu1, mu0 = m.theta[2], m.theta[3]
cov2 = m.variance[np.ix_([2, 3], [2, 3])]
record("ccw_ipw", contrasts_from_means(mu1, mu0, cov2))

# ── AIPW: theta = [b0, b1, b2, a0, a1, mu1, mu0] ─────────────────────────────
def psi_aipw(theta):
    b = np.asarray(theta[0:3])
    al = np.asarray(theta[3:5])
    mu1, mu0 = theta[5], theta[6]
    p = expit(Xout @ b)
    out_score = wt * (y - p) * Xout.T  # (3, n)
    g = expit(Zps @ al)
    ps_score = wt * (a - g) * Zps.T  # (2, n)
    q1 = expit(Xout1 @ b)
    q0_ = expit(Xout0 @ b)
    ee_mu1 = wt * (a / g * (y - q1) + q1 - mu1)
    ee_mu0 = wt * ((1 - a) / (1 - g) * (y - q0_) + q0_ - mu0)
    return np.vstack([out_score, ps_score, ee_mu1, ee_mu0])


m = MEstimator(psi_aipw, init=[0.0, 0.0, 0.0, 0.0, 0.0, q0, q0])
m.estimate(solver="lm")
mu1, mu0 = m.theta[5], m.theta[6]
cov2 = m.variance[np.ix_([5, 6], [5, 6])]
record("ccw_aipw", contrasts_from_means(mu1, mu0, cov2))

# ── write results ────────────────────────────────────────────────────────────
res = pd.DataFrame.from_records(records)
res.to_csv("ccw_marginal_results.csv", index=False)
print(f"Written ccw_marginal_results.csv ({len(res)} rows)")
print(res.to_string(index=False))
