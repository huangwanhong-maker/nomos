"""Reference values for nomos's conditional logit and willingness to pay, from statsmodels 0.14 (BSD 3-Clause) and numpy.

    python -m venv env && env/bin/pip install statsmodels
    env/bin/python make_choice_fixtures.py

Writes choice_fixtures.hpp beside this file. Three hundred situations of three alternatives drawn by a seeded generator, each with
a price, a time and whether it has a protected lane; the choices drawn from a conditional logit with β = (−0.8, −0.05, 0.6).
statsmodels' ConditionalLogit (McFadden's model, one choice per group) gives β̂ and its standard errors; numpy gives the willingness
to pay for a lane, −β_lane/β_price, and its delta-method standard error from statsmodels' covariance.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from statsmodels.discrete.conditional_models import ConditionalLogit

OUT = Path(__file__).with_name("choice_fixtures.hpp")


def lit(a) -> str:
    return "{" + ", ".join(f"{float(x):.17g}" for x in np.ravel(a)) + "}"


def main():
    rng = np.random.default_rng(20261109)
    S, J = 300, 3
    price = rng.uniform(1, 5, size=(S, J))
    time = rng.uniform(5, 40, size=(S, J))
    lane = rng.integers(0, 2, size=(S, J)).astype(float)
    beta = np.array([-0.8, -0.05, 0.6])
    V = beta[0] * price + beta[1] * time + beta[2] * lane
    U = V + rng.gumbel(size=(S, J))
    chosen = U.argmax(axis=1)
    X = np.column_stack([price.ravel(), time.ravel(), lane.ravel()])
    y = np.zeros(S * J)
    y[np.arange(S) * J + chosen] = 1
    groups = np.repeat(np.arange(S), J)
    fit = ConditionalLogit(y, X, groups=groups).fit(method="newton", maxiter=200, tol=1e-14, disp=False)
    b, cov = fit.params, fit.cov_params()
    wtp = -b[2] / b[0]
    g = np.array([b[2] / b[0] ** 2, 0, -1 / b[0]])
    wtp_se = float(np.sqrt(g @ cov @ g))
    text = ["// Reference values for nomos's conditional logit, written by make_choice_fixtures.py (statsmodels ConditionalLogit, numpy).",
            "#pragma once", "#include <cstddef>", "#include <vector>", "", "namespace fixtures {",
            f"  inline const std::size_t choice_situations = {S}, choice_alternatives = {J};",
            f"  inline const std::vector<double> choice_price = {lit(price)};",
            f"  inline const std::vector<double> choice_time = {lit(time)};",
            f"  inline const std::vector<double> choice_lane = {lit(lane)};",
            f"  inline const std::vector<double> choice_chosen = {lit(chosen)};",
            f"  inline const std::vector<double> choice_beta = {lit(b)};",
            f"  inline const std::vector<double> choice_se = {lit(fit.bse)};",
            f"  inline const double choice_loglik = {float(fit.llf):.17g};",
            f"  inline const double choice_wtp_lane = {float(wtp):.17g}, choice_wtp_lane_se = {wtp_se:.17g};",
            "}  // namespace fixtures", ""]
    OUT.write_text("\n".join(text), encoding="utf-8")
    print("wrote", OUT, b, fit.bse)


if __name__ == "__main__":
    main()
