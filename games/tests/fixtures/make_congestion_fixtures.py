"""Reference values for nomos's routing games with congestion, from scipy 1.18 (BSD 3-Clause), networkx 3.7 (BSD 3-Clause)
and numpy.

    python -m venv env && env/bin/pip install scipy networkx
    env/bin/python make_congestion_fixtures.py

Writes congestion_fixtures.hpp beside this file. A grid of three by three crossings, each street both ways with the
Bureau of Public Roads' time t0 (1 + 0.15 (x/C)^4), its free-flow times and capacities drawn by a seeded generator, and
five demands. Every simple path of every demand is enumerated (networkx), and Beckmann's potential (for the equilibrium)
and the total time (for the system optimum) are minimised over the path flows by scipy's SLSQP with exact gradients: a
path-based solution, independent of the link-based Frank-Wolfe it checks, polished by gradient projection on the paths
(the path flows are not unique, so SLSQP alone stops near 1e-6). The script asserts that every path in use is within
1e-10 of its demand's least cost.
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import scipy
from scipy.optimize import minimize

OUT = Path(__file__).with_name("congestion_fixtures.hpp")


def lit(a) -> str:
    return "{" + ", ".join(f"{float(x):.17g}" for x in np.ravel(a)) + "}"


def main():
    rng = np.random.default_rng(20261111)
    n = 3
    arcs = []
    for r in range(n):
        for c in range(n):
            v = r * n + c
            for dr, dc in ((0, 1), (1, 0)):
                rr, cc = r + dr, c + dc
                if rr < n and cc < n:
                    w = rr * n + cc
                    for a, b in ((v, w), (w, v)):
                        arcs.append((a, b, rng.uniform(1, 3), rng.uniform(400, 1200)))
    t0 = np.array([a[2] for a in arcs])
    cap = np.array([a[3] for a in arcs])
    k = t0 * 0.15 / cap ** 4
    demands = [(0, 8, 1500.0), (8, 0, 900.0), (2, 6, 1200.0), (1, 7, 700.0), (3, 5, 1000.0)]

    G = nx.DiGraph()
    index = {}
    for i, (a, b, _, _) in enumerate(arcs):
        G.add_edge(a, b)
        index[(a, b)] = i
    paths, owner = [], []
    for w, (o, d, _) in enumerate(demands):
        for p in nx.all_simple_paths(G, o, d):
            paths.append([index[(p[q], p[q + 1])] for q in range(len(p) - 1)])
            owner.append(w)
    P = len(paths)
    Delta = np.zeros((len(arcs), P))
    for j, p in enumerate(paths):
        Delta[p, j] = 1
    owner = np.array(owner)
    T = np.array([d[2] for d in demands])

    def time(x):
        return t0 + k * x ** 4

    def solve(system):
        # In thousands of trips and of the objective's units, which SLSQP's tolerances suit.
        s = 1000.0

        def f(u):
            x = Delta @ (s * u)
            return (float(np.sum(x * time(x))) if system else float(np.sum(t0 * x + k * x ** 5 / 5))) / s

        def g(u):
            x = Delta @ (s * u)
            c = t0 + 5 * k * x ** 4 if system else time(x)
            return Delta.T @ c

        cons = [{"type": "eq", "fun": (lambda u, w=w: np.sum(u[owner == w]) - T[w] / s), "jac": (lambda u, w=w: (owner == w).astype(float))} for w in range(len(demands))]
        u0 = np.array([T[owner[j]] / s / np.sum(owner == owner[j]) for j in range(P)])
        for _ in range(20):   # restarted from its own answer until it no longer moves
            res = minimize(f, u0, jac=g, bounds=[(0, None)] * P, constraints=cons, method="SLSQP", options={"ftol": 1e-16, "maxiter": 10000})
            moved = np.max(np.abs(res.x - u0))
            u0 = np.maximum(res.x, 0)
            if moved < 1e-13:
                break
        h = u0 * s
        # Polished by gradient projection on the paths (Jayakrishnan et al. 1994, Transp. Res. Rec. 1443, 20-29): within each
        # demand, flow moves from every path to the cheapest by a Newton step on their difference, until the spread of the
        # costs of the paths in use is below 1e-12 of the least.
        def arc_cost(x):
            return t0 + 5 * k * x ** 4 if system else time(x)

        def arc_slope(x):
            return 20 * k * x ** 3 if system else 4 * k * x ** 3

        x = Delta @ h
        for _ in range(100000):
            spread = 0.0
            for w in range(len(demands)):
                js = np.flatnonzero(owner == w)
                for j in js:
                    c = Delta[:, js].T @ arc_cost(x)
                    best = js[np.argmin(c)]
                    if j == best or h[j] <= 0:
                        continue
                    cj, cb = float(Delta[:, j] @ arc_cost(x)), float(Delta[:, best] @ arc_cost(x))
                    spread = max(spread, (cj - cb) / cb)
                    differ = Delta[:, j] != Delta[:, best]
                    step = min(h[j], (cj - cb) / float(np.sum(arc_slope(x)[differ])))
                    h[j] -= step
                    h[best] += step
                    x = x + step * (Delta[:, best] - Delta[:, j])
            if spread < 1e-12:
                break
        x = Delta @ h
        c = Delta.T @ (t0 + 5 * k * x ** 4 if system else time(x))
        for w in range(len(demands)):
            used = (owner == w) & (h > 0)
            least = c[owner == w].min()
            assert np.max(np.abs(c[used] - least)) < 1e-10 * least, (system, w, np.max(np.abs(c[used] - least)))
        return x

    ue = solve(False)
    so = solve(True)
    ue_total = float(np.sum(ue * time(ue)))
    so_total = float(np.sum(so * time(so)))
    beck = float(np.sum(t0 * ue + k * ue ** 5 / 5))

    text = f"""// Generated by make_congestion_fixtures.py with scipy {scipy.__version__}, networkx {nx.__version__} and numpy {np.__version__} (2026-11-11). Do not edit.
#pragma once
#include <cstddef>
namespace congestion_fixtures {{
constexpr std::size_t nodes = {n * n}, arcs = {len(arcs)}, demands = {len(demands)};
constexpr std::size_t tail[] = {{{", ".join(str(a[0]) for a in arcs)}}};
constexpr std::size_t head[] = {{{", ".join(str(a[1]) for a in arcs)}}};
constexpr double free_time[] = {lit(t0)};
constexpr double capacity[] = {lit(cap)};
constexpr std::size_t origin[] = {{{", ".join(str(d[0]) for d in demands)}}};
constexpr std::size_t destination[] = {{{", ".join(str(d[1]) for d in demands)}}};
constexpr double trips[] = {lit(T)};
constexpr double equilibrium[] = {lit(ue)};      // arc flows at Wardrop's equilibrium
constexpr double optimum[] = {lit(so)};          // arc flows at the system optimum
constexpr double equilibrium_total = {ue_total:.17g}, optimum_total = {so_total:.17g}, beckmann = {beck:.17g};
}}  // namespace congestion_fixtures
"""
    OUT.write_text(text, encoding="utf-8", newline="\n")
    print(f"{P} paths; total time {ue_total:.10g} at equilibrium, {so_total:.10g} at the optimum (ratio {ue_total / so_total:.6f}); wrote {OUT.name}")


if __name__ == "__main__":
    main()
