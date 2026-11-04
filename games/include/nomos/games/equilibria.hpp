// All Nash equilibria of a nondegenerate bimatrix game, by support enumeration.
//
// A bimatrix game is nondegenerate when no mixed strategy with support of size k has more than k
// pure best responses (von Stengel 2007, §3.2); random payoffs give one with
// probability one. In a nondegenerate game every equilibrium (x, y) has supports I and J of equal
// size k, and x, y are the unique solutions of the indifference equations
//   Σ_{j∈J} A_ij y_j = u   (i ∈ I),   Σ_{j∈J} y_j = 1,
//   Σ_{i∈I} x_i B_ij = v   (j ∈ J),   Σ_{i∈I} x_i = 1,
// each a square linear system of order k + 1. Support enumeration solves them for every pair of
// supports of equal size and keeps the solutions that are positive on their supports and satisfy
// the best response condition, (Ay)_i ≤ u for every row and (xᵀB)_j ≤ v for every column (von
// Stengel 2007, §3.2, the algorithm of support enumeration; Dickhaut and Kaplan 1991). The number
// of equilibria of a nondegenerate game is odd (Lemke and Howson 1964; Wilson 1971).
//
// Cost: Σ_k C(m, k) C(n, k) systems of order k + 1, solved by calx's Gaussian elimination;
// exponential in min(m, n), for small games (up to about 10 × 10). On a degenerate game the
// method still returns equilibria with supports of equal size and nonsingular systems, but may
// miss others (those with supports of unequal size, and continua of equilibria).

#pragma once

#include <calx/linalg/dense.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

#include "bimatrix.hpp"

namespace nomos::games {

// An equilibrium in mixed strategies and its payoffs.
struct Equilibrium {
  std::vector<double> x, y;               // the row player's and the column player's strategies
  double row_value = 0, column_value = 0;  // xᵀAy and xᵀBy
};

// Is (x, y) a Nash equilibrium, within `tolerance`? Both must be probability vectors, and neither
// player may gain more than `tolerance` by a pure deviation.
inline bool is_equilibrium(const Bimatrix& g, std::span<const double> x, std::span<const double> y, double tolerance = 1e-9) {
  if (x.size() != g.rows() || y.size() != g.cols()) return false;
  double sx = 0, sy = 0;
  for (double v : x) {
    if (v < -tolerance) return false;
    sx += v;
  }
  for (double v : y) {
    if (v < -tolerance) return false;
    sy += v;
  }
  if (std::fabs(sx - 1) > tolerance || std::fabs(sy - 1) > tolerance) return false;
  const auto Ay = calx::linalg::multiply(g.a, y);
  const auto xB = calx::linalg::multiply_transpose(g.b, x);
  const double u = row_payoff(g, x, y), v = column_payoff(g, x, y);
  for (double t : Ay)
    if (t > u + tolerance) return false;
  for (double t : xB)
    if (t > v + tolerance) return false;
  return true;
}

namespace detail {

// Every k-subset of {0, …, n − 1}, in lexicographic order, passed to f.
inline void each_subset(std::size_t n, std::size_t k, const std::function<void(const std::vector<std::size_t>&)>& f) {
  std::vector<std::size_t> c(k);
  for (std::size_t i = 0; i < k; i++) c[i] = i;
  while (true) {
    f(c);
    std::size_t i = k;
    while (i > 0 && c[i - 1] == n - k + (i - 1)) i--;
    if (i == 0) return;
    c[i - 1]++;
    for (std::size_t j = i; j < k; j++) c[j] = c[j - 1] + 1;
  }
}

// The mixture on `support` (of the player with payoffs P(own, other)) that makes the opponent's
// strategies `against` indifferent: Σ_s z_s P(s, t) = value for t in `against`, Σ z = 1.
// Returns false if the system is singular or the solution not positive.
template <class Payoff>
bool indifferent(const std::vector<std::size_t>& support, const std::vector<std::size_t>& against, Payoff P, double tolerance,
                 std::vector<double>& z, double& value) {
  const std::size_t k = support.size();
  calx::linalg::Matrix M(k + 1), rhs(k + 1);
  for (std::size_t r = 0; r < k; r++) {
    for (std::size_t c = 0; c < k; c++) M(r, c) = P(support[c], against[r]);
    M(r, k) = -1;
  }
  for (std::size_t c = 0; c < k; c++) M(k, c) = 1;
  rhs(k, 0) = 1;
  calx::linalg::Matrix X;
  try {
    X = calx::linalg::solve(M, rhs);
  } catch (const std::runtime_error&) {
    return false;
  }
  z.assign(k, 0.0);
  for (std::size_t c = 0; c < k; c++) {
    z[c] = X(c, 0);
    if (!(z[c] > tolerance)) return false;
  }
  value = X(k, 0);
  return std::isfinite(value);
}

}  // namespace detail

// All equilibria of a nondegenerate game, ordered by support size and then lexicographically by
// supports (pure equilibria first).
inline std::vector<Equilibrium> support_enumeration(const Bimatrix& g, double tolerance = 1e-9) {
  const std::size_t m = g.rows(), n = g.cols();
  std::vector<Equilibrium> out;
  // y on J makes the rows I indifferent: the column player's mixture, read through A (rows as "against").
  const auto by_A = [&](std::size_t j, std::size_t i) { return g.a(i, j); };
  const auto by_B = [&](std::size_t i, std::size_t j) { return g.b(i, j); };
  for (std::size_t k = 1; k <= std::min(m, n); k++)
    detail::each_subset(m, k, [&](const std::vector<std::size_t>& I) {
      detail::each_subset(n, k, [&](const std::vector<std::size_t>& J) {
        std::vector<double> yJ, xI;
        double u = 0, v = 0;
        if (!detail::indifferent(J, I, by_A, tolerance, yJ, u)) return;
        if (!detail::indifferent(I, J, by_B, tolerance, xI, v)) return;
        Equilibrium e;
        e.x.assign(m, 0.0);
        e.y.assign(n, 0.0);
        for (std::size_t t = 0; t < k; t++) e.x[I[t]] = xI[t], e.y[J[t]] = yJ[t];
        // The best response condition outside the supports.
        const auto Ay = calx::linalg::multiply(g.a, e.y);
        const auto xB = calx::linalg::multiply_transpose(g.b, e.x);
        for (double t : Ay)
          if (t > u + tolerance) return;
        for (double t : xB)
          if (t > v + tolerance) return;
        e.row_value = u;
        e.column_value = v;
        out.push_back(std::move(e));
      });
    });
  return out;
}

}  // namespace nomos::games
