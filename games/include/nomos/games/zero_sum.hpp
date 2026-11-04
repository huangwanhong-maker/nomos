// Two-player zero-sum games: their value and optimal strategies, by linear programming.
//
// In a zero-sum game the row player receives A_ij and the column player −A_ij. The row player
// can guarantee max_x min_j (xᵀA)_j, the column player can hold the row player to
// min_y max_i (Ay)_i, and the two are equal (the minimax theorem, von Neumann 1928): their common
// value is the value of the game, and the strategies attaining them are the optimal (maximin and
// minimax) strategies, whose pairs are exactly the Nash equilibria of (A, −A). Each is a linear
// programme (Dantzig 1951; Chvátal 1983, ch. 15),
//
//   row player:     maximise v  subject to  Σ_i x_i A_ij ≥ v  (every j),  Σ x = 1,  x ≥ 0;
//   column player:  minimise w  subject to  Σ_j A_ij y_j ≤ w  (every i),  Σ y = 1,  y ≥ 0,
//
// the dual of each other, solved here separately by calx's simplex method; v* = w* is checked by
// the caller if wanted (both are returned). The optimal strategies need not be unique; the simplex
// method returns a vertex of each optimal set.

#pragma once

#include <calx/optimisation/linear_programming.hpp>

#include <limits>
#include <stdexcept>
#include <vector>

#include "bimatrix.hpp"

namespace nomos::games {

struct ZeroSumSolution {
  double lower_value = 0;  // max_x min_j (xᵀA)_j, the row player's guarantee
  double upper_value = 0;  // min_y max_i (Ay)_i, the column player's
  double value = 0;        // their mean; equal to both within the solver's tolerance
  std::vector<double> x, y;
};

namespace detail {
// max v s.t. Σ_i z_i M_ij ≥ v for every column j of M, Σ z = 1, z ≥ 0.
inline std::pair<double, std::vector<double>> maximin(const Rectangular& M) {
  const std::size_t m = M.rows, n = M.cols, nv = m + 1;
  calx::optimisation::LinearProgram lp;
  lp.c.assign(nv, 0.0);
  lp.c.back() = -1;
  lp.a_ub = Rectangular(n, nv);
  lp.b_ub.assign(n, 0.0);
  for (std::size_t j = 0; j < n; j++) {
    for (std::size_t i = 0; i < m; i++) lp.a_ub(j, i) = -M(i, j);
    lp.a_ub(j, m) = 1;
  }
  lp.a_eq = Rectangular(1, nv);
  for (std::size_t i = 0; i < m; i++) lp.a_eq(0, i) = 1;
  lp.b_eq = {1};
  lp.lower.assign(nv, 0.0);
  lp.lower.back() = -std::numeric_limits<double>::infinity();
  lp.upper.assign(nv, std::numeric_limits<double>::infinity());
  const auto r = calx::optimisation::simplex(lp);
  if (r.status != calx::optimisation::LpStatus::optimal) throw std::runtime_error("nomos::games::solve_zero_sum: the linear programme did not reach an optimum");
  return {r.x[m], std::vector<double>(r.x.begin(), r.x.begin() + static_cast<std::ptrdiff_t>(m))};
}
}  // namespace detail

// The value and optimal strategies of the zero-sum game with payoffs A to the row player.
inline ZeroSumSolution solve_zero_sum(const Rectangular& a) {
  detail::require(a.rows > 0 && a.cols > 0, "a nonempty payoff matrix");
  ZeroSumSolution s;
  auto [v, x] = detail::maximin(a);
  // The column player maximises min_i (−A y)_i: the maximin of −Aᵀ.
  Rectangular neg = a.transpose();
  for (double& t : neg.a) t = -t;
  auto [w, y] = detail::maximin(neg);
  s.lower_value = v;
  s.upper_value = -w;
  s.value = 0.5 * (s.lower_value + s.upper_value);
  s.x = std::move(x);
  s.y = std::move(y);
  return s;
}

}  // namespace nomos::games
