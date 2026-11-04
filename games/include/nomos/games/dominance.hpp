// Iterated elimination of strictly dominated strategies.
//
// A pure strategy i of the row player is strictly dominated when another strategy gives a
// strictly higher payoff against every remaining strategy of the opponent: another pure strategy
// k (A_kj > A_ij for every remaining j), or, more generally, a mixture σ of the other remaining
// rows (Σ_k σ_k A_kj > A_ij for every remaining j). A strictly dominated strategy is never a best
// response, so it is used in no Nash equilibrium, and removing it changes no equilibrium of what
// remains; strategies are removed for both players until none is dominated. For strict dominance
// the result does not depend on the order of removal (Gilboa, Kalai and Zemel 1990; Osborne and
// Rubinstein 1994, §4.2). A game is dominance solvable when one cell survives.
//
// Mixed dominance is decided by the linear programme
//   maximise ε  subject to  Σ_k σ_k A_kj ≥ A_ij + ε  (j remaining),  Σ_k σ_k = 1,  σ ≥ 0,
// k ranging over the remaining rows other than i: i is dominated by a mixture exactly when the
// optimum ε* is positive (in two-player games, exactly when i is never a best response; Pearce
// 1984; the programme as in Shoham and Leyton-Brown 2009, ch. 4). It is solved by calx's simplex
// method. Columns are treated alike with B.
//
// Strictness: a payoff must exceed the other by more than `tolerance` (default 10⁻⁹) to count.

#pragma once

#include <calx/optimisation/linear_programming.hpp>

#include <limits>
#include <vector>

#include "bimatrix.hpp"

namespace nomos::games {

enum class Dominance { pure, mixed };

// What survives, and what was removed in which order.
struct Reduction {
  struct Step {
    bool row_player = true;
    std::size_t strategy = 0;
  };
  std::vector<std::size_t> rows, cols;  // the surviving strategies, increasing
  std::vector<Step> removed;            // in the order of removal
};

namespace detail {

// Is strategy `i` of the player whose payoffs are P(own, other) dominated, among `own`, against `other`?
template <class Payoff>
bool dominated(std::size_t i, const std::vector<std::size_t>& own, const std::vector<std::size_t>& other, Payoff P, Dominance kind,
               double tolerance) {
  for (std::size_t k : own) {
    if (k == i) continue;
    bool better = true;
    for (std::size_t j : other) better = better && P(k, j) > P(i, j) + tolerance;
    if (better) return true;
  }
  if (kind == Dominance::pure || own.size() < 3) return false;  // with two rows, mixed dominance is pure dominance
  // Variables: σ_k for the other rows, then ε (free). Minimise −ε.
  std::vector<std::size_t> others;
  for (std::size_t k : own)
    if (k != i) others.push_back(k);
  const std::size_t nv = others.size() + 1;
  calx::optimisation::LinearProgram lp;
  lp.c.assign(nv, 0.0);
  lp.c.back() = -1;
  lp.a_ub = calx::linalg::Rectangular(other.size(), nv);
  lp.b_ub.assign(other.size(), 0.0);
  for (std::size_t r = 0; r < other.size(); r++) {
    for (std::size_t c = 0; c < others.size(); c++) lp.a_ub(r, c) = -P(others[c], other[r]);
    lp.a_ub(r, nv - 1) = 1;
    lp.b_ub[r] = -P(i, other[r]);
  }
  lp.a_eq = calx::linalg::Rectangular(1, nv);
  for (std::size_t c = 0; c < others.size(); c++) lp.a_eq(0, c) = 1;
  lp.b_eq = {1};
  lp.lower.assign(nv, 0.0);
  lp.lower.back() = -std::numeric_limits<double>::infinity();
  lp.upper.assign(nv, std::numeric_limits<double>::infinity());
  const auto r = calx::optimisation::simplex(lp);
  return r.status == calx::optimisation::LpStatus::optimal && -r.value > tolerance;
}

}  // namespace detail

inline Reduction eliminate_dominated(const Bimatrix& g, Dominance kind = Dominance::pure, double tolerance = 1e-9) {
  Reduction red;
  for (std::size_t i = 0; i < g.rows(); i++) red.rows.push_back(i);
  for (std::size_t j = 0; j < g.cols(); j++) red.cols.push_back(j);
  const auto row_payoff = [&](std::size_t i, std::size_t j) { return g.a(i, j); };
  const auto col_payoff = [&](std::size_t j, std::size_t i) { return g.b(i, j); };
  for (bool changed = true; changed;) {
    changed = false;
    for (std::size_t t = 0; t < red.rows.size() && red.rows.size() > 1; t++)
      if (detail::dominated(red.rows[t], red.rows, red.cols, row_payoff, kind, tolerance)) {
        red.removed.push_back({true, red.rows[t]});
        red.rows.erase(red.rows.begin() + static_cast<std::ptrdiff_t>(t));
        changed = true;
        break;
      }
    if (changed) continue;
    for (std::size_t t = 0; t < red.cols.size() && red.cols.size() > 1; t++)
      if (detail::dominated(red.cols[t], red.cols, red.rows, col_payoff, kind, tolerance)) {
        red.removed.push_back({false, red.cols[t]});
        red.cols.erase(red.cols.begin() + static_cast<std::ptrdiff_t>(t));
        changed = true;
        break;
      }
  }
  return red;
}

}  // namespace nomos::games
