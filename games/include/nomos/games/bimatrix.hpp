// Finite two-player games in normal form: payoffs, best responses and pure equilibria.
//
// A bimatrix game (A, B) gives the row player, who chooses a row i of m, the payoff A_ij, and the
// column player, who chooses a column j of n, the payoff B_ij (von Neumann and Morgenstern 1944;
// von Stengel 2007, §3.2). A mixed strategy is a probability vector x on the rows or y on the
// columns; the expected payoffs are xᵀAy and xᵀBy. Against y, the row player's best responses are
// the rows maximising (Ay)_i; against x, the column player's maximise (xᵀB)_j. A pair (x, y) is a
// Nash equilibrium (Nash 1951) when each is a best response to the other, which holds exactly
// when every pure strategy in the support of x is a best response to y and every pure strategy
// in the support of y a best response to x: the best response condition (von Stengel 2007,
// Proposition 3.1). A pure equilibrium (i, j) is a cell in which A_ij is largest in its column
// and B_ij largest in its row.
//
// A zero-sum game has B = −A; a symmetric game has B = Aᵀ.
//
// Ties: payoffs within `tolerance` of the best count as best (default 10⁻⁹, for payoffs of order
// one; exact integer payoffs tie exactly).

#pragma once

#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::games {

using calx::linalg::Rectangular;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::games: ") + what);
}
}  // namespace detail

// The payoff matrices of the row player (a) and of the column player (b), both m × n.
struct Bimatrix {
  Rectangular a, b;

  Bimatrix() = default;
  Bimatrix(Rectangular row_payoffs, Rectangular column_payoffs) : a(std::move(row_payoffs)), b(std::move(column_payoffs)) {
    detail::require(a.rows > 0 && a.cols > 0 && a.rows == b.rows && a.cols == b.cols, "two payoff matrices of one shape");
  }
  static Bimatrix zero_sum(const Rectangular& a) {
    Rectangular b = a;
    for (double& v : b.a) v = -v;
    return {a, b};
  }
  static Bimatrix symmetric(const Rectangular& a) {
    detail::require(a.rows == a.cols, "a square matrix for a symmetric game");
    return {a, a.transpose()};
  }
  std::size_t rows() const { return a.rows; }
  std::size_t cols() const { return a.cols; }
};

// xᵀAy and xᵀBy.
inline double row_payoff(const Bimatrix& g, std::span<const double> x, std::span<const double> y) {
  detail::require(x.size() == g.rows() && y.size() == g.cols(), "a strategy for each player");
  const auto Ay = calx::linalg::multiply(g.a, y);
  double v = 0;
  for (std::size_t i = 0; i < x.size(); i++) v += x[i] * Ay[i];
  return v;
}
inline double column_payoff(const Bimatrix& g, std::span<const double> x, std::span<const double> y) {
  detail::require(x.size() == g.rows() && y.size() == g.cols(), "a strategy for each player");
  const auto By = calx::linalg::multiply(g.b, y);
  double v = 0;
  for (std::size_t i = 0; i < x.size(); i++) v += x[i] * By[i];
  return v;
}

namespace detail {
inline std::vector<std::size_t> argmax(const std::vector<double>& v, double tolerance) {
  double best = v[0];
  for (double t : v) best = std::max(best, t);
  std::vector<std::size_t> out;
  for (std::size_t i = 0; i < v.size(); i++)
    if (v[i] >= best - tolerance) out.push_back(i);
  return out;
}
}  // namespace detail

// The rows maximising (Ay)_i.
inline std::vector<std::size_t> best_responses_row(const Bimatrix& g, std::span<const double> y, double tolerance = 1e-9) {
  detail::require(y.size() == g.cols(), "a column strategy");
  return detail::argmax(calx::linalg::multiply(g.a, y), tolerance);
}

// The columns maximising (xᵀB)_j.
inline std::vector<std::size_t> best_responses_column(const Bimatrix& g, std::span<const double> x, double tolerance = 1e-9) {
  detail::require(x.size() == g.rows(), "a row strategy");
  return detail::argmax(calx::linalg::multiply_transpose(g.b, x), tolerance);
}

// A cell (row, col).
struct Cell {
  std::size_t row = 0, col = 0;
  bool operator==(const Cell&) const = default;
};

// The pure equilibria, in row-major order. O(mn).
inline std::vector<Cell> pure_equilibria(const Bimatrix& g, double tolerance = 1e-9) {
  const std::size_t m = g.rows(), n = g.cols();
  std::vector<double> col_best(n, -INFINITY), row_best(m, -INFINITY);
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < n; j++) {
      col_best[j] = std::max(col_best[j], g.a(i, j));
      row_best[i] = std::max(row_best[i], g.b(i, j));
    }
  std::vector<Cell> out;
  for (std::size_t i = 0; i < m; i++)
    for (std::size_t j = 0; j < n; j++)
      if (g.a(i, j) >= col_best[j] - tolerance && g.b(i, j) >= row_best[i] - tolerance) out.push_back({i, j});
  return out;
}

}  // namespace nomos::games
