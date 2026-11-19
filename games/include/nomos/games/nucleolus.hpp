// The core, the least core and the nucleolus of a coalitional game.
//
// A coalitional game on n players gives a worth v(S) to every coalition, v(∅) = 0 (as in claims.hpp). A pre-imputation is an x with
// x(N) = v(N); the excess of a coalition at x is e(S, x) = v(S) − x(S), what S could gain by leaving. Then:
//
//   core          the pre-imputations at which no coalition has a positive excess (Gillies 1959); it is non-empty exactly when the
//                 game is balanced (Bondareva 1963; Shapley 1967), which a linear programme decides;
//   least core    the pre-imputations minimising the largest excess ε₁ over the proper coalitions (Maschler, Peleg and Shapley 1979);
//                 the core is non-empty exactly when ε₁ ≤ 0;
//   prenucleolus  the pre-imputation whose vector of excesses, sorted from the largest, is lexicographically least (Schmeidler 1969;
//                 Sobolev 1975). It exists and is unique. When it is an imputation (x_i ≥ v({i}) for every i), as it is for every
//                 game in which a coalition is worth at least the sum of its parts, it is the nucleolus, the same least point over
//                 the imputations: a lexicographic minimum over a set that lies in a smaller set is the minimum over the smaller one.
//
// The prenucleolus is found by a sequence of linear programmes (Kopelowitz 1967; Maschler, Peleg and Shapley 1979): minimise ε
// subject to e(S, x) ≤ ε over the coalitions still free and the excesses already fixed held at their levels. Every coalition with a
// positive dual in the programme's optimum has its excess at ε in every optimal solution (complementary slackness), and at least one
// has, since the duals of the free coalitions sum to one; those are fixed at that level. Free coalitions whose excess the fixed ones
// determine (their incidence vectors in the span of the fixed ones and of N) are dropped. The rounds stop when the fixed coalitions
// span ℝⁿ, which determines x. The programmes have 2ⁿ rows; n ≤ 10 here, on calx's simplex.
//
// The nucleolus of a claims problem's game w(S) = max(0, E − Σ_{j ∉ S} c_j) is the Talmud rule (Aumann and Maschler 1985), which
// the checks hold on random problems.
//
// References. Gillies, D. B. (1959), Solutions to general non-zero-sum games, Contributions to the Theory of Games IV, 47–85.
// Bondareva, O. N. (1963), Some applications of linear programming methods to the theory of cooperative games, Problemy Kibernetiki
// 10, 119–139. Shapley, L. S. (1967), On balanced sets and cores, Naval Research Logistics Quarterly 14, 453–460. Schmeidler, D.
// (1969), The nucleolus of a characteristic function game, SIAM J. Applied Math. 17, 1163–1170. Kopelowitz, A. (1967), Computation of
// the kernels of simple games and the nucleolus of n-person games, Research Memorandum 31, Hebrew University. Maschler, M., Peleg, B.
// and Shapley, L. S. (1979), Geometric properties of the kernel, nucleolus and related solution concepts, Math. Oper. Res. 4,
// 303–338. Aumann, R. J. and Maschler, M. (1985), Game theoretic analysis of a bankruptcy problem from the Talmud, J. Econ. Theory 36,
// 195–213. Peleg, B. and Sudhölter, P. (2007), Introduction to the Theory of Cooperative Games, 2nd ed., Springer, ch. 5–6.

#pragma once

#include <calx/linalg/rectangular.hpp>
#include <calx/optimisation/linear_programming.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nomos::games {

namespace detail {

inline std::vector<double> worths(std::size_t n, const std::function<double(std::uint32_t)>& v) {
  if (n == 0 || n > 10) throw std::invalid_argument("nomos::games: the core and the nucleolus are computed for 1 to 10 players");
  const std::uint32_t sets = std::uint32_t{1} << n;
  std::vector<double> w(sets);
  for (std::uint32_t S = 0; S < sets; S++) w[S] = v(S);
  if (w[0] != 0) throw std::invalid_argument("nomos::games: the empty coalition is worth 0");
  return w;
}

inline double sum(const std::vector<double>& x, std::uint32_t S) {
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++)
    if (S >> i & 1u) s += x[i];
  return s;
}

// Rows of 0/1 incidence vectors kept in reduced echelon form, to tell whether a coalition's vector lies in their span.
struct Span {
  std::size_t n;
  std::vector<std::vector<double>> rows;
  std::vector<std::size_t> pivots;
  std::vector<double> reduce(std::uint32_t S) const {
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; i++) r[i] = (S >> i & 1u) ? 1.0 : 0.0;
    for (std::size_t k = 0; k < rows.size(); k++) {
      const double f = r[pivots[k]];
      if (f != 0)
        for (std::size_t i = 0; i < n; i++) r[i] -= f * rows[k][i];
    }
    return r;
  }
  bool contains(std::uint32_t S) const {
    for (double x : reduce(S))
      if (std::abs(x) > 1e-9) return false;
    return true;
  }
  void add(std::uint32_t S) {
    auto r = reduce(S);
    std::size_t p = n;
    for (std::size_t i = 0; i < n; i++)
      if (std::abs(r[i]) > 1e-9) {
        p = i;
        break;
      }
    if (p == n) return;
    const double f = r[p];
    for (double& x : r) x /= f;
    for (auto& row : rows) {
      const double g = row[p];
      if (g != 0)
        for (std::size_t i = 0; i < n; i++) row[i] -= g * r[i];
    }
    rows.push_back(r), pivots.push_back(p);
  }
  std::size_t rank() const { return rows.size(); }
};

}  // namespace detail

struct LeastCore {
  double epsilon = 0;     // the least largest excess over the proper coalitions
  std::vector<double> x;  // a point of the least core
  bool core_nonempty() const { return epsilon <= 1e-9; }
};

// The least core's value and a point of it: minimise ε subject to v(S) − x(S) ≤ ε for every proper coalition and x(N) = v(N).
inline LeastCore least_core(std::size_t n, const std::function<double(std::uint32_t)>& v) {
  const auto w = detail::worths(n, v);
  const std::uint32_t all = (std::uint32_t{1} << n) - 1;
  calx::optimisation::LinearProgram lp;
  lp.c.assign(n + 1, 0.0), lp.c[n] = 1;
  const std::size_t rows = all - 1;
  lp.a_ub = calx::linalg::Rectangular(rows, n + 1);
  lp.b_ub.assign(rows, 0.0);
  std::size_t r = 0;
  for (std::uint32_t S = 1; S < all; S++, r++) {
    for (std::size_t i = 0; i < n; i++)
      if (S >> i & 1u) lp.a_ub(r, i) = -1;
    lp.a_ub(r, n) = -1, lp.b_ub[r] = -w[S];
  }
  lp.a_eq = calx::linalg::Rectangular(1, n + 1);
  for (std::size_t i = 0; i < n; i++) lp.a_eq(0, i) = 1;
  lp.b_eq = {w[all]};
  const double inf = std::numeric_limits<double>::infinity();
  lp.lower.assign(n + 1, -inf), lp.upper.assign(n + 1, inf);
  LeastCore out;
  if (n == 1) return out.x = {w[all]}, out.epsilon = -inf, out;
  const auto res = calx::optimisation::simplex(lp);
  if (res.status != calx::optimisation::LpStatus::optimal) throw std::runtime_error("nomos::games::least_core: the programme was not solved");
  out.epsilon = res.value;
  out.x.assign(res.x.begin(), res.x.begin() + static_cast<std::ptrdiff_t>(n));
  return out;
}

// Whether x is in the core: efficient, and no proper coalition with a positive excess (to the tolerance).
inline bool in_core(const std::vector<double>& x, const std::function<double(std::uint32_t)>& v, double tolerance = 1e-9) {
  const std::size_t n = x.size();
  const auto w = detail::worths(n, v);
  const std::uint32_t all = (std::uint32_t{1} << n) - 1;
  if (std::abs(detail::sum(x, all) - w[all]) > tolerance * (1 + std::abs(w[all]))) return false;
  for (std::uint32_t S = 1; S < all; S++)
    if (w[S] - detail::sum(x, S) > tolerance * (1 + std::abs(w[S]))) return false;
  return true;
}

// The excesses at x of every proper coalition, sorted from the largest: the vector the nucleolus makes lexicographically least.
inline std::vector<double> excesses(const std::vector<double>& x, const std::function<double(std::uint32_t)>& v) {
  const std::size_t n = x.size();
  const auto w = detail::worths(n, v);
  const std::uint32_t all = (std::uint32_t{1} << n) - 1;
  std::vector<double> e;
  for (std::uint32_t S = 1; S < all; S++) e.push_back(w[S] - detail::sum(x, S));
  std::sort(e.begin(), e.end(), std::greater<>());
  return e;
}

struct Nucleolus {
  std::vector<double> x;
  std::vector<double> levels;  // the excesses fixed round by round, decreasing; levels[0] is the least core's ε
  std::size_t rounds = 0;
  bool imputation = false;     // x_i ≥ v({i}) for every i: then x is the nucleolus, not only the prenucleolus
};

inline Nucleolus prenucleolus(std::size_t n, const std::function<double(std::uint32_t)>& v) {
  const auto w = detail::worths(n, v);
  const std::uint32_t all = (std::uint32_t{1} << n) - 1;
  Nucleolus out;
  if (n == 1) {
    out.x = {w[all]}, out.imputation = true;
    return out;
  }
  detail::Span span{n, {}, {}};
  span.add(all);
  std::vector<std::uint32_t> free;
  for (std::uint32_t S = 1; S < all; S++) free.push_back(S);
  std::vector<std::pair<std::uint32_t, double>> fixed;  // coalition and its excess
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> x;
  while (span.rank() < n) {
    // drop free coalitions whose excess the fixed ones determine
    free.erase(std::remove_if(free.begin(), free.end(), [&](std::uint32_t S) { return span.contains(S); }), free.end());
    if (free.empty()) break;
    calx::optimisation::LinearProgram lp;
    lp.c.assign(n + 1, 0.0), lp.c[n] = 1;
    lp.a_ub = calx::linalg::Rectangular(free.size(), n + 1);
    lp.b_ub.assign(free.size(), 0.0);
    for (std::size_t r = 0; r < free.size(); r++) {
      for (std::size_t i = 0; i < n; i++)
        if (free[r] >> i & 1u) lp.a_ub(r, i) = -1;
      lp.a_ub(r, n) = -1, lp.b_ub[r] = -w[free[r]];
    }
    lp.a_eq = calx::linalg::Rectangular(1 + fixed.size(), n + 1);
    lp.b_eq.assign(1 + fixed.size(), 0.0);
    for (std::size_t i = 0; i < n; i++) lp.a_eq(0, i) = 1;
    lp.b_eq[0] = w[all];
    for (std::size_t k = 0; k < fixed.size(); k++) {
      for (std::size_t i = 0; i < n; i++)
        if (fixed[k].first >> i & 1u) lp.a_eq(1 + k, i) = 1;
      lp.b_eq[1 + k] = w[fixed[k].first] - fixed[k].second;  // x(S) = v(S) − e
    }
    lp.lower.assign(n + 1, -inf), lp.upper.assign(n + 1, inf);
    const auto res = calx::optimisation::simplex(lp);
    if (res.status != calx::optimisation::LpStatus::optimal) throw std::runtime_error("nomos::games::prenucleolus: a programme was not solved");
    const double eps = res.value;
    out.levels.push_back(eps), out.rounds++;
    x.assign(res.x.begin(), res.x.begin() + static_cast<std::ptrdiff_t>(n));
    // the coalitions with a positive dual (a negative marginal ∂value/∂b) are at ε in every optimum
    double scale = 1;
    for (double d : res.duals_ub) scale = std::max(scale, std::abs(d));
    std::vector<std::uint32_t> now;
    for (std::size_t r = 0; r < free.size(); r++)
      if (res.duals_ub[r] < -1e-9 * scale) now.push_back(free[r]);
    if (now.empty()) {
      // degenerate duals: fall back on the coalition that is tightest at the optimum
      std::size_t best = 0;
      double most = -inf;
      for (std::size_t r = 0; r < free.size(); r++) {
        const double e = w[free[r]] - detail::sum(x, free[r]);
        if (e > most) most = e, best = r;
      }
      now.push_back(free[best]);
    }
    for (std::uint32_t S : now) {
      fixed.push_back({S, eps});
      span.add(S);
    }
    free.erase(std::remove_if(free.begin(), free.end(), [&](std::uint32_t S) { return std::find(now.begin(), now.end(), S) != now.end(); }), free.end());
  }
  out.x = x;
  out.imputation = true;
  for (std::size_t i = 0; i < n; i++)
    if (x[i] < w[std::uint32_t{1} << i] - 1e-9 * (1 + std::abs(w[std::uint32_t{1} << i]))) out.imputation = false;
  return out;
}

}  // namespace nomos::games
