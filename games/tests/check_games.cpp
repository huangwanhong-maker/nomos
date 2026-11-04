// Checks of finite games (bimatrix.hpp, dominance.hpp, equilibria.hpp, zero_sum.hpp,
// normal_form.hpp).
//
// The prisoner's dilemma (Gibbons 1992, §1.1.A): one pure equilibrium, reached by dominance.
// Gibbons' dominance-solvable game (§1.1.B): Right, then Down, then Left removed, (Up, Middle)
// left. Matching pennies: no pure equilibrium, the mixed one (½, ½), value 0. The battle of the
// sexes: two pure equilibria and the mixed one ((⅔, ⅓), (⅓, ⅔)) with payoffs ⅔ each. Rock, paper,
// scissors: value 0 and uniform strategies. A 2 × 2 zero-sum game without a saddle point: value
// (ad − bc)/(a − b − c + d) = ⅕ and x = (⅖, ⅗). The 3 × 2 game of von Stengel (2007, §3.2) with
// its three equilibria. A row dominated only by a mixture. Forty random nondegenerate games: every
// equilibrium found satisfies the best response condition, their number is odd (Wilson 1971), and
// those of support one are the pure equilibria. Random zero-sum games: the two linear programmes
// meet (the minimax theorem) and agree with support enumeration. Three players: a coordination
// game with two pure equilibria, and Jordan's (1993) game of matching and mismatching with none.

#include <nomos/games/bimatrix.hpp>
#include <nomos/games/dominance.hpp>
#include <nomos/games/equilibria.hpp>
#include <nomos/games/normal_form.hpp>
#include <nomos/games/zero_sum.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

namespace ng = nomos::games;
using calx::linalg::Rectangular;

namespace {

double distance(const std::vector<double>& a, const std::vector<double>& b) {
  double d = 0;
  for (std::size_t i = 0; i < a.size(); i++) d = std::max(d, std::fabs(a[i] - b[i]));
  return a.size() == b.size() ? d : INFINITY;
}

// The largest distance from each expected equilibrium to the nearest one found (∞ if counts differ).
double match(const std::vector<ng::Equilibrium>& found, const std::vector<std::pair<std::vector<double>, std::vector<double>>>& expected) {
  if (found.size() != expected.size()) return INFINITY;
  double worst = 0;
  for (const auto& [x, y] : expected) {
    double best = INFINITY;
    for (const auto& e : found) best = std::min(best, std::max(distance(e.x, x), distance(e.y, y)));
    worst = std::max(worst, best);
  }
  return worst;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream z(20261024, 0x67616d65u, 0, 0);  // "game"

  // The prisoner's dilemma: strategies Mum (0) and Fink (1).
  {
    const ng::Bimatrix g(Rectangular::from_rows({{-1, -9}, {0, -6}}), Rectangular::from_rows({{-1, 0}, {-9, -6}}));
    const auto pure = ng::pure_equilibria(g);
    s.require("the prisoner's dilemma: one pure equilibrium, (Fink, Fink)", pure.size() == 1 && pure[0] == ng::Cell{1, 1});
    const auto r = ng::eliminate_dominated(g);
    s.require("... reached by eliminating the dominated Mum of each", r.rows == std::vector<std::size_t>{1} && r.cols == std::vector<std::size_t>{1});
    const auto all = ng::support_enumeration(g);
    s.require("... and support enumeration finds no other", all.size() == 1 && all[0].row_value == -6 && all[0].column_value == -6);
  }
  // Gibbons' Figure 1.1.1: rows Up, Down; columns Left, Middle, Right.
  {
    const ng::Bimatrix g(Rectangular::from_rows({{1, 1, 0}, {0, 0, 2}}), Rectangular::from_rows({{0, 2, 1}, {3, 1, 0}}));
    const auto r = ng::eliminate_dominated(g);
    s.require("a dominance-solvable game: (Up, Middle) survives", r.rows == std::vector<std::size_t>{0} && r.cols == std::vector<std::size_t>{1});
    s.require("... after Right, Down and Left, in that order",
              r.removed.size() == 3 && !r.removed[0].row_player && r.removed[0].strategy == 2 && r.removed[1].row_player && r.removed[1].strategy == 1 &&
                  !r.removed[2].row_player && r.removed[2].strategy == 0);
    const auto pure = ng::pure_equilibria(g);
    s.require("... which is its only pure equilibrium", pure.size() == 1 && pure[0] == ng::Cell{0, 1});
  }
  // Matching pennies.
  {
    const Rectangular A = Rectangular::from_rows({{1, -1}, {-1, 1}});
    const auto g = ng::Bimatrix::zero_sum(A);
    s.require("matching pennies: no pure equilibrium", ng::pure_equilibria(g).empty());
    s.check("... one mixed equilibrium, (½, ½) against (½, ½)", match(ng::support_enumeration(g), {{{0.5, 0.5}, {0.5, 0.5}}}), 1e-14);
    const auto v = ng::solve_zero_sum(A);
    s.check("... of value 0, by linear programming", std::fabs(v.value) + distance(v.x, {0.5, 0.5}) + distance(v.y, {0.5, 0.5}), 1e-12);
    s.require("... and the best responses to ½ are both strategies", ng::best_responses_row(g, std::vector<double>{0.5, 0.5}).size() == 2);
  }
  // The battle of the sexes.
  {
    const ng::Bimatrix g(Rectangular::from_rows({{2, 0}, {0, 1}}), Rectangular::from_rows({{1, 0}, {0, 2}}));
    const auto all = ng::support_enumeration(g);
    s.check("the battle of the sexes: three equilibria, two pure and ((⅔, ⅓), (⅓, ⅔))",
            match(all, {{{1, 0}, {1, 0}}, {{0, 1}, {0, 1}}, {{2.0 / 3, 1.0 / 3}, {1.0 / 3, 2.0 / 3}}}), 1e-14);
    double payoff = INFINITY;
    for (const auto& e : all)
      if (e.x[0] > 0 && e.x[0] < 1) payoff = std::max(std::fabs(e.row_value - 2.0 / 3), std::fabs(e.column_value - 2.0 / 3));
    s.check("... the mixed one paying ⅔ to each", payoff, 1e-14);
    const auto pure = ng::pure_equilibria(ng::normal_form(g));
    s.require("... and its normal form has the same two pure equilibria",
              pure.size() == 2 && pure[0] == std::vector<std::size_t>{0, 0} && pure[1] == std::vector<std::size_t>{1, 1});
  }
  // Rock, paper, scissors.
  {
    const Rectangular A = Rectangular::from_rows({{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}});
    const auto v = ng::solve_zero_sum(A);
    const std::vector<double> third(3, 1.0 / 3);
    s.check("rock, paper, scissors: value 0, uniform strategies", std::fabs(v.value) + distance(v.x, third) + distance(v.y, third), 1e-12);
    s.check("... and support enumeration finds the uniform equilibrium alone", match(ng::support_enumeration(ng::Bimatrix::zero_sum(A)), {{third, third}}), 1e-14);
  }
  // A 2 × 2 zero-sum game without a saddle point.
  {
    const auto v = ng::solve_zero_sum(Rectangular::from_rows({{2, -1}, {-1, 1}}));
    s.check("a 2 × 2 zero-sum game: value ⅕", v.value - 0.2, 1e-12);
    s.check("... x = y = (⅖, ⅗)", std::max(distance(v.x, {0.4, 0.6}), distance(v.y, {0.4, 0.6})), 1e-12);
    s.check("... the two programmes meet (the minimax theorem)", v.upper_value - v.lower_value, 1e-12);
  }
  // von Stengel's 3 × 2 game with three equilibria.
  {
    const ng::Bimatrix g(Rectangular::from_rows({{3, 3}, {2, 5}, {0, 6}}), Rectangular::from_rows({{3, 2}, {2, 6}, {3, 1}}));
    s.check("von Stengel's 3 × 2 game: its three equilibria ((1, 0, 0), (1, 0)), ((⅘, ⅕, 0), (⅔, ⅓)), ((0, ⅓, ⅔), (⅓, ⅔))",
            match(ng::support_enumeration(g), {{{1, 0, 0}, {1, 0}}, {{0.8, 0.2, 0}, {2.0 / 3, 1.0 / 3}}, {{0, 1.0 / 3, 2.0 / 3}, {1.0 / 3, 2.0 / 3}}}), 1e-14);
  }
  // A row dominated only by a mixture.
  {
    const ng::Bimatrix g(Rectangular::from_rows({{3, 0}, {0, 3}, {1, 1}}), Rectangular::from_rows({{1, 0}, {0, 1}, {1, 1}}));
    const auto pure = ng::eliminate_dominated(g, ng::Dominance::pure);
    const auto mixed = ng::eliminate_dominated(g, ng::Dominance::mixed);
    s.require("a row dominated only by ½(first + second): no pure strategy dominates it", pure.removed.empty());
    s.require("... the linear programme finds the mixture, and only that row goes",
              mixed.removed.size() == 1 && mixed.removed[0].row_player && mixed.removed[0].strategy == 2 && mixed.cols.size() == 2);
  }
  // Random nondegenerate games.
  {
    bool all_ok = true, odd = true, pure_agree = true;
    std::size_t total = 0;
    for (int t = 0; t < 40; t++) {
      const std::size_t m = t < 20 ? 4 : 3, n = t < 20 ? 4 : 5;
      Rectangular A(m, n), B(m, n);
      for (double& v : A.a) v = z.uniform(-1, 1);
      for (double& v : B.a) v = z.uniform(-1, 1);
      const ng::Bimatrix g(A, B);
      const auto all = ng::support_enumeration(g);
      total += all.size();
      odd = odd && all.size() % 2 == 1;
      std::vector<ng::Cell> from_supports;
      for (const auto& e : all) {
        all_ok = all_ok && ng::is_equilibrium(g, e.x, e.y, 1e-10);
        std::size_t kx = 0, ky = 0, i0 = 0, j0 = 0;
        for (std::size_t i = 0; i < m; i++)
          if (e.x[i] > 0) kx++, i0 = i;
        for (std::size_t j = 0; j < n; j++)
          if (e.y[j] > 0) ky++, j0 = j;
        if (kx == 1 && ky == 1) from_supports.push_back({i0, j0});
      }
      pure_agree = pure_agree && from_supports == ng::pure_equilibria(g);
    }
    s.require("forty random games (4 × 4 and 3 × 5): every equilibrium found satisfies the best response condition", all_ok, std::to_string(total) + " equilibria");
    s.require("... each game has an odd number of them (Wilson 1971)", odd);
    s.require("... and those of support one are exactly the pure equilibria", pure_agree);
  }
  // Random zero-sum games.
  {
    double gap = 0, guarantee = 0, agree = 0;
    for (int t = 0; t < 10; t++) {
      Rectangular A(5, 7);
      for (double& v : A.a) v = z.uniform(-2, 2);
      const auto v = ng::solve_zero_sum(A);
      gap = std::max(gap, std::fabs(v.upper_value - v.lower_value));
      const auto xA = calx::linalg::multiply_transpose(A, v.x), Ay = calx::linalg::multiply(A, v.y);
      for (double u : xA) guarantee = std::max(guarantee, v.value - u);  // x guarantees at least v
      for (double u : Ay) guarantee = std::max(guarantee, u - v.value);  // y concedes at most v
      for (const auto& e : ng::support_enumeration(ng::Bimatrix::zero_sum(A))) agree = std::max(agree, std::fabs(e.row_value - v.value));
    }
    s.check("ten random 5 × 7 zero-sum games: max-min = min-max (the minimax theorem)", gap, 1e-10);
    s.check("... x guarantees the value against every column, y concedes no more against any row", guarantee, 1e-10);
    s.check("... and every equilibrium found by support enumeration pays the value", agree, 1e-10);
  }
  // Three players.
  {
    ng::NormalForm g({2, 2, 2});
    for (std::size_t k = 0; k < g.profiles(); k++) {
      const auto p = g.profile(k);
      const double all_same = p[0] == p[1] && p[1] == p[2] ? 1 : 0;
      for (std::size_t q = 0; q < 3; q++) g.payoffs[k * 3 + q] = all_same;
    }
    const auto eq = ng::pure_equilibria(g);
    s.require("a three-player coordination game: (0, 0, 0) and (1, 1, 1)",
              eq.size() == 2 && eq[0] == std::vector<std::size_t>{0, 0, 0} && eq[1] == std::vector<std::size_t>{1, 1, 1});
    ng::NormalForm j({2, 2, 2});
    for (std::size_t k = 0; k < j.profiles(); k++) {
      const auto p = j.profile(k);
      j.payoffs[k * 3 + 0] = p[0] == p[1] ? 1 : 0;  // the first matches the second
      j.payoffs[k * 3 + 1] = p[1] == p[2] ? 1 : 0;  // the second matches the third
      j.payoffs[k * 3 + 2] = p[2] != p[0] ? 1 : 0;  // the third mismatches the first
    }
    s.require("Jordan's game of matching and mismatching: no pure equilibrium", ng::pure_equilibria(j).empty());
    bool round = true;
    ng::NormalForm h({3, 4, 2, 5});
    for (std::size_t k = 0; k < h.profiles(); k++) round = round && h.index(h.profile(k)) == k;
    s.require("profiles and their numbers correspond one to one (3 × 4 × 2 × 5)", round && h.profiles() == 120);
  }
  return s.done();
}
