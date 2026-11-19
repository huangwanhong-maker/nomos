// Checks of bargaining solutions and Adjusted Winner (bargaining.hpp).
//
// - Closed forms: a symmetric triangle (all three at (½, ½)); a linear frontier u₁ + 2u₂ = 4 (Nash and Kalai–Smorodinsky at its
//   midpoint (2, 1), the egalitarian solution at (4/3, 4/3)); a disagreement point inside (each gets d plus half the surplus).
// - Axioms: Nash and Kalai–Smorodinsky invariant under positive affine changes of utility, the egalitarian solution not; Nash
//   independent of irrelevant alternatives, Kalai–Smorodinsky not; Kalai–Smorodinsky individually monotone where Nash is not, on a
//   problem of the kind Kalai and Smorodinsky (1975) gave: T reaches higher than S at every level of player 1, yet Nash gives
//   player 2 less in T (0.7 against 0.75).
// - On 200 random problems: Nash's point maximises the product over 2,000 points of the frontier, and Kalai–Smorodinsky's lies on
//   the segment from d to the ideal point.
// - Adjusted Winner: a division worked by hand; on 300 random problems equitable, envy-free, with at most one issue split, and
//   efficient, held to a linear programme that maximises one party's points while the other keeps its own.

#include <nomos/games/bargaining.hpp>

#include <calx/optimisation/linear_programming.hpp>
#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::games;

namespace {

double gap(const Pair& a, const Pair& b) { return std::max(std::abs(a[0] - b[0]), std::abs(a[1] - b[1])); }

Bargaining transformed(const Bargaining& B, double a1, double b1, double a2, double b2) {
  Bargaining T;
  for (const auto& u : B.outcomes) T.outcomes.push_back({a1 * u[0] + b1, a2 * u[1] + b2});
  T.disagreement = {a1 * B.disagreement[0] + b1, a2 * B.disagreement[1] + b2};
  return T;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- Closed forms ------------------------------------------------------------------------------------------------------------
  {
    Bargaining tri{{{0, 0}, {1, 0}, {0, 1}}, {0, 0}};
    s.require("a symmetric triangle: Nash, Kalai–Smorodinsky and the egalitarian solution all at (½, ½)",
              gap(nash(tri).point, {0.5, 0.5}) < 1e-12 && gap(kalai_smorodinsky(tri).point, {0.5, 0.5}) < 1e-12 && gap(egalitarian(tri).point, {0.5, 0.5}) < 1e-12);
    Bargaining line{{{0, 0}, {4, 0}, {0, 2}}, {0, 0}};
    s.require("a frontier u₁ + 2u₂ = 4: Nash and Kalai–Smorodinsky at its midpoint (2, 1), the egalitarian solution at (4/3, 4/3)",
              gap(nash(line).point, {2, 1}) < 1e-12 && gap(kalai_smorodinsky(line).point, {2, 1}) < 1e-12 && gap(egalitarian(line).point, {4.0 / 3, 4.0 / 3}) < 1e-12);
    s.require("its ideal point is (4, 2)", gap(nash(line).ideal, {4, 2}) < 1e-12);
    Bargaining inside{{{0, 0}, {4, 0}, {0, 4}}, {1, 2}};
    s.require("with d = (1, 2) inside u₁ + u₂ ≤ 4, each gets d and half the surplus: (1.5, 2.5) by all three",
              gap(nash(inside).point, {1.5, 2.5}) < 1e-12 && gap(kalai_smorodinsky(inside).point, {1.5, 2.5}) < 1e-12 && gap(egalitarian(inside).point, {1.5, 2.5}) < 1e-12);
    s.require("the frontier is clipped to u ≥ d: it runs from (1, 3) to (2, 2)", gap(nash(inside).frontier.front(), {1, 3}) < 1e-12 && gap(nash(inside).frontier.back(), {2, 2}) < 1e-12);
  }

  // ---- Axioms ------------------------------------------------------------------------------------------------------------------
  {
    Bargaining B{{{0, 0}, {3, 0}, {2.5, 1.5}, {1, 2.2}, {0, 2.5}}, {0.2, 0.1}};
    const auto t = transformed(B, 2, 3, 0.5, -1);
    auto image = [](const Pair& u) { return Pair{2 * u[0] + 3, 0.5 * u[1] - 1}; };
    s.require("Nash's solution moves with a positive affine change of either party's utility", gap(nash(t).point, image(nash(B).point)) < 1e-12);
    s.require("so does Kalai–Smorodinsky's", gap(kalai_smorodinsky(t).point, image(kalai_smorodinsky(B).point)) < 1e-12);
    s.require("the egalitarian solution does not: it compares the parties' utilities", gap(egalitarian(t).point, image(egalitarian(B).point)) > 1e-3);

    Bargaining big{{{0, 0}, {2, 0}, {0, 2}}, {0, 0}};
    Bargaining sub{{{0, 0}, {1.5, 0.5}, {1, 1}, {0, 1.2}}, {0, 0}};
    s.require("Nash is independent of irrelevant alternatives: (1, 1) in the triangle, and still (1, 1) in a part of it holding (1, 1)",
              gap(nash(big).point, {1, 1}) < 1e-12 && gap(nash(sub).point, {1, 1}) < 1e-12);
    s.require("Kalai–Smorodinsky is not: the part's ideal point moves, and so does its solution", gap(kalai_smorodinsky(big).point, {1, 1}) < 1e-12 && gap(kalai_smorodinsky(sub).point, {1, 1}) > 1e-3);

    Bargaining S{{{0, 0}, {0, 1}, {0.75, 0.75}, {1, 0}}, {0, 0}};
    Bargaining T{{{0, 0}, {0, 1}, {1, 0.7}, {1, 0}}, {0, 0}};
    s.require("T reaches higher than S at every level of player 1, with the same ideal point (1, 1)", gap(nash(S).ideal, {1, 1}) < 1e-12 && gap(nash(T).ideal, {1, 1}) < 1e-12);
    s.require("yet Nash gives player 2 less in T: (1, 0.7) against (0.75, 0.75)", gap(nash(S).point, {0.75, 0.75}) < 1e-12 && gap(nash(T).point, {1, 0.7}) < 1e-12);
    s.require("Kalai–Smorodinsky is individually monotone: (0.75, 0.75) in S, (1/1.3, 1/1.3) in T",
              gap(kalai_smorodinsky(S).point, {0.75, 0.75}) < 1e-12 && gap(kalai_smorodinsky(T).point, {1 / 1.3, 1 / 1.3}) < 1e-12);
  }

  // ---- Random problems --------------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261119, 15, 0, 0);
    double worst_product = 0, worst_ray = 0;
    for (int k = 0; k < 200; k++) {
      Bargaining B;
      const std::size_t m = 3 + rng.below(8);
      for (std::size_t i = 0; i < m; i++) B.outcomes.push_back({rng.uniform(0, 10), rng.uniform(0, 10)});
      B.outcomes.push_back({0, 0});
      B.disagreement = {rng.uniform(0, 2), rng.uniform(0, 2)};
      Bargaining probe = B;
      probe.outcomes.push_back({11, 11});  // so that every problem has room above d
      const auto n = nash(probe);
      const Pair& d = probe.disagreement;
      const double best = (n.point[0] - d[0]) * (n.point[1] - d[1]);
      for (std::size_t i = 0; i + 1 < n.frontier.size(); i++)
        for (int j = 0; j <= 2000 / static_cast<int>(n.frontier.size()); j++) {
          const double t = rng.uniform();
          const Pair u{n.frontier[i][0] + t * (n.frontier[i + 1][0] - n.frontier[i][0]), n.frontier[i][1] + t * (n.frontier[i + 1][1] - n.frontier[i][1])};
          worst_product = std::max(worst_product, (u[0] - d[0]) * (u[1] - d[1]) - best);
        }
      const auto ks = kalai_smorodinsky(probe);
      worst_ray = std::max(worst_ray, std::abs((ks.point[0] - d[0]) * (ks.ideal[1] - d[1]) - (ks.point[1] - d[1]) * (ks.ideal[0] - d[0])));
    }
    s.check("on 200 random problems no point of the frontier has a larger Nash product than Nash's point", std::max(0.0, worst_product), 1e-9);
    s.check("and Kalai–Smorodinsky's point lies on the segment from d to the ideal point", worst_ray, 1e-9);
  }

  // ---- Adjusted Winner -----------------------------------------------------------------------------------------------------------
  {
    const auto D = adjusted_winner({50, 30, 20}, {20, 30, 50});
    s.require("A values the issues 50, 30, 20 and B 20, 30, 50: the tied issue passes half to B, and each holds 65 points",
              D.split == 1 && std::abs(D.share[1] - 0.5) < 1e-12 && D.share[0] == 1 && D.share[2] == 0 && std::abs(D.points_a - 65) < 1e-12 && std::abs(D.points_b - 65) < 1e-12);
    calx::random::Stream rng(20261119, 16, 0, 0);
    double unequal = 0, envy = 0, inefficient = 0;
    std::size_t splits_over = 0;
    for (int k = 0; k < 300; k++) {
      const std::size_t m = 2 + rng.below(9);
      std::vector<double> a(m), b(m);
      double ta = 0, tb = 0;
      for (std::size_t g = 0; g < m; g++) a[g] = rng.uniform(), b[g] = rng.uniform(), ta += a[g], tb += b[g];
      for (std::size_t g = 0; g < m; g++) a[g] *= 100 / ta, b[g] *= 100 / tb;
      const auto R = adjusted_winner(a, b);
      unequal = std::max(unequal, std::abs(R.points_a - R.points_b));
      envy = std::max({envy, 50 - R.points_a, 50 - R.points_b});
      std::size_t split = 0;
      for (double x : R.share) split += (x > 1e-12 && x < 1 - 1e-12);
      splits_over += split > 1;
      // maximise A's points over every division, split or not, while B keeps its own
      calx::optimisation::LinearProgram lp;
      lp.c.assign(m, 0.0);
      for (std::size_t g = 0; g < m; g++) lp.c[g] = -a[g];
      lp.a_ub = calx::linalg::Rectangular(1, m);
      double bt = 0;
      for (std::size_t g = 0; g < m; g++) lp.a_ub(0, g) = b[g], bt += b[g];
      lp.b_ub = {bt - R.points_b};  // B's points Σ (1 − s) b ≥ points_b
      lp.lower.assign(m, 0.0), lp.upper.assign(m, 1.0);
      const auto res = calx::optimisation::simplex(lp);
      inefficient = std::max(inefficient, -res.value - R.points_a);
    }
    s.check("on 300 random problems the division is equitable", unequal, 1e-9);
    s.check("envy-free: each holds at least half of its own 100 points", std::max(0.0, envy), 1e-9);
    s.require("with at most one issue split", splits_over == 0);
    s.check("and efficient: no division gives A more while B keeps its points", std::max(0.0, inefficient), 1e-7);
  }
  return s.done();
}
