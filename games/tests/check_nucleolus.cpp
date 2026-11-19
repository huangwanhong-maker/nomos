// Checks of the core, the least core and the nucleolus (nucleolus.hpp).
//
// - The glove game (one left glove, two right): the core is the single point (1, 0, 0), and so is the nucleolus.
// - The majority game of three: the core is empty, the least core's value is 1/3, the nucleolus (1/3, 1/3, 1/3).
// - A claims problem's game: its nucleolus is the Talmud rule (Aumann and Maschler 1985), on their estate of 200 against claims of
//   100, 200 and 300, and on 150 random problems of 3 to 6 claimants.
// - Convex games: the core is non-empty, and holds the nucleolus and the Shapley value (Shapley 1971).
// - The nucleolus is lexicographically least: on 40 random games of four players, no nearby pre-imputation has a lexicographically
//   smaller vector of excesses.
// - Symmetric players receive the same.

#include <nomos/games/claims.hpp>
#include <nomos/games/nucleolus.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::games;

namespace {

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::abs(a[i] - b[i]));
  return m;
}

// a < b lexicographically, beyond a tolerance
bool lex_less(const std::vector<double>& a, const std::vector<double>& b, double tol) {
  for (std::size_t i = 0; i < a.size(); i++) {
    if (a[i] < b[i] - tol) return true;
    if (a[i] > b[i] + tol) return false;
  }
  return false;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  {
    auto glove = [](std::uint32_t S) { return (S & 1u) && (S & 6u) ? 1.0 : 0.0; };
    const auto lc = least_core(3, glove);
    s.require("the glove game: its core is not empty (the least core's value is 0)", lc.core_nonempty() && std::abs(lc.epsilon) < 1e-12);
    const auto nu = prenucleolus(3, glove);
    s.require("and its nucleolus is (1, 0, 0), the core's one point", max_difference(nu.x, {1, 0, 0}) < 1e-12 && in_core(nu.x, glove) && nu.imputation);
  }
  {
    auto majority = [](std::uint32_t S) { return __builtin_popcount(S) >= 2 ? 1.0 : 0.0; };
    const auto lc = least_core(3, majority);
    s.require("the majority game of three: its core is empty, the least core's value 1/3", !lc.core_nonempty() && std::abs(lc.epsilon - 1.0 / 3) < 1e-12);
    s.require("its nucleolus is (1/3, 1/3, 1/3)", max_difference(prenucleolus(3, majority).x, {1.0 / 3, 1.0 / 3, 1.0 / 3}) < 1e-12);
  }

  // ---- Claims and the Talmud ---------------------------------------------------------------------------------------------------
  {
    Claims P{200, {100, 200, 300}};
    const auto nu = prenucleolus(3, [&](std::uint32_t S) { return secured(P, S); });
    s.require("an estate of 200 against claims of 100, 200 and 300: the nucleolus is the Talmud's (50, 75, 75)", max_difference(nu.x, {50, 75, 75}) < 1e-9);
    calx::random::Stream rng(20261119, 17, 0, 0);
    double worst = 0;
    std::size_t imputations = 0;
    for (int k = 0; k < 150; k++) {
      const std::size_t n = 3 + rng.below(4);
      Claims Q;
      double total = 0;
      for (std::size_t i = 0; i < n; i++) Q.claims.push_back(rng.uniform(1, 100)), total += Q.claims.back();
      Q.estate = rng.uniform(0.05, 0.95) * total;
      const auto r = prenucleolus(n, [&](std::uint32_t S) { return secured(Q, S); });
      worst = std::max(worst, max_difference(r.x, talmud(Q)) / Q.estate);
      imputations += r.imputation;
    }
    s.check("on 150 random problems of 3 to 6 claimants the nucleolus of the claims game is the Talmud rule", worst, 1e-9);
    s.require("and it is always an imputation, so the prenucleolus is the nucleolus", imputations == 150);
  }

  // ---- Convex games -------------------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261119, 18, 0, 0);
    std::size_t good = 0;
    for (int k = 0; k < 30; k++) {
      const std::size_t n = 3 + rng.below(4);
      std::vector<double> a(n);
      for (auto& x : a) x = rng.uniform(0.5, 2);
      auto v = [&](std::uint32_t S) {
        double t = 0;
        for (std::size_t i = 0; i < n; i++)
          if (S >> i & 1u) t += a[i];
        return t * t;  // a convex function of an additive measure: a convex game
      };
      const auto lc = least_core(n, v);
      const auto nu = prenucleolus(n, v);
      const auto phi = shapley(n, v);
      good += lc.core_nonempty() && in_core(nu.x, v, 1e-8) && in_core(phi, v, 1e-8) && nu.imputation;
    }
    s.require("on 30 convex games the core is not empty and holds both the nucleolus and the Shapley value", good == 30, std::to_string(good) + " of 30");
  }

  // ---- Lexicographic least ------------------------------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261119, 19, 0, 0);
    std::size_t beaten = 0, games = 0;
    for (int k = 0; k < 40; k++) {
      const std::size_t n = 4;
      std::vector<double> w(16);
      for (std::uint32_t S = 1; S < 16; S++) w[S] = rng.uniform(0, 1) * __builtin_popcount(S);
      auto v = [&](std::uint32_t S) { return w[S]; };
      const auto nu = prenucleolus(n, v);
      const auto e0 = excesses(nu.x, v);
      for (int t = 0; t < 300; t++) {
        std::vector<double> x = nu.x;
        double mean = 0;
        std::vector<double> delta(n);
        for (auto& d : delta) d = rng.uniform(-1, 1), mean += d / n;
        const double size = std::pow(10.0, -1 - 3 * rng.uniform());
        for (std::size_t i = 0; i < n; i++) x[i] += size * (delta[i] - mean);
        if (lex_less(excesses(x, v), e0, 1e-12)) beaten++;
      }
      games++;
    }
    s.require("on 40 random games no nearby pre-imputation has a lexicographically smaller vector of excesses than the nucleolus", beaten == 0,
              std::to_string(games) + " games, 12,000 perturbations" + (beaten ? ", " + std::to_string(beaten) + " beaten" : ""));
  }

  // ---- Symmetry ----------------------------------------------------------------------------------------------------------------
  {
    auto v = [](std::uint32_t S) {
      const int a = (S & 1u) + ((S >> 1) & 1u), b = ((S >> 2) & 1u) + ((S >> 3) & 1u) + ((S >> 4) & 1u);
      return a >= 1 && b >= 2 ? 6.0 : (a + b >= 3 ? 2.0 : 0.0);
    };
    const auto nu = prenucleolus(5, v);
    s.require("players 1 and 2 are symmetric, and so are 3, 4 and 5: each group's members receive the same",
              std::abs(nu.x[0] - nu.x[1]) < 1e-9 && std::abs(nu.x[2] - nu.x[3]) < 1e-9 && std::abs(nu.x[3] - nu.x[4]) < 1e-9);
    bool decreasing = nu.rounds == nu.levels.size();
    std::string levels;
    for (std::size_t r = 0; r < nu.levels.size(); r++) {
      levels += (r ? ", " : "") + std::to_string(nu.levels[r]);
      if (r && nu.levels[r] > nu.levels[r - 1] + 1e-9) decreasing = false;
    }
    s.require("its rounds fix levels of excess that do not increase", decreasing, levels);
  }
  return s.done();
}
