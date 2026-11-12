// Checks of nomos/games/claims against the published examples and closed forms.
//
// References for the values. The Talmud's three estates (Kethubot 93a) with claims 100, 200 and 300, as Aumann and Maschler
// (1985, Table 1) give them: 100 → (33⅓, 33⅓, 33⅓); 200 → (50, 75, 75); 300 → (50, 100, 150). The contested garment
// (Baba Metzia 2a): claims of the whole and of half give (3/4, 1/4). O'Neill (1982): random arrival is the Shapley value of
// the game of what a set secures after the others are paid. The glove game with one left and two right gloves has the value
// (2/3, 1/6, 1/6). The rest are identities held on random problems: efficiency, the bounds, the mean over all orders by brute
// force, and the limits of payment by turns.

#include <nomos/games/claims.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace nomos::games;

namespace {

double worst(const std::vector<double>& a, const std::vector<double>& b) {
  double w = 0;
  for (std::size_t i = 0; i < a.size(); i++) w = std::max(w, std::fabs(a[i] - b[i]));
  return w;
}
double sum(const std::vector<double>& a) { return std::accumulate(a.begin(), a.end(), 0.0); }

// The mean of the priority awards over every order, by brute force.
std::vector<double> mean_over_orders(const Claims& P) {
  const std::size_t n = P.claims.size();
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::vector<double> mean(n, 0.0);
  double count = 0;
  do {
    const auto x = priority(P, order);
    for (std::size_t i = 0; i < n; i++) mean[i] += x[i];
    count++;
  } while (std::next_permutation(order.begin(), order.end()));
  for (double& m : mean) m /= count;
  return mean;
}

}  // namespace

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(20261112, 2, 0, 0);

  // ---- the Shapley value ----
  {
    // One left glove (player 0) and two right gloves: a pair is worth one.
    const auto glove = shapley(3, [](std::uint32_t S) { return (S & 1) && (S & 6) ? 1.0 : 0.0; });
    s.check("Shapley: the glove game with one left and two right gloves, (2/3, 1/6, 1/6)", worst(glove, {2.0 / 3, 1.0 / 6, 1.0 / 6}), 1e-15);
    // An additive game: each player its own worth.
    const std::vector<double> w{3, 1, 4, 1, 5};
    const auto additive = shapley(5, [&](std::uint32_t S) { double v = 0; for (std::size_t i = 0; i < 5; i++) if (S & (1u << i)) v += w[i]; return v; });
    s.check("Shapley: in an additive game each player takes its own worth", worst(additive, w), 1e-14);
    // Efficiency, symmetry, and additivity in the game, on random games of six players.
    double efficient = 0, linear = 0;
    for (int trial = 0; trial < 20; trial++) {
      std::vector<double> u(64), v(64);
      for (std::uint32_t S = 1; S < 64; S++) u[S] = rng.uniform(), v[S] = rng.uniform();
      const auto a = shapley(6, [&](std::uint32_t S) { return u[S]; }), b = shapley(6, [&](std::uint32_t S) { return v[S]; });
      const auto both = shapley(6, [&](std::uint32_t S) { return u[S] + 2 * v[S]; });
      efficient = std::max(efficient, std::fabs(sum(a) - u[63]));
      for (std::size_t i = 0; i < 6; i++) linear = std::max(linear, std::fabs(both[i] - a[i] - 2 * b[i]));
    }
    s.check("Shapley: the values add up to the worth of all, on random games", efficient, 1e-14);
    s.check("Shapley: the value is additive in the game", linear, 1e-14);
    const auto majority = shapley(5, [](std::uint32_t S) { return __builtin_popcount(S) >= 3 ? 1.0 : 0.0; });
    s.check("Shapley: symmetric players share alike (simple majority of five, 1/5 each)", worst(majority, std::vector<double>(5, 0.2)), 1e-15);
  }

  // ---- the Talmud's examples ----
  {
    const std::vector<double> c{100, 200, 300};
    s.check("Talmud: an estate of 100 against claims of 100, 200 and 300 is divided equally", worst(talmud({100, c}), {100.0 / 3, 100.0 / 3, 100.0 / 3}), 1e-13);
    s.check("Talmud: an estate of 200 gives (50, 75, 75)", worst(talmud({200, c}), {50, 75, 75}), 1e-13);
    s.check("Talmud: an estate of 300 gives (50, 100, 150), the proportional division", worst(talmud({300, c}), {50, 100, 150}) + worst(proportional({300, c}), {50, 100, 150}), 1e-13);
    s.check("contested garment: claims of the whole and of half give (3/4, 1/4)", worst(contested_garment(1, 1, 0.5), {0.75, 0.25}), 1e-15);
    double pairwise = 0, two = 0;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 2 + rng.below(4);
      Claims P;
      for (std::size_t i = 0; i < n; i++) P.claims.push_back(rng.uniform(1, 10));
      P.estate = rng.uniform() * sum(P.claims);
      const auto x = talmud(P);
      // Consistency: every pair divides what the two receive as a contested garment.
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = i + 1; j < n; j++) {
          const auto g = contested_garment(x[i] + x[j], std::min(P.claims[i], x[i] + x[j]), std::min(P.claims[j], x[i] + x[j]));
          pairwise = std::max(pairwise, std::fabs(g[0] - x[i]));
        }
      if (n == 2) two = std::max(two, worst(random_arrival(P), x) + worst(contested_garment(P.estate, std::min(P.claims[0], P.estate), std::min(P.claims[1], P.estate)), x));
    }
    s.check("Talmud: every pair divides what it receives as a contested garment (consistency), on 200 problems", pairwise, 1e-12);
    s.check("two claimants: random arrival, the Talmud rule and the contested garment coincide", two, 1e-12);
  }

  // ---- priority and its mean over orders ----
  {
    const Claims P{10, {8, 6}};
    s.check("priority: an estate of 10 against 8 and 6 gives (8, 2) and (4, 6) by the two orders", worst(priority(P, {0, 1}), {8, 2}) + worst(priority(P, {1, 0}), {4, 6}), 1e-15);
    s.check("random arrival: their mean, (6, 4)", worst(random_arrival(P), {6, 4}), 1e-15);
    s.check("random arrival: 7, 5 and 4 against 10 give (4.5, 3, 2.5)", worst(random_arrival({10, {7, 5, 4}}), {4.5, 3, 2.5}), 1e-14);
    double brute = 0, range = 0;
    for (int trial = 0; trial < 60; trial++) {
      const std::size_t n = 2 + rng.below(5);
      Claims Q;
      for (std::size_t i = 0; i < n; i++) Q.claims.push_back(rng.uniform(0, 10));
      Q.estate = rng.uniform() * sum(Q.claims);
      brute = std::max(brute, worst(random_arrival(Q), mean_over_orders(Q)));
      const auto r = priority_range(Q);
      std::vector<std::size_t> order(n);
      std::iota(order.begin(), order.end(), std::size_t{0});
      std::vector<double> lo(n, 1e300), hi(n, -1e300);
      do {
        const auto x = priority(Q, order);
        for (std::size_t i = 0; i < n; i++) lo[i] = std::min(lo[i], x[i]), hi[i] = std::max(hi[i], x[i]);
      } while (std::next_permutation(order.begin(), order.end()));
      range = std::max(range, worst(r.least, lo) + worst(r.greatest, hi));
    }
    s.check("random arrival: the Shapley value of what a set secures is the mean of priority over every order, on 60 problems", brute, 1e-12);
    s.check("priority: the least and the greatest over every order, in closed form", range, 1e-12);
  }

  // ---- every rule divides the estate, within the claims ----
  {
    double efficient = 0, bounded = 0, monotone = 0;
    for (int trial = 0; trial < 200; trial++) {
      const std::size_t n = 2 + rng.below(6);
      Claims P;
      for (std::size_t i = 0; i < n; i++) P.claims.push_back(rng.uniform(0, 10));
      P.estate = rng.uniform() * sum(P.claims);
      for (const auto& x : {random_arrival(P), proportional(P), constrained_equal_awards(P), constrained_equal_losses(P), talmud(P)}) {
        efficient = std::max(efficient, std::fabs(sum(x) - P.estate));
        for (std::size_t i = 0; i < n; i++) {
          bounded = std::max({bounded, -x[i], x[i] - P.claims[i]});
          for (std::size_t j = 0; j < n; j++)
            if (P.claims[i] <= P.claims[j]) monotone = std::max({monotone, x[i] - x[j], (P.claims[i] - x[i]) - (P.claims[j] - x[j])});
        }
      }
    }
    s.check("every rule awards the whole estate, on 200 problems", efficient, 1e-11);
    s.check("every rule awards between nothing and the claim", bounded, 1e-12);
    s.check("every rule gives a greater claim no less, and leaves it no less unpaid", monotone, 1e-12);
    s.check("equal awards: 10 against 8 and 6 gives (5, 5); equal losses gives (6, 4)", worst(constrained_equal_awards({10, {8, 6}}), {5, 5}) + worst(constrained_equal_losses({10, {8, 6}}), {6, 4}), 1e-15);
  }

  // ---- payment by turns ----
  {
    const Claims P{10, {8, 6}};
    const std::vector<double> byClaim{8, 6}, alike{1, 1};
    s.check("turns: whole steps in the two orders are the two priorities", worst(interleaved(P, {0, 1}, byClaim, 1), {8, 2}) + worst(interleaved(P, {1, 0}, byClaim, 1), {4, 6}), 1e-15);
    double gap = 0, last = 0;
    bool falling = true;
    for (double h : {1.0, 1.0 / 4, 1.0 / 16, 1.0 / 64, 1.0 / 256}) {
      const auto a = interleaved(P, {0, 1}, byClaim, h), b = interleaved(P, {1, 0}, byClaim, h);
      gap = std::fabs(a[0] - b[0]);
      if (h < 1 && gap > last) falling = false;
      last = gap;
    }
    s.require("turns: the gap between the two orders of turns falls with the step, from 4 at whole steps", falling);
    s.check("turns: at steps of 1/256 the gap is below a step's payment", gap, 8.0 / 256 + 1e-12);
    const double h = 1e-5;
    s.check("turns: paid at rates in proportion to the claims, the awards tend to the proportional rule (5.714, 4.286)", worst(interleaved(P, {0, 1}, byClaim, h), proportional(P)), 2e-4);
    s.check("turns: paid at equal rates, they tend to the constrained equal awards, another rule", worst(interleaved(P, {0, 1}, alike, h), constrained_equal_awards(P)), 2e-4);
    const Claims Q{12, {2, 9, 7}};
    s.check("turns: with three claimants and equal rates, the small claim is paid in full and the rest shared", worst(interleaved(Q, {2, 0, 1}, {1, 1, 1}, h), constrained_equal_awards(Q)), 2e-4);
  }

  // ---- what is refused ----
  {
    bool threw = false;
    try { (void)proportional({20, {8, 6}}); } catch (const std::invalid_argument&) { threw = true; }
    s.require("an estate that exceeds the claims is refused: there is nothing to divide by a rule", threw);
    threw = false;
    try { (void)shapley(30, [](std::uint32_t) { return 0.0; }); } catch (const std::invalid_argument&) { threw = true; }
    s.require("a game of thirty players is refused", threw);
  }
  return s.done();
}
