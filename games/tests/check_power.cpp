// Checks of simple games and a priori voting power (power.hpp).
//
// - The UN Security Council's rule (nine of fifteen votes, the five permanent members among them), whose swings are counted in
//   closed form: 848 for a permanent member, C(9, 3) = 84 for another, so Banzhaf's index 848/5080 and 84/5080 and Shapley–Shubik's
//   C(9, 3) 8! 6!/15! for a non-permanent member; the same rule as the weighted game [39; 7, 7, 7, 7, 7, 1, …, 1].
// - Nassau County's board of 1964, weights 9, 9, 7, 3, 1, 1 and quota 16, in which Banzhaf (1965) found three members with no power.
// - The two routes held to each other on 300 random weighted games: enumeration of the coalitions against generating functions.
// - A large body: 101 and 301 equal members, whose Penrose measure is C(n − 1, ⌊n/2⌋)/2^(n−1) and Shapley–Shubik 1/n.
// - The properties of a game: vetoes, dummies, minimal winning coalitions, proper, strong, monotone.
// - The square-root rule: in a two-tier system of ten units of sizes k² million, council weights ∝ √m spread citizens' indirect
//   power within 10%, while under weights ∝ m the smallest unit is a dummy; decisive_in_majority against its asymptote √(2/(π m)).

#include <nomos/games/power.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace nomos::games;

namespace {

double factorial(int k) {
  double f = 1;
  for (int i = 2; i <= k; i++) f *= i;
  return f;
}

double max_difference(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (std::size_t i = 0; i < a.size(); i++) m = std::max(m, std::abs(a[i] - b[i]));
  return m;
}

}  // namespace

int main() {
  calx::verification::Suite s;

  // ---- The Security Council --------------------------------------------------------------------------------------------------
  SimpleGame council{15, [](std::uint32_t S) { return (S & 0x1Fu) == 0x1Fu && __builtin_popcount(S) >= 9; }};
  {
    const Power p = power(council);
    s.require("a permanent member has 848 swings, another member C(9, 3) = 84", p.swings[0] == 848 && p.swings[14] == 84);
    s.check("Banzhaf's index of a permanent member: 848/5080", std::abs(p.banzhaf[0] - 848.0 / 5080.0), 1e-15);
    s.check("of a non-permanent member: 84/5080", std::abs(p.banzhaf[7] - 84.0 / 5080.0), 1e-15);
    s.check("the Penrose measure of a permanent member: 848/2¹⁴", std::abs(p.penrose[2] - 848.0 / 16384.0), 1e-15);
    const double phi = 84.0 * factorial(8) * factorial(6) / factorial(15);
    s.check("Shapley–Shubik of a non-permanent member: C(9, 3) 8! 6!/15! = 0.0018648", std::abs(p.shapley_shubik[10] - phi), 1e-15);
    s.check("of a permanent member: (1 − 10 × that)/5 = 0.19627", std::abs(p.shapley_shubik[0] - (1 - 10 * phi) / 5), 1e-14);
    const double W = 848.0 + 0;  // winning coalitions: all five permanent and at least four others, Σ_{k≥4} C(10, k) = 848
    s.check("Coleman's power of the collectivity to act: 848/2¹⁵", std::abs(p.power_to_act - W / 32768.0), 1e-15);
    s.check("a permanent member's power to prevent action is complete: it is in every winning coalition", std::abs(p.coleman_prevent[0] - 1), 1e-15);
    const Properties pr = properties(council);
    s.require("the five permanent members are the vetoes, no member is a dummy, the rule is proper and not strong",
              pr.vetoes == std::vector<std::size_t>{0, 1, 2, 3, 4} && pr.dummies.empty() && pr.proper && !pr.strong && pr.monotone);
    s.require("its minimal winning coalitions are the five with four others: C(10, 4) = 210", pr.minimal_winning.size() == 210);
    WeightedGame weighted{{7, 7, 7, 7, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 39};
    const Power q = power(weighted);
    s.require("as the weighted game [39; 7 × 5, 1 × 10], by generating functions, every measure is the same",
              max_difference(p.penrose, q.penrose) < 1e-15 && max_difference(p.shapley_shubik, q.shapley_shubik) < 1e-15 && max_difference(p.banzhaf, q.banzhaf) < 1e-15 &&
                  std::abs(p.power_to_act - q.power_to_act) < 1e-15 && max_difference(p.coleman_initiate, q.coleman_initiate) < 1e-15);
  }

  // ---- Nassau County, 1964 ------------------------------------------------------------------------------------------------------
  {
    WeightedGame nassau{{9, 9, 7, 3, 1, 1}, 16};
    const Power p = power(nassau);
    const Properties pr = properties(simple(nassau));
    s.require("Nassau County's board: the members of weight 3, 1 and 1 are dummies (Banzhaf 1965)", pr.dummies == std::vector<std::size_t>{3, 4, 5} && p.swings[3] == 0);
    s.require("and the three others share all power equally, by Banzhaf's index and by Shapley–Shubik's",
              std::abs(p.banzhaf[0] - 1.0 / 3) < 1e-15 && std::abs(p.banzhaf[2] - 1.0 / 3) < 1e-15 && std::abs(p.shapley_shubik[2] - 1.0 / 3) < 1e-15);
  }

  // ---- Enumeration against generating functions ----------------------------------------------------------------------------------
  {
    calx::random::Stream rng(20261119, 14, 0, 0);
    double worst = 0;
    for (int k = 0; k < 300; k++) {
      const std::size_t n = 1 + rng.below(14);
      WeightedGame g;
      long total = 0;
      for (std::size_t i = 0; i < n; i++) g.weights.push_back(static_cast<long>(rng.below(12))), total += g.weights.back();
      if (total == 0) g.weights[0] = 1, total = 1;
      g.quota = 1 + static_cast<long>(rng.below(static_cast<std::uint64_t>(total)));
      const Power a = power(simple(g)), b = power(g);
      worst = std::max({worst, max_difference(a.penrose, b.penrose), max_difference(a.shapley_shubik, b.shapley_shubik), max_difference(a.banzhaf, b.banzhaf),
                        max_difference(a.coleman_prevent, b.coleman_prevent), max_difference(a.coleman_initiate, b.coleman_initiate), std::abs(a.power_to_act - b.power_to_act)});
    }
    s.check("on 300 random weighted games the two routes agree on every measure", worst, 1e-13);
  }

  // ---- A large body --------------------------------------------------------------------------------------------------------
  for (std::size_t n : {101u, 301u}) {
    WeightedGame g{std::vector<long>(n, 1), static_cast<long>(n / 2 + 1)};
    const Power p = power(g);
    const double exact = std::exp(calx::special_functions::log_choose(static_cast<double>(n - 1), static_cast<double>(n / 2)) - static_cast<double>(n - 1) * std::log(2.0));
    s.check(std::to_string(n) + " equal members by simple majority: each one's Penrose measure is C(n − 1, ⌊n/2⌋)/2^(n−1)", std::abs(p.penrose[n / 3] - exact) / exact, 1e-12);
    s.check("and its Shapley–Shubik index 1/n", std::abs(p.shapley_shubik[n / 2] - 1.0 / static_cast<double>(n)), 1e-14);
  }

  // ---- Properties ----------------------------------------------------------------------------------------------------------------
  {
    SimpleGame two_of_three{3, [](std::uint32_t S) { return __builtin_popcount(S) >= 2; }};
    s.require("a rule under which any two of three win is monotone", properties(two_of_three).monotone);
    SimpleGame exactly_two{4, [](std::uint32_t S) { return __builtin_popcount(S) == 2 || S == 15u; }};
    s.require("a rule under which exactly two of four win (and all four) is found not monotone", !properties(exactly_two).monotone);
    SimpleGame majority{5, [](std::uint32_t S) { return __builtin_popcount(S) >= 3; }};
    s.require("simple majority of five is decisive (proper and strong), with no veto and no dummy",
              properties(majority).decisive() && properties(majority).vetoes.empty() && properties(majority).dummies.empty());
    // a double majority: a majority of members holding a majority of a declared weight
    SimpleGame members{5, [](std::uint32_t S) { return __builtin_popcount(S) >= 3; }};
    const SimpleGame weight = simple(WeightedGame{{45, 20, 15, 12, 8}, 60});
    const SimpleGame both = all_of({members, weight});
    const Power p = power(both);
    double sum = 0;
    for (double x : p.shapley_shubik) sum += x;
    s.check("a double majority (three of five members holding 60 of 100): Shapley–Shubik's index sums to 1", std::abs(sum - 1), 1e-14);
    s.require("and it is monotone, and the largest member is a veto: the others together hold 55", properties(both).monotone && properties(both).vetoes == std::vector<std::size_t>{0});
  }

  // ---- The square-root rule -----------------------------------------------------------------------------------------------------
  {
    s.check("decisive_in_majority(3) = 1/2", std::abs(decisive_in_majority(3) - 0.5), 1e-15);
    s.check("decisive_in_majority(4) = C(3, 2)/2³ = 3/8", std::abs(decisive_in_majority(4) - 0.375), 1e-15);
    const double m = 1e6 + 1;
    s.check("for a million voters it is within 1/(4m) of √(2/(π m))", std::abs(decisive_in_majority(m) / std::sqrt(2 / (std::numbers::pi * m)) - 1), 1e-6);
    std::vector<double> sizes;
    std::vector<long> root, proportional;
    long troot = 0, tprop = 0;
    for (int k = 1; k <= 10; k++) {
      sizes.push_back(1e6 * k * k);
      root.push_back(k), troot += k;
      proportional.push_back(k * k), tprop += k * k;
    }
    auto indirect = [&](const std::vector<long>& w, long total) { return indirect_power(sizes, power(WeightedGame{w, total / 2 + 1})); };
    const auto by_root = indirect(root, troot), by_size = indirect(proportional, tprop);
    const double sr = *std::max_element(by_root.begin(), by_root.end()) / *std::min_element(by_root.begin(), by_root.end());
    s.require("ten units of 1 to 100 million: with council weights ∝ √m citizens' indirect power is nearly equal (Penrose 1946)", sr < 1.15,
              "max/min " + std::to_string(sr));
    s.require("with weights ∝ m the smallest unit is a dummy, and its citizens have no indirect power at all", by_size[0] == 0 && by_root[0] > 0);
  }
  return s.done();
}
