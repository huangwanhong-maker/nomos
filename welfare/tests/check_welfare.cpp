// Checks of nomos/welfare: principles as orders over allocations, and what can be said when they disagree.
//
// - Principles: the prioritarian value at ε = 0 the mean and rising towards the minimum as ε grows; the shortfall's share,
//   mean gap and severity by hand; leximin comparing the worst off first; the least mean over types.
// - Dominance: on 300 random pairs, generalised Lorenz dominance implies a higher prioritarian value for every ε tried;
//   a transfer from richer to poorer gives Lorenz dominance; sequential dominance by hand when the neediest group gains.
// - Comparison: an option worse under every principle left out of the maximal set; two undominated options and the pivot
//   between them in words; the inequality aversion at which (1, 9) and (4, 5) exchange places, where their values are equal;
//   rank acceptabilities summing to one for each option and each rank, and an option better under every principle first in
//   every draw.

#include <nomos/welfare/comparison.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::welfare;

int main() {
  calx::verification::Suite s;

  // ---- principles ----
  {
    const Outcome o{"o", {1, 4, 9}, {}};
    const double e0 = evaluate(Principle{"mean", Kind::prioritarian, 0}, o), e1 = evaluate(Principle{"ε=1", Kind::prioritarian, 1}, o),
                 e5 = evaluate(Principle{"ε=5", Kind::prioritarian, 5}, o);
    s.require("the prioritarian value at ε = 0 is the mean, and falls towards the minimum as ε grows",
              std::abs(e0 - 14.0 / 3) < 1e-12 && e1 < e0 && e5 < e1 && e5 > 1 && std::abs(e1 - std::cbrt(36.0)) < 1e-9);
    const Outcome p{"p", {2, 6, 10, 20}, {}};
    const double head = evaluate(Principle{"share below", Kind::shortfall, 8, 0}, p), gap = evaluate(Principle{"gap", Kind::shortfall, 8, 1}, p),
                 sev = evaluate(Principle{"severity", Kind::shortfall, 8, 2}, p);
    s.require("the shortfall below τ = 8: half below, a mean gap of (6/8 + 2/8)/4, a severity of ((6/8)² + (2/8)²)/4",
              std::abs(head - 0.5) < 1e-12 && std::abs(gap - 0.25) < 1e-12 && std::abs(sev - (0.5625 + 0.0625) / 4) < 1e-12);
    s.require("leximin compares the worst off first, then the next",
              compare(Principle{"leximin", Kind::leximin}, Outcome{"a", {2, 2, 10}, {}}, Outcome{"b", {1, 9, 9}, {}}) == 1 &&
                  compare(Principle{"leximin", Kind::leximin}, Outcome{"a", {2, 3, 3}, {}}, Outcome{"b", {2, 2, 9}, {}}) == 1);
    Principle opp{"opportunity", Kind::opportunity};
    opp.types = {0, 0, 1, 1};
    s.require("the opportunity principle takes the least mean over the types", std::abs(evaluate(opp, Outcome{"t", {1, 3, 10, 20}, {}}) - 2) < 1e-12);
  }

  // ---- dominance ----
  {
    calx::random::Stream rng(20261120, 136, 11, 0);
    bool implies = true;
    std::size_t dominations = 0;
    for (int k = 0; k < 300; k++) {
      Outcome a{"a", {}, {}}, b{"b", {}, {}};
      for (int i = 0; i < 5; i++) a.value.push_back(rng.uniform(0.5, 10)), b.value.push_back(rng.uniform(0.5, 10));
      const auto r = generalised_lorenz(a, b);
      if (r != Relation::dominates) continue;
      dominations++;
      for (double eps : {0.0, 0.5, 1.0, 2.0, 5.0})
        implies = implies && compare(Principle{"p", Kind::prioritarian, eps}, a, b) >= 0;
    }
    s.require("generalised Lorenz dominance implies a higher prioritarian value for every ε, on the random pairs where it holds", implies && dominations > 5,
              std::to_string(dominations) + " dominations");
    const Outcome before{"before", {1, 5, 10}, {}}, after{"after", {3, 5, 8}, {}};
    s.require("a transfer from the richest to the poorest, keeping the total, gives Lorenz dominance", lorenz(after, before) == Relation::dominates);
    s.require("and two outcomes whose curves cross are not comparable", lorenz(Outcome{"x", {1, 6, 6}, {}}, Outcome{"y", {2, 2, 9}, {}}) == Relation::unordered);
    // Sequential dominance: group 0 the neediest. a gives the neediest more and the others the same total.
    const Outcome sa{"a", {4, 4, 6, 6}, {}}, sb{"b", {2, 2, 8, 8}, {}};
    s.require("sequential dominance when the neediest group gains and the pooled distribution is less unequal", sequential(sa, sb, {0, 0, 1, 1}) == Relation::dominates);
  }

  // ---- comparison ----
  {
    const std::vector<Outcome> opts{{"A", {1, 9}, {}}, {"B", {4, 5}, {}}, {"C", {1, 5}, {}}};
    const std::vector<Principle> ps{{"total", Kind::total}, {"leximin", Kind::leximin}, {"prioritarian ε=2", Kind::prioritarian, 2}};
    const auto c = compare_all(opts, ps);
    s.require("an option worse than another under every principle is left out of the maximal set", c.maximal == std::vector<std::size_t>{0, 1});
    const std::string pv = pivot(c, 0, 1);
    s.require("the pivot between the two undominated options says which principles favour which",
              pv.find("A is better under total") != std::string::npos && pv.find("B under leximin") != std::string::npos, pv);
    const auto f = flips(Principle{"prioritarian", Kind::prioritarian, 0}, opts[0], opts[1], 0, 10);
    bool equal_at_flip = !f.empty();
    if (!f.empty()) {
      const double va = evaluate(Principle{"p", Kind::prioritarian, f[0]}, opts[0]), vb = evaluate(Principle{"p", Kind::prioritarian, f[0]}, opts[1]);
      equal_at_flip = std::abs(va - vb) < 1e-6;
    }
    s.require("(1, 9) and (4, 5) exchange places at one inequality aversion, where their values are equal", f.size() == 1 && equal_at_flip,
              f.empty() ? "none" : std::to_string(f[0]));
    const auto acc = acceptability(c, {false, false, false}, 4000, 12);
    bool sums = true;
    for (std::size_t i = 0; i < 3; i++) {
      double row = 0, col = 0;
      for (std::size_t r = 0; r < 3; r++) row += acc[i][r], col += acc[r][i];
      sums = sums && std::abs(row - 1) < 1e-9 && std::abs(col - 1) < 1e-9;
    }
    s.require("rank acceptabilities sum to one for each option and for each rank", sums);
    s.require("the option worse under every principle is never first", acc[2][0] == 0);
  }

  return s.done();
}
