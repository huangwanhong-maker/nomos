// Checks of nomos/allocation: who first, what and how much when there is not enough.
//
// - Rationing: the TAL family at θ = ½ is the Talmud rule on the three estates of the Talmud, at 0 and 1 the equal losses and
//   awards; Piniles's rule by hand; adjusted proportional equal to the contested garment for two claimants; the duals of equal
//   awards and equal losses each other, proportional and Talmud self-dual, and every rule feasible, on random problems;
//   weighted proportional capped at the claims; classes served in order.
// - Formulas: sufficiency first closing every shortfall in plenty and rationing them in scarcity (the acequias' switch); banded
//   points shared in proportion; equalisation's grants summing to zero, by hand for two units, and an envelope's floor kept.
// - Reserves: the sequential reserve matching eligible, non-wasteful and respecting priorities on random instances, and one of
//   the matchings that satisfy the three, found by enumerating them all; hard reserves processed open first leaving the
//   reserved unit idle (Pathak et al.'s Example A.2), reserve first using it; the beneficiaries under every precedence order.
// - Assignment: probabilistic serial giving two agents of one preference half of each object; rows and columns of random
//   assignments within their bounds; probabilistic serial envy-free in stochastic dominance on random instances; the random
//   serial dictatorship's exact probabilities against sampling; a random assignment written as a lottery that sums back.
// - Division: round-robin and envy-cycle elimination envy-free up to one good, and the maximum Nash welfare EF1, on random
//   additive instances; a division that is not proportional named.
// - Budgeting: a majority of six and a minority of four with ten unit projects each: utilitarian greedy funds only the
//   majority's and breaks EJR; the Method of Equal Shares funds six and four, satisfies EJR, and spends what it is given.

#include <nomos/allocation/budgeting.hpp>
#include <nomos/allocation/division.hpp>
#include <nomos/allocation/rationing.hpp>
#include <nomos/allocation/reserves.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace nomos::allocation;
namespace g = nomos::games;

static bool close(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-9) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); i++)
    if (std::abs(a[i] - b[i]) > tol * std::max(1.0, std::abs(b[i]))) return false;
  return true;
}

static bool feasible(const Claims& P, const std::vector<double>& x) {
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    if (x[i] < -1e-9 || x[i] > P.claims[i] + 1e-9) return false;
    s += x[i];
  }
  return std::abs(s - P.estate) <= 1e-7 * std::max(1.0, P.estate);
}

int main() {
  calx::verification::Suite s;

  // ---- rationing ----
  {
    const std::vector<double> c{100, 200, 300};
    bool talmud = true;
    for (double E : {100.0, 200.0, 300.0, 450.0}) talmud = talmud && close(tal(Claims{E, c}, 0.5), g::talmud(Claims{E, c}), 1e-7);
    s.require("the TAL family at θ = ½ is the Talmud rule, on the Talmud's three estates and beyond", talmud);
    s.require("the Talmud's estates of 100, 200 and 300 give (33⅓, 33⅓, 33⅓), (50, 75, 75) and (50, 100, 150)",
              close(tal(Claims{100, c}, 0.5), {100.0 / 3, 100.0 / 3, 100.0 / 3}, 1e-7) && close(tal(Claims{200, c}, 0.5), {50, 75, 75}, 1e-7) &&
                  close(tal(Claims{300, c}, 0.5), {50, 100, 150}, 1e-7));
    s.require("at θ = 0 it is constrained equal losses, at θ = 1 constrained equal awards",
              close(tal(Claims{250, c}, 0), g::constrained_equal_losses(Claims{250, c}), 1e-7) && close(tal(Claims{250, c}, 1), g::constrained_equal_awards(Claims{250, c}), 1e-7));
    s.require("Piniles's rule: equal awards on the half-claims below C/2, the half-claims plus equal awards above",
              close(piniles(Claims{100, c}), {100.0 / 3, 100.0 / 3, 100.0 / 3}, 1e-7) && close(piniles(Claims{400, c}), {50 + 100.0 / 3, 100 + 100.0 / 3, 150 + 100.0 / 3}, 1e-7));
    s.require("adjusted proportional is the contested garment for two claimants",
              close(adjusted_proportional(Claims{1, {1, 0.5}}), g::contested_garment(1, 1, 0.5), 1e-9) && close(adjusted_proportional(Claims{0.8, {1, 0.5}}), g::contested_garment(0.8, 1, 0.5), 1e-9));
    calx::random::Stream rng(20261120, 136, 1, 0);
    bool duals = true, self_dual = true, all_feasible = true;
    const Rule cea = g::constrained_equal_awards, cel = g::constrained_equal_losses, pro = g::proportional, tlm = g::talmud;
    for (int k = 0; k < 200; k++) {
      Claims P{0, {}};
      const std::size_t n = 2 + rng.below(5);
      double C = 0;
      for (std::size_t i = 0; i < n; i++) P.claims.push_back(rng.uniform(0.1, 10)), C += P.claims.back();
      P.estate = rng.uniform(0, C);
      duals = duals && close(dual(cea)(P), cel(P), 1e-6) && close(dual(cel)(P), cea(P), 1e-6);
      self_dual = self_dual && close(dual(pro)(P), pro(P), 1e-9) && close(dual(tlm)(P), tlm(P), 1e-6);
      for (const auto& x : {piniles(P), adjusted_proportional(P), tal(P, rng.uniform()), weighted_proportional(P, std::vector<double>(n, 1.0))})
        all_feasible = all_feasible && feasible(P, x);
    }
    s.require("on 200 random problems the dual of equal awards is equal losses and the reverse", duals);
    s.require("proportional and Talmud are their own duals", self_dual);
    s.require("every rule's awards lie between zero and the claims and sum to the estate", all_feasible);
    const auto w = weighted_proportional(Claims{10, {2, 20, 20}}, {10, 1, 1});
    s.require("weighted proportional caps a claimant at its claim and passes the excess on", close(w, {2, 4, 4}, 1e-9));
    const auto pc = priority_classes(Claims{25, {10, 10, 20, 20}}, {0, 0, 1, 1}, g::proportional);
    s.require("classes served in order: the first in full, the second shares what is left", close(pc, {10, 10, 2.5, 2.5}, 1e-9));
  }

  // ---- formulas ----
  {
    const auto plenty = sufficiency_first(100, {10, 30, 0}, {1, 1, 2}, Rule(g::proportional));
    s.require("in plenty every shortfall is closed and the rest shared by the weights", close(plenty.first, {10, 30, 0}) && close(plenty.second, {15, 15, 30}) && !plenty.scarce);
    const auto scarce = sufficiency_first(20, {10, 30, 0}, {1, 1, 2}, Rule(g::constrained_equal_awards));
    s.require("in scarcity the shortfalls are rationed (equal awards here) and nothing goes to the weights: the acequias' switch",
              close(scarce.first, {10, 10, 0}) && close(scarce.second, {0, 0, 0}) && scarce.scarce);
    const auto pts = banded_points(1000, {Indicator{"need", {5, 50}, {10, 40}, true, 4}, Indicator{"population", {1000, 10000}, {5000}, true, 2}});
    s.require("banded points: grades from the bands, points the weighted grades, shares in proportion",
              close(pts.points, {4 * 1 + 2 * 1, 4 * 3 + 2 * 2}) && std::abs(pts.award[0] + pts.award[1] - 1000) < 1e-9 && std::abs(pts.share[1] - 16.0 / 22) < 1e-12);
    std::vector<Unit> units{{"a", 100, 10, {50}, 0}, {"b", 300, 6, {30}, 0}};
    const auto eq = equalisation(units, {0.2}, 1.0);
    // ē = (100·10 + 300·6)/400 = 7; b̄ = (100·50 + 300·30)/400 = 35; a: 100[(10 − 7) + 0.2(35 − 50)] = 0; b: 300[(6 − 7) + 0.2(35 − 30)] = 0.
    s.require("equalisation by hand: needs above the standard offset by capacity above it", std::abs(eq.first[0]) < 1e-9 && std::abs(eq.first[1]) < 1e-9);
    std::vector<Unit> u2{{"a", 100, 12, {50}, 0}, {"b", 300, 6, {30}, 0}, {"c", 200, 8, {40}, 0}};
    const auto e2 = equalisation(u2, {0.2}, 0.8);
    s.require("equalisation's grants sum to zero: a redistribution among the units", std::abs(e2.first[0] + e2.first[1] + e2.first[2]) < 1e-9);
    for (auto& u : u2) u.previous = 50;
    const auto e3 = equalisation(u2, {0.2}, 1.0, Envelope{0.0, 1e300, 1.0});
    bool floor = true;
    for (std::size_t i = 0; i < 3; i++) floor = floor && e3.total[i] >= 50 - 1e-9 && std::abs(e3.total[i] - e3.first[i] - e3.second[i]) < 1e-9;
    s.require("an envelope with a floor of no loss keeps every unit at its previous grant or above, the adjustment reported apart", floor);
  }

  // ---- reserves ----
  {
    calx::random::Stream rng(20261120, 136, 2, 0);
    bool axioms = true, among = true;
    for (int k = 0; k < 150; k++) {
      const std::size_t n = 3 + rng.below(3);
      std::vector<Category> cats;
      const std::size_t C = 2 + rng.below(2);
      for (std::size_t c = 0; c < C; c++) {
        Category cat{"c" + std::to_string(c), 1 + rng.below(2), {}};
        std::vector<std::size_t> people(n);
        for (std::size_t i = 0; i < n; i++) people[i] = i;
        for (std::size_t i = n; i > 1; i--) std::swap(people[i - 1], people[rng.below(i)]);
        for (std::size_t i : people)
          if (rng.uniform() < 0.75) cat.priority.push_back(i);
        cats.push_back(cat);
      }
      std::vector<std::size_t> order(C);
      for (std::size_t c = 0; c < C; c++) order[c] = c;
      const auto m = sequential_reserve(n, cats, order);
      const auto r = check_axioms(n, cats, m.category_of);
      axioms = axioms && r.eligible && r.non_wasteful && r.respects_priorities;
      // Enumerate every matching and keep those satisfying the three axioms: the sequential one must be among them.
      std::vector<std::size_t> mu(n, 0);
      bool found = false;
      for (;;) {
        std::vector<std::size_t> as(n);
        for (std::size_t i = 0; i < n; i++) as[i] = mu[i] == C ? unserved : mu[i];
        const auto rr = check_axioms(n, cats, as);
        if (rr.eligible && rr.non_wasteful && rr.respects_priorities && as == m.category_of) found = true;
        std::size_t q = 0;
        while (q < n && ++mu[q] == C + 1) mu[q++] = 0;
        if (q == n) break;
      }
      among = among && found;
    }
    s.require("the sequential reserve matching is eligible, non-wasteful and respects priorities on 150 random instances", axioms);
    s.require("and is one of the matchings satisfying the three axioms, found by enumerating them all", among);
    // Hard reserves (Example A.2): patients 0 ≻ 1; an open unit for both, a reserved unit for patient 0 alone.
    std::vector<Category> hard{{"open", 1, {0, 1}}, {"reserve", 1, {0}}};
    const auto open_first = sequential_reserve(2, hard, {0, 1}), reserve_first = sequential_reserve(2, hard, {1, 0});
    s.require("hard reserves processed open first leave the reserved unit idle and patient 1 unserved",
              open_first.category_of[0] == 0 && open_first.category_of[1] == unserved && open_first.holders[1].empty());
    s.require("processed reserve first, both are served", reserve_first.category_of[0] == 1 && reserve_first.category_of[1] == 0);
    s.require("yet the open-first matching also satisfies the three axioms: they do not exclude waste of this kind",
              [&] { const auto r = check_axioms(2, hard, open_first.category_of); return r.eligible && r.non_wasteful && r.respects_priorities; }());
    const auto every = every_order(2, hard);
    s.require("the beneficiaries under every precedence order are reported", every.size() == 2 && every[0].served.size() + every[1].served.size() == 3);
  }

  // ---- assignment ----
  {
    const Preferences same{{0, 1}, {0, 1}};
    const auto ps = probabilistic_serial(same, {1, 1});
    s.require("two agents of one preference each get half of each object under probabilistic serial", close(ps.p[0], {0.5, 0.5}) && close(ps.p[1], {0.5, 0.5}));
    calx::random::Stream rng(20261120, 136, 3, 0);
    bool bounds = true, sd_envy_free = true, sums_back = true, exact_vs_sampled = true;
    for (int k = 0; k < 60; k++) {
      const std::size_t n = 3 + rng.below(3), m = 3;
      Preferences prefs(n);
      for (auto& p : prefs) {
        p = {0, 1, 2};
        for (std::size_t i = m; i > 1; i--) std::swap(p[i - 1], p[rng.below(i)]);
      }
      std::vector<std::size_t> copies{1 + rng.below(2), 1, 1};
      const auto a = probabilistic_serial(prefs, copies), b = random_serial(prefs, copies);
      for (const auto* r : {&a, &b}) {
        for (const auto& row : r->p) {
          double t = 0;
          for (double v : row) t += v;
          bounds = bounds && t <= 1 + 1e-9;
        }
        for (std::size_t o = 0; o < m; o++) {
          double t = 0;
          for (const auto& row : r->p) t += row[o];
          bounds = bounds && t <= static_cast<double>(copies[o]) + 1e-9;
        }
      }
      // sd envy-freeness: for each agent i and other j, i's cumulative share of its top objects ≥ j's of the same objects.
      for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++) {
          double ci = 0, cj = 0;
          for (std::size_t o : prefs[i]) {
            ci += a.p[i][o], cj += a.p[j][o];
            if (ci < cj - 1e-9) sd_envy_free = false;
          }
        }
      const auto sq = square(a, copies);
      const auto lot = birkhoff(sq.matrix);
      std::vector<std::vector<double>> back(sq.matrix.size(), std::vector<double>(sq.matrix.size(), 0));
      double total = 0;
      for (const auto& d : lot) {
        total += d.weight;
        for (std::size_t r = 0; r < d.object_of.size(); r++) back[r][d.object_of[r]] += d.weight;
      }
      for (std::size_t r = 0; r < back.size(); r++) sums_back = sums_back && close(back[r], sq.matrix[r], 1e-7);
      sums_back = sums_back && std::abs(total - 1) < 1e-7;
    }
    {
      // Exact against sampled: the sampler forced by more than eight agents is held to an enumeration on a cloned eight.
      Preferences prefs(8, {0, 1, 2});
      for (std::size_t a = 0; a < 8; a++)
        if (a % 2) prefs[a] = {1, 0, 2};
      const std::vector<std::size_t> copies{2, 2, 2};
      const auto exact = random_serial(prefs, copies);
      Preferences nine = prefs;
      nine.push_back({2, 1, 0});
      calx::random::Stream rs(20261120, 136, 5, 0);
      const auto sampled = random_serial(nine, {2, 2, 2}, &rs, 40000);
      // Every agent accepts every object and there are six copies for eight: an agent is served unless it is among the last two
      // of a uniformly random order, a chance of 6/8 exactly.
      exact_vs_sampled = exact.exact && !sampled.exact;
      for (const auto& row : exact.p) {
        double t = 0;
        for (double v : row) t += v;
        exact_vs_sampled = exact_vs_sampled && std::abs(t - 0.75) < 1e-12;
      }
      double rowsum = 0;
      for (double v : sampled.p[8]) rowsum += v;
      // The ninth agent loses its first choice only when six or more agents precede it (both copies of objects 0 and 1, then
      // both of 2, gone): a chance of 3/9, so it receives object 2 with probability 2/3 and nothing otherwise.
      exact_vs_sampled = exact_vs_sampled && std::abs(sampled.p[8][2] - 2.0 / 3) < 0.01 && std::abs(rowsum - 2.0 / 3) < 0.01;
    }
    s.require("random assignments keep each agent's chances within one and each object's within its copies", bounds);
    s.require("probabilistic serial is envy-free in stochastic dominance on 60 random instances", sd_envy_free);
    s.require("a random assignment made square and written as a lottery over assignments sums back to it", sums_back);
    s.require("the random serial dictatorship is exact over all orders for eight agents (each served with chance 6/8), and sampled beyond gives the ninth agent its first choice with the chance 2/3 that counting gives", exact_vs_sampled);
  }

  // ---- division ----
  {
    calx::random::Stream rng(20261120, 136, 6, 0);
    bool rr = true, ec = true, mnw = true;
    for (int k = 0; k < 80; k++) {
      const std::size_t n = 2 + rng.below(2), m = 4 + rng.below(3);
      Values v(n, std::vector<double>(m));
      for (auto& row : v)
        for (double& x : row) x = std::floor(rng.uniform(0, 10));
      rr = rr && check_fairness(v, round_robin(v), false).ef1;
      ec = ec && check_fairness(v, envy_cycle(v), false).ef1;
      mnw = mnw && check_fairness(v, max_nash_welfare(v), false).ef1;
    }
    s.require("round-robin is envy-free up to one good on 80 random additive instances", rr);
    s.require("envy-cycle elimination is envy-free up to one good", ec);
    s.require("the maximum Nash welfare is envy-free up to one good", mnw);
    const Values v{{1, 1, 1, 1}, {1, 1, 1, 1}};
    const auto rep = check_fairness(v, {0, 0, 0, 1});
    s.require("a division of three goods against one is named not proportional and not EF1, with its witnesses",
              !rep.proportional && !rep.ef1 && !rep.mms && !rep.witnesses.empty() && std::abs(rep.mms_share[0] - 2) < 1e-12);
  }

  // ---- budgeting ----
  {
    Election e;
    e.budget = 10;
    e.costs.assign(20, 1.0);
    for (int i = 0; i < 6; i++) e.approvals.push_back({0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    for (int i = 0; i < 4; i++) e.approvals.push_back({10, 11, 12, 13, 14, 15, 16, 17, 18, 19});
    const auto greedy = utilitarian_greedy(e), mes = equal_shares(e, Satisfaction::cardinality, false), ph = phragmen(e);
    auto count_minority = [](const Selection& sel) {
      std::size_t k = 0;
      for (std::size_t j : sel.funded) k += j >= 10;
      return k;
    };
    s.require("utilitarian greedy gives the majority of six every project and the minority of four none", greedy.funded.size() == 10 && count_minority(greedy) == 0);
    s.require("which breaks extended justified representation: the four could have afforded four projects", !check_ejr(e, greedy.funded).satisfied && check_ejr(e, greedy.funded).projects.size() == 4 && check_ejr(e, greedy.funded).group.size() == 4);
    s.require("the Method of Equal Shares funds six for the six and four for the four, spending the budget", mes.funded.size() == 10 && count_minority(mes) == 4 && std::abs(mes.spent - 10) < 1e-9);
    s.require("and satisfies EJR", check_ejr(e, mes.funded).satisfied);
    s.require("sequential Phragmén also funds four for the four", count_minority(ph) == 4 && ph.funded.size() == 10);
    double paid = 0;
    for (double x : mes.paid) paid += x;
    s.require("what the voters paid under equal shares is what was spent", std::abs(paid - mes.spent) < 1e-9);
  }

  return s.done();
}
