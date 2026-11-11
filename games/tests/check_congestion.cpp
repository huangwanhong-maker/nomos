// Checks of routing games with congestion (congestion.hpp).
//
// - Pigou's two roads (a constant hour, and one whose time is its flow): everyone takes the second at equilibrium (total 1),
//   half do at the optimum (total 3/4), the price of anarchy 4/3; with the second road's time x⁴, the price of anarchy
//   1 / (1 − 4 · 5^(−5/4)) (Roughgarden and Tardos 2002).
// - Braess's network of 4,000 drivers: 2,000 on each route at 65 minutes; with a road of no time added, all 4,000 on the
//   new route at 80 minutes, every journey longer.
// - Two parallel roads of the Bureau of Public Roads: the equilibrium where their times are equal, found by bisection.
// - A grid of three by three crossings with five demands: the arc flows at Wardrop's equilibrium and at the system optimum,
//   the total times and Beckmann's potential, against scipy's path-based solution (fixtures/make_congestion_fixtures.py);
//   the gap of scipy's flows is zero by the code here; Frank–Wolfe at a gap of 10⁻⁵ is held within the distance its own
//   bound allows, its lower bound below the least potential; and the skim at the equilibrium's times is the least of all
//   the paths' times.
// - Forty random networks with affine times: the optimum's total time at most the equilibrium's, and their ratio at most 4/3.

#include <nomos/games/congestion.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "fixtures/congestion_fixtures.hpp"

namespace ng = nomos::games;

int main() {
  calx::verification::Suite s;
  {
    ng::Network pigou{2, {{0, 1, 1, 0, 1}, {0, 1, 0, 1, 1}}, {{0, 1, 1}}};
    const auto ue = ng::assign(pigou, {ng::Objective::user, 1e-12});
    const auto so = ng::assign(pigou, {ng::Objective::system, 1e-12, 1000000});
    s.check("Pigou: everyone takes the road whose time is its flow, a total of 1", std::fabs(ue.flow[1] - 1) + std::fabs(ue.total_time - 1), 1e-12);
    s.check("Pigou: half take it at the optimum, a total of 3/4", std::fabs(so.flow[1] - 0.5) + std::fabs(so.total_time - 0.75), 1e-9);
    s.check("Pigou: the price of anarchy is 4/3", ng::price_of_anarchy(pigou) - 4.0 / 3, 1e-9);
    ng::Network quartic{2, {{0, 1, 1, 0, 1}, {0, 1, 0, 1, 4}}, {{0, 1, 1}}};
    s.check("Pigou with x⁴: the price of anarchy is 1 / (1 − 4 · 5^(−5/4))", ng::price_of_anarchy(quartic, 1e-12) - 1 / (1 - 4 * std::pow(5.0, -1.25)), 1e-8);
  }
  {
    // S = 0, A = 1, B = 2, E = 3.
    ng::Network braess{4, {{0, 1, 0, 0.01, 1}, {1, 3, 45, 0, 1}, {0, 2, 45, 0, 1}, {2, 3, 0, 0.01, 1}}, {{0, 3, 4000}}};
    const auto before = ng::assign(braess, {ng::Objective::user, 1e-12});
    s.check("Braess: 2,000 drivers on each route", std::fabs(before.flow[0] - 2000) + std::fabs(before.flow[2] - 2000), 1e-6);
    s.check("Braess: each journey takes 65 minutes", std::fabs(before.time[0] + before.time[1] - 65) + std::fabs(before.time[2] + before.time[3] - 65), 1e-9);
    braess.arcs.push_back({1, 2, 0, 0, 1});
    const auto after = ng::assign(braess, {ng::Objective::user, 1e-10, 1000000});
    s.check("Braess: with the new road, all 4,000 take it", std::fabs(after.flow[4] - 4000) / 4000, 1e-6);
    s.check("Braess: and every journey takes 80 minutes", std::fabs(after.total_time / 4000 - 80), 1e-4);
    s.require("Braess: the road added lengthens the total time from 260,000 to 320,000 minutes", before.total_time < after.total_time);
  }
  {
    ng::Network two{2, {ng::bpr(0, 1, 10, 1000), ng::bpr(0, 1, 15, 3000)}, {{0, 1, 4000}}};
    double lo = 0, hi = 4000;
    for (int k = 0; k < 200; k++) {
      const double mid = 0.5 * (lo + hi);
      (ng::arc_time(two.arcs[0], mid) > ng::arc_time(two.arcs[1], 4000 - mid) ? hi : lo) = mid;
    }
    const auto ue = ng::assign(two, {ng::Objective::user, 1e-13});
    s.check("Two roads of the Bureau of Public Roads: the equilibrium where their times are equal (by bisection)", std::fabs(ue.flow[0] - 0.5 * (lo + hi)) / 4000, 1e-9);
  }
  {
    namespace f = congestion_fixtures;
    ng::Network grid{f::nodes, {}, {}};
    for (std::size_t a = 0; a < f::arcs; a++) grid.arcs.push_back(ng::bpr(f::tail[a], f::head[a], f::free_time[a], f::capacity[a]));
    for (std::size_t w = 0; w < f::demands; w++) grid.demands.push_back({f::origin[w], f::destination[w], f::trips[w]});
    const std::vector<double> ref_ue(f::equilibrium, f::equilibrium + f::arcs), ref_so(f::optimum, f::optimum + f::arcs);
    s.check("The grid: scipy's equilibrium has a relative gap of zero by the code here", ng::relative_gap(grid, ref_ue), 1e-11);
    s.check("The grid: scipy's optimum has a relative gap of zero at the marginal costs", ng::relative_gap(grid, ref_so, ng::Objective::system), 1e-11);
    const auto ue = ng::assign(grid, {ng::Objective::user, 1e-5});
    const auto so = ng::assign(grid, {ng::Objective::system, 1e-5});
    s.require("The grid: Frank–Wolfe reaches a relative gap of 10⁻⁵ at the equilibrium and at the optimum", ue.converged && so.converged,
              std::to_string(ue.iterations) + " and " + std::to_string(so.iterations) + " iterations");
    double eu = 0, es = 0, top = 0;
    for (std::size_t a = 0; a < f::arcs; a++) eu = std::max(eu, std::fabs(ue.flow[a] - ref_ue[a])), es = std::max(es, std::fabs(so.flow[a] - ref_so[a])), top = std::max(top, ref_ue[a]);
    s.check("The grid: the equilibrium arc flows against scipy's path-based solution, relative to the largest", eu / top, 5e-3);
    s.check("The grid: the optimum arc flows against scipy's", es / top, 5e-3);
    s.check("The grid: the total time at equilibrium against scipy's", (ue.total_time - f::equilibrium_total) / f::equilibrium_total, 1e-4);
    s.check("The grid: the total time at the optimum against scipy's", (so.total_time - f::optimum_total) / f::optimum_total, 1e-4);
    s.check("The grid: Beckmann's potential at equilibrium against scipy's, within what the bound allows", (ue.objective - f::beckmann) / (ue.objective - ue.lower_bound), 1);
    s.check("The grid: the total time at the optimum against scipy's, within what the bound allows", (so.objective - f::optimum_total) / (so.objective - so.lower_bound), 1);
    s.require("The grid: the lower bound lies below the least potential, and the potential above",
              ue.lower_bound <= f::beckmann * (1 + 1e-12) && ue.objective >= f::beckmann * (1 - 1e-12));
    const auto sk = ng::skims(grid, ng::times(grid, ref_ue), {0}, {8});
    // The least time from 0 to 8 by brute force: every simple path of the grid, by depth-first search.
    const auto t = ng::times(grid, ref_ue);
    double path = 1e300;
    std::vector<char> seen(f::nodes, 0);
    std::function<void(std::size_t, double)> walk = [&](std::size_t v, double sum) {
      if (v == 8) { path = std::min(path, sum); return; }
      seen[v] = 1;
      for (std::size_t a = 0; a < f::arcs; a++)
        if (f::tail[a] == v && !seen[f::head[a]]) walk(f::head[a], sum + t[a]);
      seen[v] = 0;
    };
    walk(0, 0);
    s.check("The grid: the skim from 0 to 8 at equilibrium is the least time of all its simple paths", (sk[0][0] - path) / path, 1e-12);
  }
  {
    calx::random::Stream st(20261111, 0x63676e73u, 0, 0);  // "cgns"
    double worst = 0;
    bool ordered = true;
    for (int g = 0; g < 40; g++) {
      ng::Network net{6, {}, {}};
      for (std::size_t u = 0; u < 6; u++)
        for (std::size_t v = 0; v < 6; v++)
          if (u != v && st.uniform(0, 1) < 0.45) net.arcs.push_back({u, v, st.uniform(0, 5), st.uniform(0, 0.02), 1});
      for (std::size_t v = 0; v + 1 < 6; v++) net.arcs.push_back({v, v + 1, st.uniform(5, 10), st.uniform(0, 0.02), 1});   // a route for every pair forward
      for (int w = 0; w < 3; w++) {
        const auto o = static_cast<std::size_t>(st.uniform(0, 3)), d = 3 + static_cast<std::size_t>(st.uniform(0, 3));
        net.demands.push_back({std::min<std::size_t>(o, 2), std::min<std::size_t>(d, 5), st.uniform(50, 500)});
      }
      const auto ue = ng::assign(net, {ng::Objective::user, 1e-10});
      const auto so = ng::assign(net, {ng::Objective::system, 1e-10});
      ordered = ordered && so.total_time <= ue.total_time * (1 + 1e-9);
      worst = std::max(worst, ue.total_time / so.total_time);
    }
    s.require("Forty random networks with affine times: the optimum's total time at most the equilibrium's", ordered);
    s.require("Forty random networks with affine times: the price of anarchy at most 4/3 (Roughgarden and Tardos)", worst <= 4.0 / 3 + 1e-6, "worst " + std::to_string(worst));
  }
  return s.done();
}
