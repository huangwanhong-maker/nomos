// Checks of accessibility (accessibility.hpp).
//
// - One clinic of 100 places shared by two neighbourhoods of 50 and 150 people: its ratio is 0.5 a person, and each has access 0.5.
// - The two-step floating catchment conserves supply: Σ_i P_i A_i equals the supply some demand reaches, on random places with the
//   binary, exponential and Gaussian decays alike.
// - Cumulative opportunities count what lies within the limit; gravity with a vanishing rate is cumulative opportunities.
// - The decays: f(0) = 1, and the Gaussian's f(T) = 0.
// - Unmet need: the people with nothing in reach, the people below half the mean, and a Gini of zero when every place has the same
//   access.
// - The Rural Access Index from people, rurality and distances.

#include <nomos/measurement/accessibility.hpp>

#include <calx/random/stream.hpp>
#include <calx/verification/check.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace nm = nomos::measurement;

int main() {
  calx::verification::Suite s;
  {
    const std::vector<nm::Reachable> pairs{{0, 0, 5}, {1, 0, 9}};
    const auto f = nm::two_step(pairs, std::vector<double>{50, 150}, std::vector<double>{100}, nm::Catchment{15});
    s.check("2SFCA: one clinic of 100 shared by 50 and 150 people: a ratio of 0.5 a person, and access 0.5 for both",
            std::fabs(f.ratio[0] - 0.5) + std::fabs(f.access[0] - 0.5) + std::fabs(f.access[1] - 0.5), 1e-15);
  }
  {
    calx::random::Stream st(20261108, 0x61636365u, 0, 0);  // "acce"
    const std::size_t D = 60, S = 12;
    std::vector<double> P(D), C(S), dx(D), dy(D), sx(S), sy(S);
    for (std::size_t i = 0; i < D; i++) P[i] = st.uniform(10, 500), dx[i] = st.uniform(0, 10), dy[i] = st.uniform(0, 10);
    for (std::size_t j = 0; j < S; j++) C[j] = st.uniform(5, 50), sx[j] = st.uniform(0, 10), sy[j] = st.uniform(0, 10);
    std::vector<nm::Reachable> pairs;
    for (std::size_t i = 0; i < D; i++)
      for (std::size_t j = 0; j < S; j++) {
        const double t = 4 * std::hypot(dx[i] - sx[j], dy[i] - sy[j]);  // minutes, at 15 km/h over km
        if (t <= 20) pairs.push_back({i, j, t});
      }
    double worst = 0;
    for (const auto decay : {nm::Decay::binary, nm::Decay::exponential, nm::Decay::gaussian}) {
      const auto f = nm::two_step(pairs, P, C, nm::Catchment{20, decay, 0.15});
      double total = 0;
      for (std::size_t i = 0; i < D; i++) total += P[i] * f.access[i];
      worst = std::max(worst, std::fabs(total - f.supply_reached) / f.supply_reached);
    }
    s.check("2SFCA conserves supply: Σ P A is the supply reached, with binary, exponential and Gaussian decays (relative)", worst, 1e-12);
    const auto cum = nm::cumulative(pairs, std::vector<double>(S, 1.0), D, 10);
    double miss = 0;
    for (std::size_t i = 0; i < D; i++) {
      double count = 0;
      for (const auto& p : pairs)
        if (p.demand == i && p.cost <= 10) count += 1;
      miss = std::max(miss, std::fabs(cum[i] - count));
    }
    s.check("cumulative opportunities count the supplies within the limit", miss, 0);
    const auto g = nm::gravity(pairs, C, D, nm::Catchment{20, nm::Decay::exponential, 1e-12});
    const auto c20 = nm::cumulative(pairs, C, D, 20);
    double gap = 0;
    for (std::size_t i = 0; i < D; i++) gap = std::max(gap, std::fabs(g[i] - c20[i]) / std::max(1.0, c20[i]));
    s.check("gravity with a vanishing rate is cumulative opportunities", gap, 1e-9);
  }
  {
    const nm::Catchment gauss{30, nm::Decay::gaussian, 0};
    s.check("the decays: f(0) = 1 for each, and the Gaussian's f(T) = 0",
            std::fabs(nm::weight(gauss, 0) - 1) + std::fabs(nm::weight(gauss, 30)) + std::fabs(nm::weight(nm::Catchment{30, nm::Decay::exponential, 0.2}, 0) - 1), 1e-15);
  }
  {
    const std::vector<double> people{100, 200, 300, 400}, access{0, 0.2, 1.0, 1.2};
    const auto u = nm::unmet_need(people, access, 0.5);
    // Mean access (0·100 + 0.2·200 + 1·300 + 1.2·400)/1000 = 0.82; below 0.41: the first two places.
    s.check("unmet need: 100 of 1,000 with nothing in reach, 300 below half the mean (0.41)", std::fabs(u.share_without - 0.1) + std::fabs(u.share_short - 0.3) + std::fabs(u.mean_access - 0.82), 1e-15);
    const auto even = nm::unmet_need(people, std::vector<double>{0.7, 0.7, 0.7, 0.7});
    s.check("... and a Gini of zero when every place has the same access", even.gini, 1e-15);
  }
  {
    const auto r = nm::rural_access(std::vector<double>{100, 300, 600, 1000}, std::vector<int>{1, 1, 1, 0}, std::vector<double>{500, 2500, 1800, 9000});
    s.check("the Rural Access Index: 700 of 1,000 rural people within 2 km of a road", r.index - 0.7, 1e-15);
  }
  return s.done();
}
