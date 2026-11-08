// Accessibility: how much of what people need lies within their reach, and who is left without it (nomos's decision 0008; the
// infrastructure's ADR 0086).
//
// The data are places of demand i (people P_i), places of supply j (a capacity S_j: beds, places at school, floor space), and the
// travel costs t_ij of the pairs within reach, which the caller computes on a network (urban science's access) and gives as a sparse
// list. A decay f weighs a cost: binary (f = 1 for t ≤ T), exponential e^{−βt}, or Gaussian within T (Dai 2010); f(t) = 0 beyond T.
//
// - Cumulative opportunities: A_i = Σ_j S_j [t_ij ≤ T], the supply within reach (Wachs and Kumagai 1973).
// - Gravity (Hansen 1959, J. Am. Inst. Planners 25, 73–76): A_i = Σ_j S_j f(t_ij), opportunities weighted by their nearness.
// - The two-step floating catchment area (2SFCA; Luo and Wang 2003, Environ. Plann. B 30, 865–884), with a decay in both steps
//   (E2SFCA; Luo and Qi 2009): each supply's ratio to the demand that can reach it, R_j = S_j / Σ_k f(t_kj) P_k, then each place's
//   access, the ratios within its reach, A_i = Σ_j f(t_ij) R_j. Supply is shared, not counted twice: Σ_i P_i A_i = Σ_j S_j over the
//   supplies some demand reaches (a conservation law, checked), so A_i is supply per person, comparable with Σ S / Σ P.
// - Unmet need: the people with nothing within reach (A_i = 0), and those whose access falls below a share of the region's mean
//   (the "shortage areas" of the method's practice; a declared threshold). The inequality of access across people is the
//   population-weighted Gini of A (inequality.hpp).
// - The Rural Access Index (Roberts, KC and Rastogi 2006, World Bank Transport Paper TP-10; SDG indicator 9.1.1): the share of the rural
//   population living within 2 km of an all-season road.

#pragma once

#include <nomos/measurement/inequality.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::measurement {

struct Reachable {
  std::size_t demand = 0, supply = 0;
  double cost = 0;  // minutes, metres or whatever the caller measures
};

enum class Decay { binary, exponential, gaussian };

struct Catchment {
  double limit = 15;   // T: nothing beyond it is within reach
  Decay decay = Decay::binary;
  double beta = 0.1;   // the exponential's rate (per unit of cost)
};

// f(t) under the catchment.
inline double weight(const Catchment& c, double t) {
  if (!(t <= c.limit) || t < 0) return 0;
  switch (c.decay) {
    case Decay::binary: return 1;
    case Decay::exponential: return std::exp(-c.beta * t);
    case Decay::gaussian: {
      // Dai (2010): f(t) = (e^{−(t/T)²/2} − e^{−1/2}) / (1 − e^{−1/2}), so that f(0) = 1 and f(T) = 0.
      const double g = std::exp(-0.5), x = std::exp(-0.5 * (t / c.limit) * (t / c.limit));
      return (x - g) / (1 - g);
    }
  }
  return 0;
}

namespace detail_access {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::measurement (accessibility): ") + what);
}
inline void check(const std::vector<Reachable>& pairs, std::size_t demands, std::size_t supplies) {
  for (const auto& p : pairs) require(p.demand < demands && p.supply < supplies && p.cost >= 0, "pairs of a demand, a supply and a cost ≥ 0");
}
}  // namespace detail_access

// The supply within reach of each place of demand: Σ_j S_j f(t_ij), with f binary for cumulative opportunities.
inline std::vector<double> gravity(const std::vector<Reachable>& pairs, std::span<const double> supply, std::size_t demands, const Catchment& c) {
  detail_access::check(pairs, demands, supply.size());
  std::vector<double> a(demands, 0.0);
  for (const auto& p : pairs) a[p.demand] += supply[p.supply] * weight(c, p.cost);
  return a;
}

inline std::vector<double> cumulative(const std::vector<Reachable>& pairs, std::span<const double> supply, std::size_t demands, double limit) {
  return gravity(pairs, supply, demands, Catchment{limit, Decay::binary, 0});
}

struct FloatingCatchment {
  std::vector<double> ratio;   // R_j, per supply (NaN where no demand reaches it)
  std::vector<double> access;  // A_i, per place of demand
  double supply_reached = 0;   // Σ S_j over supplies some demand reaches
  double regional = 0;         // Σ S / Σ P over the places and supplies given
};

inline FloatingCatchment two_step(const std::vector<Reachable>& pairs, std::span<const double> demand, std::span<const double> supply, const Catchment& c) {
  detail_access::check(pairs, demand.size(), supply.size());
  FloatingCatchment f;
  std::vector<double> served(supply.size(), 0.0);
  for (const auto& p : pairs) served[p.supply] += weight(c, p.cost) * demand[p.demand];
  f.ratio.assign(supply.size(), std::numeric_limits<double>::quiet_NaN());
  for (std::size_t j = 0; j < supply.size(); j++)
    if (served[j] > 0) f.ratio[j] = supply[j] / served[j], f.supply_reached += supply[j];
  f.access.assign(demand.size(), 0.0);
  for (const auto& p : pairs)
    if (served[p.supply] > 0) f.access[p.demand] += weight(c, p.cost) * f.ratio[p.supply];
  double S = 0, P = 0;
  for (double s : supply) S += s;
  for (double d : demand) P += d;
  f.regional = P > 0 ? S / P : 0;
  return f;
}

struct UnmetNeed {
  double people = 0, without = 0, short_of = 0;  // all, with nothing in reach, below the share of the mean
  double share_without = 0, share_short = 0;
  double mean_access = 0;                        // population-weighted
  double gini = 0;                               // of access across people
};

// Who is left without: nothing within reach, or access below `share` of the population-weighted mean.
inline UnmetNeed unmet_need(std::span<const double> demand, std::span<const double> access, double share = 0.5) {
  detail_access::require(demand.size() == access.size() && !demand.empty(), "an access for each place of demand");
  UnmetNeed u;
  double weighted = 0;
  for (std::size_t i = 0; i < demand.size(); i++) u.people += demand[i], weighted += demand[i] * access[i];
  detail_access::require(u.people > 0, "people");
  u.mean_access = weighted / u.people;
  for (std::size_t i = 0; i < demand.size(); i++) {
    if (access[i] <= 0) u.without += demand[i];
    if (access[i] < share * u.mean_access) u.short_of += demand[i];
  }
  u.share_without = u.without / u.people, u.share_short = u.short_of / u.people;
  u.gini = u.mean_access > 0 ? gini(access, demand) : 0;
  return u;
}

// The Rural Access Index: the share of the rural population within `within` (2 km) of an all-season road, from each place's people,
// whether it is rural, and its distance to the nearest such road.
struct RuralAccess {
  double rural = 0, served = 0, index = std::numeric_limits<double>::quiet_NaN();
};
inline RuralAccess rural_access(std::span<const double> people, std::span<const int> rural, std::span<const double> distance, double within = 2000) {
  detail_access::require(people.size() == rural.size() && people.size() == distance.size(), "people, rurality and a distance for each place");
  RuralAccess r;
  for (std::size_t i = 0; i < people.size(); i++)
    if (rural[i]) r.rural += people[i], r.served += distance[i] <= within ? people[i] : 0.0;
  if (r.rural > 0) r.index = r.served / r.rural;
  return r;
}

}  // namespace nomos::measurement
