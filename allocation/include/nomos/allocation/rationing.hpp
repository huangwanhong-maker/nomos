// How much each receives when there is not enough: the rationing rules beyond those of nomos::games's claims problems,
// and the formulas by which public bodies divide a sum among places and groups.
//
// A claims problem is an estate E ≥ 0 and claims c_i ≥ 0 with Σc ≥ E (nomos::games::Claims; its rules there — priority,
// random arrival, proportional, constrained equal awards and losses, the Talmud rule — are used, not written again).
// Added here, in Thomson's names (Thomson 2003, Math. Social Sci. 45, 249–297; 2015, Math. Social Sci. 74, 41–59):
//
//   piniles               equal awards on the half-claims while E ≤ C/2, then the half-claims plus equal awards on the
//                         other halves (Piniles 1861, as Thomson reports him): the Talmud rule's first half twice over
//   adjusted_proportional each claimant first receives what the others leave, m_i = max(0, E − Σ_{j≠i} c_j); the claims
//                         are then truncated at what remains and the remainder divided in proportion (Curiel, Maschler
//                         and Tijs 1987, Z. Oper. Res. 31, A143–A159)
//   tal(θ)                the family between the losses and the awards: equal awards on θc while E ≤ θC, then θc plus
//                         equal losses on (1 − θ)c; θ = 0 is constrained equal losses, θ = ½ the Talmud rule, θ = 1
//                         constrained equal awards (Moreno-Ternero and Villar 2006, Soc. Choice Welf. 27, 231–257)
//   dual(R)               R^d(c, E) = c − R(c, C − E): what R does to awards, R^d does to losses; proportional and Talmud
//                         are their own duals, equal awards and equal losses each other's
//   weighted_proportional awards in proportion to declared weights, each capped at its claim, the excess passed on
//   priority_classes      the claimants in ordered classes, a class served in full before the next, the class where the
//                         estate runs out divided by a declared rule (Moulin 2000, Econometrica 68, 643–684: with
//                         consistency and composition, rules of this form are the only ones)
//
// And the formulas by which a sum is divided among units in practice, each a declared rule whose parts are reported apart:
//
//   sufficiency_first     the units' shortfalls below a declared threshold are closed first (by a declared rationing
//                         rule if the sum cannot close them all), and what remains is divided by declared weights. It
//                         is also the acequias' rule of the survey (Part IV §3): shares by entitlement in plenty, survival
//                         shares first in scarcity, the switch made by the size of the sum, not by anyone's choice.
//   banded_points         each unit graded on indicators by declared bands, its points the weighted sum of its grades, and
//                         the sum shared in proportion to points: the procedure of Porto Alegre's participatory budget
//                         (Santos 1998, Politics & Society 26, 461–510, pp. 474–475 as the survey read them: need,
//                         population and the community's own priority, weighted 4, 2 and 4 in 1998)
//   equalisation          a grant to each unit of λ P_i [(e_i − ē) + Σ_k t_k (b̄_k − b_ik)]: what it needs to spend per head
//                         above the standard, plus what it cannot raise per head below the standard at standard rates
//                         (the representative expenditure-need and tax system of horizontal fiscal equalisation; the
//                         survey, Part III §2), then a declared envelope of floors, caps and phasing whose effect on each
//                         unit is reported separately from the formula's.
#pragma once

#include <nomos/games/claims.hpp>

#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::allocation {

using games::Claims;
using Rule = std::function<std::vector<double>(const Claims&)>;

namespace detail {

inline double total(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0); }

inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("nomos::allocation: " + what);
}

// Equal awards on caps a with total A ≥ E: x_i = min(a_i, λ).
inline std::vector<double> equal_awards(const std::vector<double>& a, double E) {
  const double A = total(a);
  std::vector<double> x(a.size(), 0);
  if (E >= A) return a;
  double lo = 0, hi = *std::max_element(a.begin(), a.end());
  for (int k = 0; k < 200; k++) {
    const double mid = 0.5 * (lo + hi);
    double s = 0;
    for (double v : a) s += std::min(v, mid);
    (s < E ? lo : hi) = mid;
  }
  for (std::size_t i = 0; i < a.size(); i++) x[i] = std::min(a[i], hi);
  // Absorb the bisection's rounding into the uncapped claimants, so that Σx = E exactly.
  const double gap = E - total(x);
  std::size_t free = 0;
  for (std::size_t i = 0; i < a.size(); i++) free += x[i] < a[i];
  if (free > 0)
    for (std::size_t i = 0; i < a.size(); i++)
      if (x[i] < a[i]) x[i] += gap / static_cast<double>(free);
  return x;
}

inline std::vector<double> equal_losses(const std::vector<double>& a, double E) {
  // x_i = max(0, a_i − λ): equal awards on the losses.
  const double A = total(a);
  const auto loss = equal_awards(a, A - E);
  std::vector<double> x(a.size());
  for (std::size_t i = 0; i < a.size(); i++) x[i] = a[i] - loss[i];
  return x;
}

}  // namespace detail

// ---- rationing rules --------------------------------------------------------------------------------------------

inline std::vector<double> piniles(const Claims& P) {
  games::validate(P);
  const double C = detail::total(P.claims);
  std::vector<double> half(P.claims.size());
  for (std::size_t i = 0; i < half.size(); i++) half[i] = P.claims[i] / 2;
  if (P.estate <= C / 2) return detail::equal_awards(half, P.estate);
  const auto rest = detail::equal_awards(half, P.estate - C / 2);
  std::vector<double> x(half.size());
  for (std::size_t i = 0; i < x.size(); i++) x[i] = half[i] + rest[i];
  return x;
}

inline std::vector<double> adjusted_proportional(const Claims& P) {
  games::validate(P);
  const double C = detail::total(P.claims);
  const std::size_t n = P.claims.size();
  std::vector<double> m(n), revised(n);
  for (std::size_t i = 0; i < n; i++) m[i] = std::max(0.0, P.estate - (C - P.claims[i]));
  const double left = P.estate - detail::total(m);
  for (std::size_t i = 0; i < n; i++) revised[i] = std::min(P.claims[i] - m[i], left);
  const double R = detail::total(revised);
  std::vector<double> x(n);
  for (std::size_t i = 0; i < n; i++) x[i] = m[i] + (R > 0 ? left * revised[i] / R : 0.0);
  return x;
}

inline std::vector<double> tal(const Claims& P, double theta) {
  games::validate(P);
  detail::require(theta >= 0 && theta <= 1, "the TAL family's θ lies in [0, 1]");
  const double C = detail::total(P.claims);
  std::vector<double> first(P.claims.size()), second(P.claims.size());
  for (std::size_t i = 0; i < first.size(); i++) first[i] = theta * P.claims[i], second[i] = (1 - theta) * P.claims[i];
  if (P.estate <= theta * C) return detail::equal_awards(first, P.estate);
  const auto rest = detail::equal_losses(second, P.estate - theta * C);
  std::vector<double> x(first.size());
  for (std::size_t i = 0; i < x.size(); i++) x[i] = first[i] + rest[i];
  return x;
}

// The dual of a rule: R^d(c, E) = c − R(c, C − E).
inline Rule dual(Rule R) {
  return [R](const Claims& P) {
    games::validate(P);
    const double C = detail::total(P.claims);
    const auto y = R(Claims{C - P.estate, P.claims});
    std::vector<double> x(P.claims.size());
    for (std::size_t i = 0; i < x.size(); i++) x[i] = P.claims[i] - y[i];
    return x;
  };
}

inline std::vector<double> weighted_proportional(const Claims& P, const std::vector<double>& weights) {
  games::validate(P);
  const std::size_t n = P.claims.size();
  detail::require(weights.size() == n, "a weight for every claimant");
  for (double w : weights) detail::require(w >= 0, "weights are not negative");
  std::vector<double> x(n, 0);
  std::vector<bool> capped(n, false);
  double left = P.estate;
  // Water-filling: share what is left among those not yet at their claims, in proportion to their weights.
  for (std::size_t round = 0; round <= n && left > 1e-15 * std::max(1.0, P.estate); round++) {
    double W = 0;
    for (std::size_t i = 0; i < n; i++)
      if (!capped[i]) W += weights[i];
    if (!(W > 0)) break;
    double scale = left / W;
    for (std::size_t i = 0; i < n; i++)
      if (!capped[i] && weights[i] > 0) scale = std::min(scale, (P.claims[i] - x[i]) / weights[i]);
    for (std::size_t i = 0; i < n; i++)
      if (!capped[i]) x[i] += scale * weights[i];
    left = P.estate - detail::total(x);
    for (std::size_t i = 0; i < n; i++)
      if (!capped[i] && x[i] >= P.claims[i] * (1 - 1e-12)) capped[i] = true, x[i] = P.claims[i];
  }
  return x;
}

// The claimants in classes (class[i] = 0 first), each class served in full before the next; the class where the estate runs
// out divided by `within`.
inline std::vector<double> priority_classes(const Claims& P, const std::vector<std::size_t>& cls, const Rule& within) {
  games::validate(P);
  const std::size_t n = P.claims.size();
  detail::require(cls.size() == n, "a class for every claimant");
  std::vector<double> x(n, 0);
  double left = P.estate;
  const std::size_t last = *std::max_element(cls.begin(), cls.end());
  for (std::size_t k = 0; k <= last && left > 0; k++) {
    std::vector<std::size_t> members;
    double demand = 0;
    for (std::size_t i = 0; i < n; i++)
      if (cls[i] == k) members.push_back(i), demand += P.claims[i];
    if (members.empty()) continue;
    if (demand <= left) {
      for (std::size_t i : members) x[i] = P.claims[i];
      left -= demand;
    } else {
      Claims sub{left, {}};
      for (std::size_t i : members) sub.claims.push_back(P.claims[i]);
      const auto y = within(sub);
      for (std::size_t j = 0; j < members.size(); j++) x[members[j]] = y[j];
      left = 0;
    }
  }
  return x;
}

// ---- formulas of practice ---------------------------------------------------------------------------------------

// The parts of an award, reported apart so that each can be contested on its own ground.
struct Parts {
  std::vector<double> total;
  std::vector<double> first;     // sufficiency_first: towards the threshold; equalisation: the formula
  std::vector<double> second;    // sufficiency_first: the remainder by weights; equalisation: the envelope's adjustment
  bool scarce = false;           // sufficiency_first: the sum could not close every shortfall
  std::string said;
};

// Close each unit's shortfall (`gaps`, what it lacks to a declared threshold) first; if the sum cannot close them all,
// divide it among the shortfalls by `rationing` (applied to the gaps as claims). What remains is divided in proportion to
// `weights`.
inline Parts sufficiency_first(double sum, const std::vector<double>& gaps, const std::vector<double>& weights, const Rule& rationing) {
  const std::size_t n = gaps.size();
  detail::require(weights.size() == n, "a weight for every unit");
  detail::require(sum >= 0, "a sum that is not negative");
  for (double g : gaps) detail::require(g >= 0, "shortfalls are not negative");
  Parts p;
  p.first.assign(n, 0), p.second.assign(n, 0), p.total.assign(n, 0);
  const double G = detail::total(gaps);
  if (sum <= G) {
    p.scarce = true;
    if (sum > 0) p.first = rationing(Claims{sum, gaps});
    p.said = "the sum (" + std::to_string(sum) + ") does not close the shortfalls (" + std::to_string(G) +
             "): it is divided among them by the declared rationing rule, and nothing is left for the second rule";
  } else {
    p.first = gaps;
    const double W = detail::total(weights), left = sum - G;
    detail::require(W > 0, "weights for the remainder that are not all zero");
    for (std::size_t i = 0; i < n; i++) p.second[i] = left * weights[i] / W;
    p.said = "every shortfall is closed; the remaining " + std::to_string(left) + " is divided by the declared weights";
  }
  for (std::size_t i = 0; i < n; i++) p.total[i] = p.first[i] + p.second[i];
  return p;
}

// A grade from declared bands: the number of band edges a value reaches (edges ascending; `rising` false when a lower value
// earns a higher grade, as for income).
inline double band_grade(double value, const std::vector<double>& edges, bool rising = true) {
  std::size_t g = 0;
  for (double e : edges) g += rising ? value >= e : value <= e;
  return static_cast<double>(g + 1);
}

struct Indicator {
  std::string name;
  std::vector<double> values;   // one per unit
  std::vector<double> edges;    // band edges, ascending
  bool rising = true;           // a higher value earns a higher grade
  double weight = 1;
};

struct Points {
  std::vector<std::vector<double>> grades;  // [indicator][unit]
  std::vector<double> points, share, award;
};

inline Points banded_points(double sum, const std::vector<Indicator>& indicators) {
  detail::require(!indicators.empty(), "at least one indicator");
  const std::size_t n = indicators.front().values.size();
  Points p;
  p.points.assign(n, 0);
  for (const auto& ind : indicators) {
    detail::require(ind.values.size() == n, "every indicator a value for every unit");
    std::vector<double> g(n);
    for (std::size_t i = 0; i < n; i++) {
      g[i] = band_grade(ind.values[i], ind.edges, ind.rising);
      p.points[i] += ind.weight * g[i];
    }
    p.grades.push_back(g);
  }
  const double T = detail::total(p.points);
  p.share.resize(n), p.award.resize(n);
  for (std::size_t i = 0; i < n; i++) p.share[i] = T > 0 ? p.points[i] / T : 0, p.award[i] = sum * p.share[i];
  return p;
}

struct Unit {
  std::string name;
  double population = 0;
  double need = 0;                  // e_i, assessed spending need per head
  std::vector<double> bases;        // b_ik, tax bases per head
  double previous = 0;              // last period's grant, for floors, caps and phasing
};

struct Envelope {
  double floor = -1e300;            // the least change from the previous grant, per head (e.g. 0: no unit loses)
  double cap = 1e300;               // the most change per head
  double phase = 1;                 // the share of the change taken this period
};

// λ P_i [(e_i − ē) + Σ_k t_k (b̄_k − b_ik)] with ē and b̄ the population-weighted means, then the envelope: the change from
// the previous grant phased in, held between the floor and the cap. The formula's grant and the envelope's adjustment are
// reported apart.
inline Parts equalisation(const std::vector<Unit>& units, const std::vector<double>& rates, double lambda, const Envelope& env = {}) {
  const std::size_t n = units.size();
  detail::require(n > 0, "at least one unit");
  double P = 0, ebar = 0;
  std::vector<double> bbar(rates.size(), 0);
  for (const auto& u : units) {
    detail::require(u.population > 0, "a population above zero for every unit");
    detail::require(u.bases.size() == rates.size(), "a tax base per head for every rate");
    P += u.population, ebar += u.population * u.need;
    for (std::size_t k = 0; k < rates.size(); k++) bbar[k] += u.population * u.bases[k];
  }
  ebar /= P;
  for (double& b : bbar) b /= P;
  Parts p;
  p.first.resize(n), p.second.resize(n), p.total.resize(n);
  for (std::size_t i = 0; i < n; i++) {
    const auto& u = units[i];
    double capacity = 0;
    for (std::size_t k = 0; k < rates.size(); k++) capacity += rates[k] * (bbar[k] - u.bases[k]);
    p.first[i] = lambda * u.population * ((u.need - ebar) + capacity);
    const double change = p.first[i] - u.previous;
    const double per_head = std::clamp(env.phase * change / u.population, env.floor, env.cap);
    p.total[i] = u.previous + per_head * u.population;
    p.second[i] = p.total[i] - p.first[i];
  }
  p.said = "the formula's grant, and the envelope's adjustment (phasing " + std::to_string(env.phase) + ", floors and caps per head) apart";
  return p;
}

}  // namespace nomos::allocation
