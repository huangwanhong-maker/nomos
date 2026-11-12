// Coalitional games and their Shapley value; and the division of an estate that is short of the claims on it.
//
// A coalitional game on n players gives a worth v(S) to every set S of them, v(∅) = 0. The Shapley
// value (Shapley 1953) gives player i the mean, over the n! orders in which the players may
// arrive, of what i adds to those who came before:
//
//   φ_i = Σ_{S ∌ i} |S|! (n − |S| − 1)! / n! · [v(S ∪ {i}) − v(S)].
//
// It is the one value that is efficient (Σ φ = v(N)), symmetric, zero for a player who adds
// nothing, and additive in the game. Computed here exactly over the 2ⁿ sets, n ≤ 24.
//
// A claims problem (O'Neill 1982; the bankruptcy problem of Aumann and Maschler 1985) is an estate
// E ≥ 0 and claims c_i ≥ 0 with Σ c ≥ E. A rule gives awards x with 0 ≤ x_i ≤ c_i and Σ x = E:
//
//   priority                  the claimants are paid in a given order, each in full while the estate lasts
//   random_arrival            the mean of the priority awards over all orders (O'Neill): the Shapley value of the
//                             game w(S) = max(0, E − Σ_{j ∉ S} c_j), what S can secure after the others are paid in full
//   proportional              x_i = E c_i / Σ c
//   constrained_equal_awards  x_i = min(c_i, λ), λ set so that Σ x = E
//   constrained_equal_losses  x_i = max(0, c_i − λ)
//   talmud                    Aumann and Maschler's rule, the one consistent with the contested-garment division of every
//                             pair: equal awards on the half-claims while E ≤ Σ c / 2, equal losses on them beyond
//   interleaved               the claimants paid by turns in small steps: at each turn claimant i receives
//                             min(h r_i, what is left of the claim, what is left of the estate). As h → 0 the awards
//                             no longer depend on the order of the turns, and tend to the rule of the rates r:
//                             r ∝ c gives the proportional rule, equal r the constrained equal awards.
//
// With two claimants random arrival, the Talmud rule and the contested-garment division coincide.
//
// References. Shapley (1953), Contributions to the Theory of Games II: 307–317. O'Neill (1982),
// Math. Social Sci. 2: 345–371. Aumann and Maschler (1985), J. Econ. Theory 36: 195–213.
// Thomson (2003), Math. Social Sci. 45: 249–297 (the survey; the rules' names follow it).

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nomos::games {

// The Shapley value of the game v on n players; v is asked once for each of the 2ⁿ sets, as a mask of bits.
inline std::vector<double> shapley(std::size_t n, const std::function<double(std::uint32_t)>& v) {
  if (n == 0 || n > 24) throw std::invalid_argument("nomos::games::shapley: between 1 and 24 players");
  const std::uint32_t sets = std::uint32_t{1} << n;
  std::vector<double> worth(sets);
  for (std::uint32_t S = 0; S < sets; S++) worth[S] = v(S);
  // weight[s] = s! (n − s − 1)! / n!, by the recurrence weight[s + 1] = weight[s] (s + 1)/(n − s − 1).
  std::vector<double> weight(n);
  weight[0] = 1.0 / static_cast<double>(n);
  for (std::size_t k = 0; k + 1 < n; k++) weight[k + 1] = weight[k] * static_cast<double>(k + 1) / static_cast<double>(n - k - 1);
  std::vector<double> phi(n, 0.0);
  for (std::uint32_t S = 0; S < sets; S++) {
    const auto size = static_cast<std::size_t>(__builtin_popcount(S));
    for (std::size_t i = 0; i < n; i++) {
      const std::uint32_t bit = std::uint32_t{1} << i;
      if (S & bit) continue;
      phi[i] += weight[size] * (worth[S | bit] - worth[S]);
    }
  }
  return phi;
}

struct Claims {
  double estate = 0;
  std::vector<double> claims;
};

inline void validate(const Claims& P) {
  if (!(P.estate >= 0) || P.claims.empty()) throw std::invalid_argument("nomos::games::Claims: an estate that is not negative, and claims on it");
  double total = 0;
  for (double c : P.claims) {
    if (!(c >= 0)) throw std::invalid_argument("nomos::games::Claims: claims are not negative");
    total += c;
  }
  if (total < P.estate * (1 - 1e-12)) throw std::invalid_argument("nomos::games::Claims: the estate exceeds the claims on it: there is nothing to divide by a rule");
}

// The claimants paid in the order given (a permutation of their indices).
inline std::vector<double> priority(const Claims& P, const std::vector<std::size_t>& order) {
  validate(P);
  if (order.size() != P.claims.size()) throw std::invalid_argument("nomos::games::priority: an order of all the claimants");
  std::vector<double> x(P.claims.size(), 0.0);
  double left = P.estate;
  for (std::size_t i : order) {
    x.at(i) = std::min(P.claims[i], left);
    left -= x[i];
  }
  return x;
}

// What a set of claimants can secure after the others are paid in full.
inline double secured(const Claims& P, std::uint32_t S) {
  double others = 0;
  for (std::size_t j = 0; j < P.claims.size(); j++)
    if (!(S & (std::uint32_t{1} << j))) others += P.claims[j];
  return std::max(0.0, P.estate - others);
}

inline std::vector<double> random_arrival(const Claims& P) {
  validate(P);
  return shapley(P.claims.size(), [&](std::uint32_t S) { return secured(P, S); });
}

inline std::vector<double> proportional(const Claims& P) {
  validate(P);
  const double total = std::accumulate(P.claims.begin(), P.claims.end(), 0.0);
  std::vector<double> x(P.claims.size(), 0.0);
  if (total > 0)
    for (std::size_t i = 0; i < x.size(); i++) x[i] = P.estate * P.claims[i] / total;
  return x;
}

// x_i = min(c_i, λ) with Σ x = E: the claims in increasing order are paid in full while an equal share of what is left covers them.
inline std::vector<double> constrained_equal_awards(const Claims& P) {
  validate(P);
  const std::size_t n = P.claims.size();
  std::vector<std::size_t> order(n);
  std::iota(order.begin(), order.end(), std::size_t{0});
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return P.claims[a] < P.claims[b]; });
  std::vector<double> x(n, 0.0);
  double left = P.estate;
  for (std::size_t k = 0; k < n; k++) {
    const double share = left / static_cast<double>(n - k);
    x[order[k]] = std::min(P.claims[order[k]], share);
    left -= x[order[k]];
  }
  return x;
}

// x_i = max(0, c_i − λ): the equal awards of the losses Σ c − E, taken from the claims.
inline std::vector<double> constrained_equal_losses(const Claims& P) {
  validate(P);
  const double total = std::accumulate(P.claims.begin(), P.claims.end(), 0.0);
  const auto loss = constrained_equal_awards({total - P.estate, P.claims});
  std::vector<double> x(P.claims.size());
  for (std::size_t i = 0; i < x.size(); i++) x[i] = P.claims[i] - loss[i];
  return x;
}

inline std::vector<double> talmud(const Claims& P) {
  validate(P);
  const double total = std::accumulate(P.claims.begin(), P.claims.end(), 0.0);
  std::vector<double> half(P.claims);
  for (double& c : half) c /= 2;
  if (P.estate <= total / 2) return constrained_equal_awards({P.estate, half});
  const auto loss = constrained_equal_awards({total - P.estate, half});
  std::vector<double> x(P.claims.size());
  for (std::size_t i = 0; i < x.size(); i++) x[i] = P.claims[i] - loss[i];
  return x;
}

// The division of a garment claimed by two (Baba Metzia 2a): each concedes to the other what it does not claim, and the
// rest is halved.
inline std::vector<double> contested_garment(double estate, double c1, double c2) {
  const double to1 = std::max(0.0, estate - c2), to2 = std::max(0.0, estate - c1), rest = estate - to1 - to2;
  return {to1 + rest / 2, to2 + rest / 2};
}

// The claimants paid by turns, in the order given and round again, each turn min(step · rate, the claim left, the estate left).
inline std::vector<double> interleaved(const Claims& P, const std::vector<std::size_t>& order, const std::vector<double>& rates, double step) {
  validate(P);
  const std::size_t n = P.claims.size();
  if (order.size() != n || rates.size() != n || !(step > 0)) throw std::invalid_argument("nomos::games::interleaved: an order and a rate for each claimant, and a positive step");
  double moving = 0;
  for (std::size_t i = 0; i < n; i++) {
    if (!(rates[i] >= 0)) throw std::invalid_argument("nomos::games::interleaved: rates are not negative");
    if (P.claims[i] > 0) moving += rates[i];
  }
  if (!(moving > 0) && P.estate > 0) throw std::invalid_argument("nomos::games::interleaved: some claimant with a claim is paid at a positive rate");
  std::vector<double> x(n, 0.0);
  double left = P.estate;
  while (left > 0) {
    double paid = 0;
    for (std::size_t i : order) {
      const double pay = std::min({step * rates.at(i), P.claims[i] - x[i], left});
      x[i] += pay;
      left -= pay;
      paid += pay;
    }
    if (!(paid > 0)) break;  // every claimant with a positive rate is paid in full
  }
  return x;
}

// Over every order of the claimants: the least and the greatest each may be awarded by priority. The least is what is left
// when all the others are paid before, the greatest what is taken when paid first.
struct OrderRange {
  std::vector<double> least, greatest;
};
inline OrderRange priority_range(const Claims& P) {
  validate(P);
  const std::size_t n = P.claims.size();
  const double total = std::accumulate(P.claims.begin(), P.claims.end(), 0.0);
  OrderRange r{std::vector<double>(n), std::vector<double>(n)};
  for (std::size_t i = 0; i < n; i++) {
    r.least[i] = std::min(P.claims[i], std::max(0.0, P.estate - (total - P.claims[i])));
    r.greatest[i] = std::min(P.claims[i], P.estate);
  }
  return r;
}

}  // namespace nomos::games
