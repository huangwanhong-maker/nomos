// Simple games and a priori voting power.
//
// A simple game on n players says which coalitions win (Shapley 1962; von Neumann and Morgenstern 1944, ch. X). It is taken to be
// monotone (a superset of a winning coalition wins), with the grand coalition winning and the empty one losing. A decision rule
// of a body (a weighted vote with a quota, a double majority, vetoes, a rule deemed met when too few would block it) is a simple
// game, and its a priori power is what the rule alone gives each member, before any preference is known: a measure of the rule,
// never of influence.
//
// A swing of player i is a coalition S ∌ i that loses while S ∪ {i} wins. With η_i the number of i's swings:
//
//   Penrose measure     β′_i = η_i / 2^(n−1): the probability that i is decisive when every other member votes yes or no with equal
//                       chance and independently (Penrose 1946; it is Banzhaf's absolute index, Banzhaf 1965; Dubey and Shapley 1979)
//   Banzhaf index       β_i = η_i / Σ_j η_j, the normalised form, which sums to 1
//   Shapley–Shubik      φ_i, the Shapley value of the game v(S) = 1 if S wins: the share of the n! orders of arrival in which i turns a
//                       losing coalition into a winning one (Shapley and Shubik 1954)
//   Coleman             the power to prevent action, η_i / W, and to initiate it, η_i / (2^n − W), with W the number of winning
//                       coalitions; and the collectivity's power to act, W / 2^n (Coleman 1971)
//
// A veto player belongs to every winning coalition; a dummy has no swing. A game is proper when no two disjoint coalitions both win,
// strong when of every coalition and its complement one wins, decisive when both (Felsenthal and Machover 1998, ch. 2).
//
// Two routes. Any simple game on n ≤ 24 players is computed exactly over its 2ⁿ coalitions. A weighted game (S wins when its weight
// reaches the quota, integer weights) is computed by generating functions on any number of players, in time O(n² q) for the
// Penrose and Banzhaf measures and O(n³ q) for Shapley–Shubik, with q the quota: pseudo-polynomial, so that large weights are the
// caller's to round (Brams and Affuso 1976; Mann and Shapley 1962). The two routes are held to each other by the checks.
//
// The square-root rule. In a two-tier system, the citizens of unit i (m_i of them) choose by simple majority how their delegate
// votes in a council. A citizen is decisive when both their unit's vote and their delegate's vote turn on them; with independent
// equiprobable votes the probability is the product of decisive_in_majority(m_i) ≈ √(2/(π m_i)) and the delegate's Penrose
// measure. Penrose (1946) observed that weights in the council proportional to √m_i make the citizens' indirect power nearly equal
// across units; Słomczyński and Życzkowski (2006) choose the quota that makes it nearly exact.
//
// References. Penrose, L. S. (1946), The elementary statistics of majority voting, J. Royal Statistical Soc. 109, 53–57. Shapley,
// L. S. and Shubik, M. (1954), A method for evaluating the distribution of power in a committee system, APSR 48, 787–792. Banzhaf,
// J. F. (1965), Weighted voting doesn't work: a mathematical analysis, Rutgers Law Review 19, 317–343. Coleman, J. S. (1971),
// Control of collectivities and the power of a collectivity to act, in Lieberman, ed., Social Choice. Dubey, P. and Shapley, L. S.
// (1979), Mathematical properties of the Banzhaf power index, Math. Oper. Res. 4, 99–131. Felsenthal, D. S. and Machover, M.
// (1998), The Measurement of Voting Power, Edward Elgar. Słomczyński, W. and Życzkowski, K. (2006), Penrose voting system and optimal
// quota, Acta Physica Polonica B 37, 3133–3143.

#pragma once

#include <calx/special_functions/gamma.hpp>

#include "claims.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace nomos::games {

struct SimpleGame {
  std::size_t players = 0;
  std::function<bool(std::uint32_t)> wins;  // a coalition as a mask of bits
};

struct WeightedGame {
  std::vector<long> weights;
  long quota = 0;  // S wins when Σ_{i ∈ S} w_i ≥ quota
};

inline void validate(const WeightedGame& g) {
  if (g.weights.empty()) throw std::invalid_argument("nomos::games::WeightedGame: at least one player");
  long total = 0;
  for (long w : g.weights) {
    if (w < 0) throw std::invalid_argument("nomos::games::WeightedGame: weights are not negative");
    total += w;
  }
  if (g.quota <= 0 || g.quota > total) throw std::invalid_argument("nomos::games::WeightedGame: a quota above 0 and at most the total weight");
}

inline SimpleGame simple(const WeightedGame& g) {
  validate(g);
  if (g.weights.size() > 24) throw std::invalid_argument("nomos::games::simple: a weighted game of more than 24 players is computed by generating functions");
  return {g.weights.size(), [g](std::uint32_t S) {
            long w = 0;
            for (std::size_t i = 0; i < g.weights.size(); i++)
              if (S >> i & 1u) w += g.weights[i];
            return w >= g.quota;
          }};
}

// A compound game: wins when it wins every component (a double majority), or any component.
inline SimpleGame all_of(const std::vector<SimpleGame>& games) {
  if (games.empty()) throw std::invalid_argument("nomos::games::all_of: at least one game");
  for (const auto& g : games)
    if (g.players != games[0].players) throw std::invalid_argument("nomos::games::all_of: games on the same players");
  return {games[0].players, [games](std::uint32_t S) {
            for (const auto& g : games)
              if (!g.wins(S)) return false;
            return true;
          }};
}

inline SimpleGame any_of(const std::vector<SimpleGame>& games) {
  if (games.empty()) throw std::invalid_argument("nomos::games::any_of: at least one game");
  for (const auto& g : games)
    if (g.players != games[0].players) throw std::invalid_argument("nomos::games::any_of: games on the same players");
  return {games[0].players, [games](std::uint32_t S) {
            for (const auto& g : games)
              if (g.wins(S)) return true;
            return false;
          }};
}

struct Properties {
  bool monotone = true, proper = true, strong = true;
  std::vector<std::size_t> vetoes, dummies;
  std::vector<std::uint32_t> minimal_winning;
  bool decisive() const { return proper && strong; }
};

inline std::vector<char> outcomes(const SimpleGame& g) {
  if (g.players == 0 || g.players > 24) throw std::invalid_argument("nomos::games: between 1 and 24 players for exact enumeration");
  const std::uint32_t sets = std::uint32_t{1} << g.players;
  std::vector<char> w(sets);
  for (std::uint32_t S = 0; S < sets; S++) w[S] = g.wins(S) ? 1 : 0;
  if (w[0]) throw std::invalid_argument("nomos::games: the empty coalition wins");
  if (!w[sets - 1]) throw std::invalid_argument("nomos::games: the grand coalition loses");
  return w;
}

inline Properties properties(const SimpleGame& g) {
  const auto w = outcomes(g);
  const std::size_t n = g.players;
  const std::uint32_t sets = std::uint32_t{1} << n, all = sets - 1;
  Properties p;
  std::vector<char> in_every(n, 1), swing(n, 0);
  for (std::uint32_t S = 0; S < sets; S++) {
    bool minimal = w[S];
    for (std::size_t i = 0; i < n; i++) {
      const std::uint32_t bit = std::uint32_t{1} << i;
      if (S & bit) {
        if (w[S] && w[S & ~bit]) minimal = false;
      } else {
        if (w[S] && !w[S | bit]) p.monotone = false;
        if (!w[S] && w[S | bit]) swing[i] = 1;
      }
      if (w[S] && !(S & bit)) in_every[i] = 0;
    }
    if (minimal) p.minimal_winning.push_back(S);
    const std::uint32_t C = all & ~S;
    if (w[S] && w[C]) p.proper = false;
    if (!w[S] && !w[C]) p.strong = false;
  }
  for (std::size_t i = 0; i < n; i++) {
    if (in_every[i]) p.vetoes.push_back(i);
    if (!swing[i]) p.dummies.push_back(i);
  }
  return p;
}

struct Power {
  std::vector<double> swings;            // η_i (as a real number, exact below 2^53)
  std::vector<double> penrose;           // η_i / 2^(n−1)
  std::vector<double> banzhaf;           // η_i / Σ η
  std::vector<double> shapley_shubik;    // φ_i
  std::vector<double> coleman_prevent;   // η_i / W
  std::vector<double> coleman_initiate;  // η_i / (2^n − W)
  double power_to_act = 0;               // W / 2^n
};

namespace detail {
inline void normalise(Power& p) {
  double total = 0;
  for (double e : p.swings) total += e;
  p.banzhaf.assign(p.swings.size(), 0.0);
  if (total > 0)
    for (std::size_t i = 0; i < p.swings.size(); i++) p.banzhaf[i] = p.swings[i] / total;
}
}  // namespace detail

// The a priori power of any simple game on n ≤ 24 players, exactly over its coalitions.
inline Power power(const SimpleGame& g) {
  const auto w = outcomes(g);
  const std::size_t n = g.players;
  const std::uint32_t sets = std::uint32_t{1} << n;
  Power p;
  p.swings.assign(n, 0.0);
  double winning = 0;
  for (std::uint32_t S = 0; S < sets; S++) {
    winning += w[S];
    for (std::size_t i = 0; i < n; i++) {
      const std::uint32_t bit = std::uint32_t{1} << i;
      if (!(S & bit) && !w[S] && w[S | bit]) p.swings[i] += 1;
    }
  }
  const double half = std::ldexp(1.0, static_cast<int>(n) - 1), all = std::ldexp(1.0, static_cast<int>(n));
  p.penrose.resize(n), p.coleman_prevent.resize(n), p.coleman_initiate.resize(n);
  for (std::size_t i = 0; i < n; i++) {
    p.penrose[i] = p.swings[i] / half;
    p.coleman_prevent[i] = p.swings[i] / winning;
    p.coleman_initiate[i] = p.swings[i] / (all - winning);
  }
  p.power_to_act = winning / all;
  detail::normalise(p);
  p.shapley_shubik = shapley(n, [&](std::uint32_t S) { return w[S] ? 1.0 : 0.0; });
  return p;
}

// The a priori power of a weighted game on any number of players, by generating functions. Counts are kept as probabilities (each
// other player in with chance ½, or each size of coalition weighted by 1/C(n − 1, k)) so that nothing overflows.
inline Power power(const WeightedGame& g) {
  validate(g);
  const std::size_t n = g.weights.size();
  const long q = g.quota;
  Power p;
  p.swings.assign(n, 0.0), p.penrose.assign(n, 0.0), p.shapley_shubik.assign(n, 0.0);
  p.coleman_prevent.assign(n, 0.0), p.coleman_initiate.assign(n, 0.0);
  // the chance that a random coalition (each player in with chance ½) wins: the whole distribution of weight, truncated at q
  {
    std::vector<long double> f(static_cast<std::size_t>(q) + 1, 0.0L);  // f[w] for w < q; f[q] collects every weight ≥ q
    f[0] = 1;
    for (long wi : g.weights) {
      std::vector<long double> h(f.size(), 0.0L);
      for (std::size_t w = 0; w < f.size(); w++) {
        if (f[w] == 0) continue;
        h[w] += f[w] / 2;
        h[std::min<std::size_t>(f.size() - 1, w + static_cast<std::size_t>(wi))] += f[w] / 2;
      }
      f.swap(h);
    }
    p.power_to_act = static_cast<double>(f[static_cast<std::size_t>(q)]);
  }
  for (std::size_t i = 0; i < n; i++) {
    const long wi = g.weights[i];
    if (wi == 0) continue;  // a player of no weight swings nothing
    // b[w]: the chance that the others (each in with chance ½) weigh w, for w < q
    std::vector<long double> b(static_cast<std::size_t>(q), 0.0L);
    // c[k][w]: the share of the C(n − 1, k) coalitions of k others that weigh w, for w < q
    std::vector<std::vector<long double>> c(n, std::vector<long double>(static_cast<std::size_t>(q), 0.0L));
    b[0] = 1, c[0][0] = 1;
    std::size_t seen = 0;
    for (std::size_t j = 0; j < n; j++) {
      if (j == i) continue;
      const long wj = g.weights[j];
      for (long w = q - 1; w >= 0; w--) {
        const long double stay = b[static_cast<std::size_t>(w)] / 2;
        const long double moved = (w - wj >= 0) ? b[static_cast<std::size_t>(w - wj)] / 2 : 0.0L;
        b[static_cast<std::size_t>(w)] = stay + moved;
      }
      // c over m = seen + 1 others: share of size k is (C(seen, k) c_k[w] + C(seen, k − 1) c_{k−1}[w − wj]) / C(seen + 1, k)
      //   = ((seen + 1 − k) c_k[w] + k c_{k−1}[w − wj]) / (seen + 1)
      const long double m = static_cast<long double>(seen + 1);
      for (std::size_t k = std::min(seen + 1, n - 1); k + 1 > 0; k--) {
        for (long w = q - 1; w >= 0; w--) {
          const long double keep = k <= seen ? c[k][static_cast<std::size_t>(w)] * (m - static_cast<long double>(k)) : 0.0L;
          const long double add = (k >= 1 && w - wj >= 0) ? c[k - 1][static_cast<std::size_t>(w - wj)] * static_cast<long double>(k) : 0.0L;
          c[k][static_cast<std::size_t>(w)] = (keep + add) / m;
        }
        if (k == 0) break;
      }
      seen++;
    }
    long double pen = 0, ss = 0;
    for (long w = std::max(0L, q - wi); w < q; w++) {
      pen += b[static_cast<std::size_t>(w)];
      for (std::size_t k = 0; k < n; k++) ss += c[k][static_cast<std::size_t>(w)];
    }
    p.penrose[i] = static_cast<double>(pen);
    p.shapley_shubik[i] = static_cast<double>(ss / static_cast<long double>(n));
    p.swings[i] = static_cast<double>(pen * std::pow(2.0L, static_cast<long double>(n - 1)));
    // Coleman: η_i / W = β′_i 2^(n−1) / (A 2^n) = β′_i / (2A), with A the power to act
    p.coleman_prevent[i] = p.power_to_act > 0 ? p.penrose[i] / (2 * p.power_to_act) : 0;
    p.coleman_initiate[i] = p.power_to_act < 1 ? p.penrose[i] / (2 * (1 - p.power_to_act)) : 0;
  }
  detail::normalise(p);
  return p;
}

// The chance that one voter of m is decisive in a simple majority vote (at least ⌊m/2⌋ + 1 yes votes), the others voting yes or no
// with equal chance and independently: C(m − 1, ⌊m/2⌋) / 2^(m−1), ≈ √(2/(π m)) for large m.
inline double decisive_in_majority(double m) {
  if (!(m >= 1)) throw std::invalid_argument("nomos::games::decisive_in_majority: at least one voter");
  const double others = m - 1, k = std::floor(m / 2);
  if (m == 1) return 1;
  return std::exp(calx::special_functions::log_choose(others, k) - others * std::log(2.0));
}

// A citizen's indirect power in a two-tier system: decisive in their unit's majority and their delegate decisive in the council.
inline std::vector<double> indirect_power(const std::vector<double>& sizes, const Power& council) {
  if (sizes.size() != council.penrose.size()) throw std::invalid_argument("nomos::games::indirect_power: one size per member of the council");
  std::vector<double> out(sizes.size());
  for (std::size_t i = 0; i < sizes.size(); i++) out[i] = decisive_in_majority(sizes[i]) * council.penrose[i];
  return out;
}

}  // namespace nomos::games
