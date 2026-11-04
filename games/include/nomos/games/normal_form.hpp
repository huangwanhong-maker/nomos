// Finite games of n players in normal form, and their pure Nash equilibria.
//
// Player p chooses one of s_p strategies; a profile is one strategy for each player, and every
// profile gives every player a payoff. A profile is a pure Nash equilibrium (Nash 1951) when no
// player gains by changing only their own strategy. Profiles are numbered in mixed radix with
// player 0 the most significant digit, so that a two-player game's profile (i, j) is i·s_1 + j,
// the row-major cell of its bimatrix.
//
// Enumeration checks every profile against every unilateral deviation: O(Π s_p · Σ s_p), for
// small games. A pure equilibrium need not exist (matching pennies); in a finite potential game
// one always does (Monderer and Shapley 1996).

#pragma once

#include <numeric>
#include <span>
#include <vector>

#include "bimatrix.hpp"

namespace nomos::games {

struct NormalForm {
  std::vector<std::size_t> strategies;  // s_p for each player
  std::vector<double> payoffs;          // payoffs[index(profile) · players + p]

  NormalForm() = default;
  explicit NormalForm(std::vector<std::size_t> counts) : strategies(std::move(counts)) {
    detail::require(!strategies.empty(), "at least one player");
    std::size_t n = 1;
    for (std::size_t s : strategies) {
      detail::require(s > 0, "at least one strategy for each player");
      n *= s;
    }
    payoffs.assign(n * strategies.size(), 0.0);
  }
  std::size_t players() const { return strategies.size(); }
  std::size_t profiles() const { return payoffs.size() / strategies.size(); }
  std::size_t index(std::span<const std::size_t> profile) const {
    detail::require(profile.size() == players(), "a strategy for each player");
    std::size_t k = 0;
    for (std::size_t p = 0; p < players(); p++) {
      detail::require(profile[p] < strategies[p], "a strategy within range");
      k = k * strategies[p] + profile[p];
    }
    return k;
  }
  std::vector<std::size_t> profile(std::size_t k) const {
    std::vector<std::size_t> s(players());
    for (std::size_t p = players(); p-- > 0;) {
      s[p] = k % strategies[p];
      k /= strategies[p];
    }
    return s;
  }
  double& payoff(std::span<const std::size_t> profile, std::size_t player) { return payoffs[index(profile) * players() + player]; }
  double payoff(std::span<const std::size_t> profile, std::size_t player) const { return payoffs[index(profile) * players() + player]; }
};

// The two-player normal form of a bimatrix game.
inline NormalForm normal_form(const Bimatrix& g) {
  NormalForm f({g.rows(), g.cols()});
  for (std::size_t i = 0; i < g.rows(); i++)
    for (std::size_t j = 0; j < g.cols(); j++) {
      const std::size_t s[2] = {i, j};
      f.payoff(s, 0) = g.a(i, j);
      f.payoff(s, 1) = g.b(i, j);
    }
  return f;
}

// The pure Nash equilibria, in increasing order of profile number.
inline std::vector<std::vector<std::size_t>> pure_equilibria(const NormalForm& g, double tolerance = 1e-9) {
  std::vector<std::vector<std::size_t>> out;
  for (std::size_t k = 0; k < g.profiles(); k++) {
    auto s = g.profile(k);
    bool stable = true;
    for (std::size_t p = 0; p < g.players() && stable; p++) {
      const double here = g.payoffs[k * g.players() + p];
      const std::size_t own = s[p];
      for (std::size_t t = 0; t < g.strategies[p] && stable; t++) {
        if (t == own) continue;
        s[p] = t;
        if (g.payoff(s, p) > here + tolerance) stable = false;
      }
      s[p] = own;
    }
    if (stable) out.push_back(s);
  }
  return out;
}

}  // namespace nomos::games
