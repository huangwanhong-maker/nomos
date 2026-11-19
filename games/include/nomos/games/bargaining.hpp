// Bargaining between two parties, and the fair division of issues between them.
//
// A bargaining problem (S, d) is a compact convex set S of pairs of utilities the parties can reach together, and a disagreement
// point d ∈ S, what each has if they do not agree (Nash 1950). Here S is the convex hull of finitely many pairs (the outcomes the
// parties can agree to, and their lotteries), and a solution is a point of its Pareto frontier at or above d:
//
//   Nash                  the point maximising (u₁ − d₁)(u₂ − d₂): the one solution that is Pareto optimal, symmetric, invariant
//                         under positive affine changes of either party's utility, and independent of irrelevant alternatives
//                         (Nash 1950)
//   Kalai–Smorodinsky     the Pareto point on the segment from d to the ideal point I (each party's best within S at or above d):
//                         it trades independence of irrelevant alternatives for individual monotonicity, that a party whose
//                         opportunities grow at every level of the other's does not lose (Kalai and Smorodinsky 1975)
//   egalitarian           the Pareto point with equal gains u₁ − d₁ = u₂ − d₂: monotone in S, but not invariant under changes of
//                         scale, so it needs utilities that are comparable between the parties (Kalai 1977)
//
// Each presupposes cardinal utilities the parties declare; the survey behind this module (the infrastructure's diplomacy domain)
// notes that they seldom do, which is why these are tools for a declared model, never a measure of what was gained.
//
// Adjusted Winner (Brams and Taylor 1996) divides issues between two parties who each spread 100 points over them. Each issue
// goes first to the party that values it more; then issues pass from the party ahead, in the increasing order of the ratio of its
// points to the other's, until the totals are equal, the last issue split. The division is equitable (equal points), envy-free
// (each values its share at least as much as the other's) and efficient (no division of the issues, split or not, is better for
// one and no worse for the other); the checks hold the last to a linear programme.
//
// References. Nash, J. F. (1950), The bargaining problem, Econometrica 18, 155–162. Kalai, E. and Smorodinsky, M. (1975), Other
// solutions to Nash's bargaining problem, Econometrica 43, 513–518. Kalai, E. (1977), Proportional solutions to bargaining
// situations: interpersonal utility comparisons, Econometrica 45, 1623–1630. Brams, S. J. and Taylor, A. D. (1996), Fair Division:
// From Cake-Cutting to Dispute Resolution, Cambridge University Press, ch. 4. Osborne, M. J. and Rubinstein, A. (1990), Bargaining
// and Markets, Academic Press, ch. 2.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nomos::games {

using Pair = std::array<double, 2>;

struct Bargaining {
  std::vector<Pair> outcomes;  // S is their convex hull
  Pair disagreement{0, 0};
};

struct BargainingSolution {
  Pair point{0, 0};
  Pair ideal{0, 0};
  std::vector<Pair> frontier;  // the Pareto frontier at or above d, from the most for party 1's partner to the most for party 1
};

namespace detail {

inline double cross(const Pair& o, const Pair& a, const Pair& b) { return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]); }

// The Pareto frontier of the convex hull of the points, from its point of greatest u₂ to its point of greatest u₁ (Andrew's
// monotone chain for the upper hull).
inline std::vector<Pair> pareto_frontier(std::vector<Pair> p) {
  std::sort(p.begin(), p.end());
  p.erase(std::unique(p.begin(), p.end()), p.end());
  std::vector<Pair> upper;
  for (auto it = p.rbegin(); it != p.rend(); ++it) {  // right to left, keeping right turns: the upper hull from the right
    while (upper.size() >= 2 && cross(upper[upper.size() - 2], upper.back(), *it) <= 0) upper.pop_back();
    upper.push_back(*it);
  }
  std::reverse(upper.begin(), upper.end());  // left to right
  // keep from the highest point (the rightmost among the highest) to the rightmost point (the highest among the rightmost)
  std::size_t top = 0;
  for (std::size_t i = 0; i < upper.size(); i++)
    if (upper[i][1] >= upper[top][1]) top = i;
  std::vector<Pair> f(upper.begin() + static_cast<std::ptrdiff_t>(top), upper.end());
  return f;
}

// The frontier clipped to u ≥ d.
inline std::vector<Pair> clip(const std::vector<Pair>& f, const Pair& d) {
  std::vector<Pair> out;
  for (std::size_t i = 0; i < f.size(); i++) {
    const Pair& a = f[i];
    const bool in = a[0] >= d[0] && a[1] >= d[1];
    if (in) out.push_back(a);
    if (i + 1 == f.size()) break;
    const Pair& b = f[i + 1];
    // where the segment a–b crosses u₁ = d₁ (u₁ rises along the frontier) or u₂ = d₂ (u₂ falls)
    std::vector<std::pair<double, Pair>> cuts;
    if ((a[0] - d[0]) * (b[0] - d[0]) < 0) {
      const double t = (d[0] - a[0]) / (b[0] - a[0]);
      cuts.push_back({t, {d[0], a[1] + t * (b[1] - a[1])}});
    }
    if ((a[1] - d[1]) * (b[1] - d[1]) < 0) {
      const double t = (d[1] - a[1]) / (b[1] - a[1]);
      cuts.push_back({t, {a[0] + t * (b[0] - a[0]), d[1]}});
    }
    std::sort(cuts.begin(), cuts.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    for (const auto& [t, c] : cuts)
      if (c[0] >= d[0] - 1e-12 && c[1] >= d[1] - 1e-12) out.push_back(c);
  }
  return out;
}

inline BargainingSolution prepare(const Bargaining& B) {
  if (B.outcomes.empty()) throw std::invalid_argument("nomos::games::Bargaining: at least one outcome");
  BargainingSolution s;
  s.frontier = clip(pareto_frontier(B.outcomes), B.disagreement);
  if (s.frontier.empty()) throw std::invalid_argument("nomos::games::Bargaining: no outcome is at or above the disagreement point");
  s.ideal = {s.frontier.back()[0], s.frontier.front()[1]};
  return s;
}

// The frontier's point on the ray from d in direction r (r ≥ 0, not zero), or the frontier's end it reaches.
inline Pair on_ray(const std::vector<Pair>& f, const Pair& d, const Pair& r) {
  if (f.size() == 1) return f[0];
  for (std::size_t i = 0; i + 1 < f.size(); i++) {
    const Pair &a = f[i], &b = f[i + 1];
    // solve d + λ r = a + t (b − a)
    const double det = r[0] * (a[1] - b[1]) - r[1] * (a[0] - b[0]);
    if (std::abs(det) < 1e-300) continue;
    const double t = (r[1] * (d[0] - a[0]) - r[0] * (d[1] - a[1])) / det;
    if (t >= -1e-12 && t <= 1 + 1e-12) return {a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])};
  }
  // the ray passes the frontier's end: the nearest end along it
  return r[1] * (f.front()[0] - d[0]) >= r[0] * (f.front()[1] - d[1]) ? f.front() : f.back();
}

}  // namespace detail

inline BargainingSolution nash(const Bargaining& B) {
  auto s = detail::prepare(B);
  const Pair& d = B.disagreement;
  auto product = [&](const Pair& u) { return (u[0] - d[0]) * (u[1] - d[1]); };
  Pair best = s.frontier[0];
  double most = product(best);
  for (std::size_t i = 0; i < s.frontier.size(); i++) {
    if (product(s.frontier[i]) > most) most = product(s.frontier[i]), best = s.frontier[i];
    if (i + 1 == s.frontier.size()) break;
    const Pair &a = s.frontier[i], &b = s.frontier[i + 1];
    const double x = a[0] - d[0], y = a[1] - d[1], dx = b[0] - a[0], dy = b[1] - a[1];
    // f(t) = (x + t dx)(y + t dy), stationary at t = −(x dy + y dx)/(2 dx dy)
    if (dx * dy != 0) {
      const double t = -(x * dy + y * dx) / (2 * dx * dy);
      if (t > 0 && t < 1) {
        const Pair u{a[0] + t * dx, a[1] + t * dy};
        if (product(u) > most) most = product(u), best = u;
      }
    }
  }
  s.point = best;
  return s;
}

inline BargainingSolution kalai_smorodinsky(const Bargaining& B) {
  auto s = detail::prepare(B);
  const Pair& d = B.disagreement;
  const Pair r{s.ideal[0] - d[0], s.ideal[1] - d[1]};
  s.point = (r[0] <= 0 && r[1] <= 0) ? d : detail::on_ray(s.frontier, d, r);
  return s;
}

inline BargainingSolution egalitarian(const Bargaining& B) {
  auto s = detail::prepare(B);
  s.point = detail::on_ray(s.frontier, B.disagreement, {1, 1});
  return s;
}

// ---- Adjusted Winner ----------------------------------------------------------------------------------------------------------------

struct Division {
  std::vector<double> share;   // per issue: the share of it party A holds, in [0, 1] (B holds the rest)
  double points_a = 0, points_b = 0;
  std::size_t split = static_cast<std::size_t>(-1);  // the issue split, if any
};

inline Division adjusted_winner(const std::vector<double>& a, const std::vector<double>& b) {
  const std::size_t m = a.size();
  if (m == 0 || b.size() != m) throw std::invalid_argument("nomos::games::adjusted_winner: the same issues for both parties");
  double ta = 0, tb = 0;
  for (std::size_t g = 0; g < m; g++) {
    if (!(a[g] >= 0) || !(b[g] >= 0)) throw std::invalid_argument("nomos::games::adjusted_winner: points are not negative");
    ta += a[g], tb += b[g];
  }
  if (std::abs(ta - 100) > 1e-9 || std::abs(tb - 100) > 1e-9) throw std::invalid_argument("nomos::games::adjusted_winner: each party spreads 100 points");
  Division D;
  D.share.assign(m, 0.0);
  for (std::size_t g = 0; g < m; g++) D.share[g] = a[g] >= b[g] ? 1.0 : 0.0;  // a tie goes first to A, and is the first to pass
  auto totals = [&] {
    D.points_a = D.points_b = 0;
    for (std::size_t g = 0; g < m; g++) D.points_a += D.share[g] * a[g], D.points_b += (1 - D.share[g]) * b[g];
  };
  totals();
  // the party ahead gives up issues in the increasing order of its points over the other's
  const bool a_ahead = D.points_a >= D.points_b;
  const std::vector<double>& mine = a_ahead ? a : b;
  const std::vector<double>& theirs = a_ahead ? b : a;
  std::vector<std::size_t> held;
  for (std::size_t g = 0; g < m; g++)
    if ((a_ahead ? D.share[g] : 1 - D.share[g]) == 1.0) held.push_back(g);
  std::sort(held.begin(), held.end(), [&](std::size_t x, std::size_t y) {
    // mine[x]/theirs[x] < mine[y]/theirs[y], without dividing by zero
    return mine[x] * theirs[y] < mine[y] * theirs[x] || (mine[x] * theirs[y] == mine[y] * theirs[x] && x < y);
  });
  for (std::size_t g : held) {
    const double lead = std::abs(D.points_a - D.points_b);
    if (lead <= 1e-12) break;
    // passing the whole issue changes the lead by mine[g] + theirs[g]
    const double swing = mine[g] + theirs[g];
    if (swing <= 0) continue;
    const double part = std::min(1.0, lead / swing);  // the part of the issue that passes
    if (a_ahead) D.share[g] = 1 - part;
    else D.share[g] = part;
    totals();
    if (part < 1) {
      D.split = g;
      break;
    }
  }
  totals();
  return D;
}

}  // namespace nomos::games
