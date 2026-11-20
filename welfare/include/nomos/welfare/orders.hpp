// Principles as orders over whole allocations, and what can be said when they disagree.
//
// An allocation's outcome is a vector u of what each unit (a person, a household, a place) ends with, with population weights
// w. A principle orders such outcomes; principles conflict, and no complete aggregation of them escapes one of Fleurbaey and
// Tungodden's two "tyrannies" (2010, Soc. Choice Welf. 35, 399–407, as the survey records them: tiny gains to the worst off
// outweighing everything, or small gains to many outweighing a large loss to the worst off). So nothing here adds principles
// into one score unless a person declares the weights, and even then the undominated set is reported beside it (Sen 1997,
// Econometrica 65, 745–779, on maximal against optimal choice). The infrastructure's survey *Allocative justice under
// scarcity* (Part I §2, §7; Part II §7, §9) collects what this rests on.
//
// The principles, each a function of the outcome (higher better unless said):
//
//   total            Σ w u: the utilitarian sum
//   prioritarian(ε)  Atkinson's equally distributed equivalent, the level that, given to all, is as good as u under the
//                    concave transform with inequality aversion ε (Atkinson 1970, J. Econ. Theory 2, 244–263); ε = 0 is the
//                    mean, ε → ∞ the minimum; it is Adler's prioritarianism with the same γ (Adler 2012)
//   leximin          the worst off first, then the next (a lexicographic order of sorted outcomes; Rawls's maximin refined)
//   shortfall(τ, α)  the Foster–Greer–Thorbecke gap below a threshold τ, Σ w max(0, (τ − u)/τ)^α / Σ w (Foster, Greer and
//                    Thorbecke 1984, Econometrica 52, 761–766): α = 0 the share below τ, 1 the mean gap, 2 its severity;
//                    lower is better; sufficientarianism's measure
//   opportunity      the least mean outcome over declared types of circumstance: the worst-off type first, a simplification of
//                    Roemer's equality of opportunity (Roemer 1998), which compares tranches of effort within types; said so
//
// And the dominance relations that give partial orders robust to a principle's parameter: Lorenz dominance (equal means,
// every share of the poorest at least as large), generalised Lorenz dominance (every cumulative total of the poorest at least
// as large; Shorrocks 1983, Economica 50, 3–17), which implies a higher equally distributed equivalent for every ε ≥ 0, and
// sequential dominance for groups ordered by need (Atkinson and Bourguignon 1987: generalised Lorenz dominance of the neediest
// group, of the two neediest pooled, and so on).
#pragma once

#include <nomos/measurement/inequality.hpp>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::welfare {

struct Outcome {
  std::string name;
  std::vector<double> value;     // per unit
  std::vector<double> weight;    // per unit (population); empty: all one
};

enum class Kind { total, prioritarian, leximin, shortfall, opportunity };

struct Principle {
  std::string name;
  Kind kind = Kind::total;
  double parameter = 0;          // prioritarian: ε; shortfall: τ
  double alpha = 1;              // shortfall: α
  std::vector<std::size_t> types;  // opportunity: the type of each unit
  bool lower_better() const { return kind == Kind::shortfall; }
};

namespace detail {

inline std::vector<double> weights(const Outcome& o) {
  if (o.weight.empty()) return std::vector<double>(o.value.size(), 1.0);
  if (o.weight.size() != o.value.size()) throw std::invalid_argument("nomos::welfare: a weight for every unit or none");
  return o.weight;
}

}  // namespace detail

// The principle's value of an outcome (for leximin, the minimum; compare leximin with `compare`).
inline double evaluate(const Principle& p, const Outcome& o) {
  const auto w = detail::weights(o);
  const double W = std::accumulate(w.begin(), w.end(), 0.0);
  if (!(W > 0)) throw std::invalid_argument("nomos::welfare: weights that are not all zero");
  switch (p.kind) {
    case Kind::total: {
      double s = 0;
      for (std::size_t i = 0; i < o.value.size(); i++) s += w[i] * o.value[i];
      return s;
    }
    case Kind::prioritarian: {
      const double mean = measurement::mean(o.value, w);
      if (p.parameter == 0) return mean;
      for (double v : o.value)
        if (!(v > 0)) throw std::invalid_argument("nomos::welfare: a prioritarian value with aversion above zero needs outcomes above zero");
      return mean * (1 - measurement::atkinson(o.value, p.parameter, w));
    }
    case Kind::leximin: return *std::min_element(o.value.begin(), o.value.end());
    case Kind::shortfall: {
      if (!(p.parameter > 0)) throw std::invalid_argument("nomos::welfare: a shortfall threshold above zero");
      double s = 0;
      for (std::size_t i = 0; i < o.value.size(); i++) {
        const double g = std::max(0.0, (p.parameter - o.value[i]) / p.parameter);
        s += w[i] * (p.alpha == 0 ? (g > 0 ? 1.0 : 0.0) : std::pow(g, p.alpha));
      }
      return s / W;
    }
    case Kind::opportunity: {
      if (p.types.size() != o.value.size()) throw std::invalid_argument("nomos::welfare: a type for every unit");
      const std::size_t T = *std::max_element(p.types.begin(), p.types.end()) + 1;
      std::vector<double> sum(T, 0), mass(T, 0);
      for (std::size_t i = 0; i < o.value.size(); i++) sum[p.types[i]] += w[i] * o.value[i], mass[p.types[i]] += w[i];
      double least = 1e300;
      for (std::size_t t = 0; t < T; t++)
        if (mass[t] > 0) least = std::min(least, sum[t] / mass[t]);
      return least;
    }
  }
  return 0;
}

// +1 when a is better than b under p, −1 when worse, 0 when equal within a relative tolerance.
inline int compare(const Principle& p, const Outcome& a, const Outcome& b, double tol = 1e-9) {
  if (p.kind == Kind::leximin) {
    auto x = a.value, y = b.value;
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    for (std::size_t i = 0; i < std::min(x.size(), y.size()); i++) {
      const double s = std::max({1.0, std::abs(x[i]), std::abs(y[i])});
      if (x[i] > y[i] + tol * s) return 1;
      if (x[i] < y[i] - tol * s) return -1;
    }
    return 0;
  }
  const double va = evaluate(p, a), vb = evaluate(p, b), s = std::max({1.0, std::abs(va), std::abs(vb)});
  if (std::abs(va - vb) <= tol * s) return 0;
  const bool a_better = p.lower_better() ? va < vb : va > vb;
  return a_better ? 1 : -1;
}

// ---- dominance --------------------------------------------------------------------------------------------------

// The generalised Lorenz curve at the population shares of both outcomes' breakpoints: the cumulative total of the poorest
// share q, per head of the whole population.
inline double generalised_lorenz_at(const Outcome& o, double q) {
  const auto w = detail::weights(o);
  std::vector<std::size_t> idx(o.value.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return o.value[a] < o.value[b]; });
  const double W = std::accumulate(w.begin(), w.end(), 0.0);
  double mass = 0, total = 0;
  for (std::size_t i : idx) {
    const double share = w[i] / W;
    if (mass + share >= q) return total + (q - mass) * o.value[i];
    mass += share, total += share * o.value[i];
  }
  return total;
}

enum class Relation { dominates, dominated, equal, unordered };

inline const char* name(Relation r) {
  switch (r) {
    case Relation::dominates: return "dominates";
    case Relation::dominated: return "is dominated by";
    case Relation::equal: return "equals";
    default: return "is not comparable with";
  }
}

namespace detail {

inline std::vector<double> breakpoints(const Outcome& a, const Outcome& b) {
  std::vector<double> qs{0, 1};
  for (const Outcome* o : {&a, &b}) {
    const auto w = weights(*o);
    const double W = std::accumulate(w.begin(), w.end(), 0.0);
    std::vector<std::size_t> idx(o->value.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t x, std::size_t y) { return o->value[x] < o->value[y]; });
    double mass = 0;
    for (std::size_t i : idx) qs.push_back(mass += w[i] / W);
  }
  std::sort(qs.begin(), qs.end());
  return qs;
}

inline Relation from_curves(const std::vector<double>& fa, const std::vector<double>& fb, double tol) {
  bool above = false, below = false;
  for (std::size_t k = 0; k < fa.size(); k++) {
    const double s = std::max({1.0, std::abs(fa[k]), std::abs(fb[k])});
    if (fa[k] > fb[k] + tol * s) above = true;
    if (fa[k] < fb[k] - tol * s) below = true;
  }
  if (above && below) return Relation::unordered;
  if (above) return Relation::dominates;
  if (below) return Relation::dominated;
  return Relation::equal;
}

}  // namespace detail

inline Relation generalised_lorenz(const Outcome& a, const Outcome& b, double tol = 1e-9) {
  const auto qs = detail::breakpoints(a, b);
  std::vector<double> fa, fb;
  for (double q : qs) fa.push_back(generalised_lorenz_at(a, q)), fb.push_back(generalised_lorenz_at(b, q));
  return detail::from_curves(fa, fb, tol);
}

// Lorenz dominance compares shares of the total: each generalised Lorenz curve divided by its mean.
inline Relation lorenz(const Outcome& a, const Outcome& b, double tol = 1e-9) {
  const auto qs = detail::breakpoints(a, b);
  const double ma = generalised_lorenz_at(a, 1), mb = generalised_lorenz_at(b, 1);
  if (!(ma > 0) || !(mb > 0)) throw std::invalid_argument("nomos::welfare::lorenz: outcomes with positive totals");
  std::vector<double> fa, fb;
  for (double q : qs) fa.push_back(generalised_lorenz_at(a, q) / ma), fb.push_back(generalised_lorenz_at(b, q) / mb);
  return detail::from_curves(fa, fb, tol);
}

// Sequential dominance for groups ordered by need (group 0 the neediest): a dominates b when, for every k, the units of groups
// 0 … k pooled are generalised-Lorenz dominated by none of b's; the relation holds at every stage or not at all.
inline Relation sequential(const Outcome& a, const Outcome& b, const std::vector<std::size_t>& group, double tol = 1e-9) {
  if (group.size() != a.value.size() || group.size() != b.value.size()) throw std::invalid_argument("nomos::welfare::sequential: a group for every unit, the same units in both");
  const std::size_t G = *std::max_element(group.begin(), group.end()) + 1;
  const auto wa = detail::weights(a), wb = detail::weights(b);
  double total_a = 0, total_b = 0;
  for (double x : wa) total_a += x;
  for (double x : wb) total_b += x;
  bool above = false, below = false;
  for (std::size_t k = 0; k < G; k++) {
    Outcome pa{a.name, {}, {}}, pb{b.name, {}, {}};
    for (std::size_t i = 0; i < group.size(); i++)
      if (group[i] <= k) pa.value.push_back(a.value[i]), pa.weight.push_back(wa[i]), pb.value.push_back(b.value[i]), pb.weight.push_back(wb[i]);
    if (pa.value.empty()) continue;
    // Pooled curves on the whole population's scale, so that stages are comparable.
    const auto qs = detail::breakpoints(pa, pb);
    double ma = 0, mb = 0;
    for (double x : pa.weight) ma += x;
    for (double x : pb.weight) mb += x;
    for (double q : qs) {
      const double fa = generalised_lorenz_at(pa, q) * ma / total_a, fb = generalised_lorenz_at(pb, q) * mb / total_b;
      const double s = std::max({1.0, std::abs(fa), std::abs(fb)});
      if (fa > fb + tol * s) above = true;
      if (fa < fb - tol * s) below = true;
    }
  }
  if (above && below) return Relation::unordered;
  if (above) return Relation::dominates;
  if (below) return Relation::dominated;
  return Relation::equal;
}

}  // namespace nomos::welfare
