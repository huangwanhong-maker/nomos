// The inequality of a distribution: the Gini coefficient, the Lorenz curve and its generalisation,
// the generalised entropy family (Theil's T and L among it), Atkinson's indices, the coefficient
// of variation, quantile shares and the Palma ratio.
//
// A distribution is a list of amounts x_i held by units with weights w_i ≥ 0 (frequency or
// sampling weights, of positive total W; none given means a weight of one each). Every index here
// is a functional of the weighted empirical distribution F = Σ w_i δ(x_i) / W. Hence each is
// unchanged when the population is replicated, or a unit of weight 2 is split into two of weight
// 1 (Dalton's principle of population; Dalton 1920), and, being relative, when every amount is
// multiplied by the same positive number (scale invariance). Each also satisfies the principle of
// transfers: a transfer from a richer unit to a poorer one that does not reverse their order
// lowers it (Pigou 1912; Dalton 1920). μ = Σ w_i x_i / W is the mean, and r_i = x_i / μ the
// relative amounts.
//
// - The Gini coefficient is half the relative mean difference (Gini 1912; Cowell 2011, ch. 2):
//     G = Σ_i Σ_j w_i w_j |x_i − x_j| / (2 W² μ),
//   computed in O(n log n) as G = 2 cov(x, F̂) / μ with F̂_i = (Σ_{j<i} w_j + w_i / 2) / W the
//   midpoint cumulative weight in increasing order of amount (Lerman and Yitzhaki 1984; with
//   weights, Lerman and Yitzhaki 1989). The two agree exactly for any weights, and
//   for tied amounts in any order: Σ_{i,j} w_i w_j |x_i − x_j| = 2W Σ_i w_i x_i (2F̂_i − 1). This
//   is the population convention: equal amounts give 0, and one unit holding everything among
//   n of equal weight gives (n − 1)/n. (Multiplying by n/(n − 1) gives the estimator of a
//   sample's parent population, as Deaton 1997, ch. 3, writes it; it is not done here.)
// - The Lorenz curve L(p) = (1/μ) ∫_0^p F⁻¹(t) dt (Gastwirth 1971) is, for a discrete
//   distribution, the polygon through (0, 0) and the points (C_k / W, S_k / (W μ)) of cumulative
//   weight C_k and cumulative amount S_k in increasing order; between them it is linear, so
//   interpolating linearly is exact. The generalised Lorenz curve is μ L(p) (Shorrocks 1983).
//   G = 1 − 2 ∫_0^1 L(p) dp, which on the polygon is exactly the trapezoidal sum.
// - Quantile shares are differences of L: the poorest fraction p holds L(p), the richest
//   fraction p holds 1 − L(1 − p). The Palma ratio is the share of the richest 10 per cent over
//   the share of the poorest 40 per cent, (1 − L(0.9)) / L(0.4) (Palma 2011; Cobham and Sumner
//   2013).
// - The generalised entropy family (Shorrocks 1980; Cowell 2011, ch. 3)
//     GE(α) = (E[r^α] − 1) / (α(α − 1))   for α ≠ 0, 1,
//     GE(0) = −E[ln r]   (the mean logarithmic deviation, Theil's L),
//     GE(1) =  E[r ln r] (Theil's T; Theil 1967),
//   with E the weighted mean over the units. GE(2) is half the squared coefficient of variation.
//   α < 0 and α = 0 need every amount positive; α > 0 allows zeros (0 ln 0 = 0).
//   E[r^α] − 1 is evaluated as E[expm1(α ln r)] for α < ½ and as E[r expm1((α − 1) ln r)] for
//   α ≥ ½ (using E[r] = 1), so that neither near α = 0 nor near α = 1 does the difference cancel;
//   the indices are continuous in α and accurate to a few units of rounding times ln(max r).
// - Atkinson's index (Atkinson 1970): one minus the ratio of the equally distributed
//   equivalent amount to the mean,
//     A(ε) = 1 − (E[r^{1−ε}])^{1/(1−ε)}   for ε ≥ 0, ε ≠ 1,      A(1) = 1 − exp(E[ln r]),
//   ε the aversion to inequality. It is an ordinal transform of GE(1 − ε):
//     A(ε) = 1 − (1 + ε(ε − 1) GE(1 − ε))^{1/(1−ε)},   A(1) = 1 − exp(−GE(0))
//   (Cowell 2011, ch. 3), which is how it is computed (with log1p and expm1). ε ≥ 1 needs
//   every amount positive.
// - The coefficient of variation is σ / μ, σ the population standard deviation (weights as
//   frequencies, divisor W).
//
// Accuracy: sums are plain, so each index carries a relative error of order n·u (u the unit
// roundoff), and an index that vanishes (equal amounts) is returned as a few u.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::measurement {

namespace detail {

inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::measurement: ") + what);
}

// The weights, checked: those given, or ones.
inline std::vector<double> weights(std::span<const double> x, std::span<const double> w) {
  if (w.empty()) return std::vector<double>(x.size(), 1.0);
  require(w.size() == x.size(), "as many weights as amounts");
  for (double v : w) require(v >= 0 && std::isfinite(v), "nonnegative, finite weights");
  return {w.begin(), w.end()};
}

// The total weight W and the mean μ.
struct Totals {
  double weight = 0, mean = 0;
};

inline Totals totals(std::span<const double> x, const std::vector<double>& w) {
  require(!x.empty(), "a nonempty distribution");
  Totals t;
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    require(std::isfinite(x[i]), "finite amounts");
    t.weight += w[i];
    s += w[i] * x[i];
  }
  require(t.weight > 0, "a positive total weight");
  t.mean = s / t.weight;
  return t;
}

// The positions of x in increasing order of amount (stable: ties keep their order).
inline std::vector<std::size_t> ascending(std::span<const double> x) {
  std::vector<std::size_t> o(x.size());
  std::iota(o.begin(), o.end(), std::size_t{0});
  std::stable_sort(o.begin(), o.end(), [&](std::size_t a, std::size_t b) { return x[a] < x[b]; });
  return o;
}

}  // namespace detail

// μ = Σ w_i x_i / W.
inline double mean(std::span<const double> x, std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  return detail::totals(x, ww).mean;
}

// The Gini coefficient G = Σ_ij w_i w_j |x_i − x_j| / (2 W² μ), by the covariance with the
// midpoint cumulative weights (Lerman and Yitzhaki 1989). Needs μ > 0. O(n log n).
inline double gini(std::span<const double> x, std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  detail::require(t.mean > 0, "a positive mean");
  const auto o = detail::ascending(x);
  double below = 0, s = 0;  // the weight below the current unit; Σ w_i x_i (2F̂_i − 1)
  for (std::size_t i : o) {
    const double F = (below + 0.5 * ww[i]) / t.weight;
    s += ww[i] * x[i] * (2 * F - 1);
    below += ww[i];
  }
  return s / (t.weight * t.mean);
}

// The Lorenz curve of a distribution: the polygon through (p[k], L[k]), from (0, 0) to (1, 1).
struct LorenzCurve {
  std::vector<double> p;  // cumulative shares of the weight, increasing from 0 to 1
  std::vector<double> L;  // cumulative shares of the total amount, from 0 to 1
  double mean = 0;        // μ; the generalised Lorenz curve is μ L(p)

  // L(q) for q in [0, 1], linear between the vertices (exact for a discrete distribution).
  double at(double q) const {
    detail::require(q >= 0 && q <= 1, "a share of the population in [0, 1]");
    const auto k = static_cast<std::size_t>(std::upper_bound(p.begin(), p.end(), q) - p.begin());
    if (k == 0) return L.front();
    if (k >= p.size()) return L.back();
    const double d = p[k] - p[k - 1];
    return d > 0 ? L[k - 1] + (L[k] - L[k - 1]) * (q - p[k - 1]) / d : L[k];
  }
  // The generalised Lorenz curve μ L(q) (Shorrocks 1983).
  double generalised_at(double q) const { return mean * at(q); }
};

// The Lorenz curve. Needs μ > 0.
inline LorenzCurve lorenz(std::span<const double> x, std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  detail::require(t.mean > 0, "a positive mean");
  const double total = t.weight * t.mean;
  LorenzCurve c;
  c.mean = t.mean;
  c.p.reserve(x.size() + 1);
  c.L.reserve(x.size() + 1);
  c.p.push_back(0);
  c.L.push_back(0);
  double cw = 0, cx = 0;
  for (std::size_t i : detail::ascending(x)) {
    cw += ww[i];
    cx += ww[i] * x[i];
    c.p.push_back(cw / t.weight);
    c.L.push_back(cx / total);
  }
  c.p.back() = 1;  // exact ends, whatever the rounding of the cumulative sums
  c.L.back() = 1;
  return c;
}

// G = 1 − 2 ∫ L, by trapezoids on the polygon (exact): a second computation of the Gini coefficient.
inline double gini(const LorenzCurve& c) {
  double area = 0;
  for (std::size_t k = 1; k < c.p.size(); k++) area += (c.p[k] - c.p[k - 1]) * (c.L[k] + c.L[k - 1]);
  return 1 - area;
}

// The share of the total held by the units between population shares p0 < p1 (from the poorest).
inline double quantile_share(const LorenzCurve& c, double p0, double p1) {
  detail::require(p0 <= p1, "p0 ≤ p1");
  return c.at(p1) - c.at(p0);
}

// The share held by the richest fraction q of the population: 1 − L(1 − q).
inline double top_share(const LorenzCurve& c, double q) { return 1 - c.at(1 - q); }

// The Palma ratio: the share of the richest 10 per cent over that of the poorest 40 per cent.
inline double palma_ratio(const LorenzCurve& c) { return (1 - c.at(0.9)) / c.at(0.4); }
inline double palma_ratio(std::span<const double> x, std::span<const double> w = {}) { return palma_ratio(lorenz(x, w)); }

// The generalised entropy index GE(α) (Shorrocks 1980). α ≤ 0 needs positive amounts; α > 0
// allows zeros; every α needs μ > 0.
inline double generalised_entropy(std::span<const double> x, double alpha, std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  detail::require(t.mean > 0, "a positive mean");
  detail::require(std::isfinite(alpha), "a finite α");
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++) {
    if (ww[i] == 0) continue;
    const double r = x[i] / t.mean;
    detail::require(r > 0 || (r == 0 && alpha > 0), "positive amounts (zeros only for α > 0)");
    if (alpha == 0) {
      s -= ww[i] * std::log(r);
    } else if (alpha == 1) {
      if (r > 0) s += ww[i] * r * std::log(r);
    } else if (alpha < 0.5) {
      s += ww[i] * (r > 0 ? std::expm1(alpha * std::log(r)) : -1.0);  // r^α − 1
    } else {
      s += ww[i] * (r > 0 ? r * std::expm1((alpha - 1) * std::log(r)) : 0.0);  // r^α − r
    }
  }
  s /= t.weight;
  return (alpha == 0 || alpha == 1) ? s : s / (alpha * (alpha - 1));
}

// Theil's T, GE(1).
inline double theil_t(std::span<const double> x, std::span<const double> w = {}) { return generalised_entropy(x, 1, w); }
// Theil's L, the mean logarithmic deviation, GE(0).
inline double theil_l(std::span<const double> x, std::span<const double> w = {}) { return generalised_entropy(x, 0, w); }
inline double mean_log_deviation(std::span<const double> x, std::span<const double> w = {}) { return generalised_entropy(x, 0, w); }

// Atkinson's index A(ε), ε ≥ 0 (Atkinson 1970), through GE(1 − ε).
inline double atkinson(std::span<const double> x, double epsilon, std::span<const double> w = {}) {
  detail::require(epsilon >= 0 && std::isfinite(epsilon), "an aversion ε ≥ 0");
  if (epsilon == 0) return 0;
  if (epsilon == 1) return -std::expm1(-generalised_entropy(x, 0, w));
  const double ge = generalised_entropy(x, 1 - epsilon, w);
  return -std::expm1(std::log1p(epsilon * (epsilon - 1) * ge) / (1 - epsilon));
}

// The coefficient of variation σ / μ (σ with divisor W).
inline double coefficient_of_variation(std::span<const double> x, std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  detail::require(t.mean > 0, "a positive mean");
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++) s += ww[i] * (x[i] - t.mean) * (x[i] - t.mean);
  return std::sqrt(s / t.weight) / t.mean;
}

}  // namespace nomos::measurement
