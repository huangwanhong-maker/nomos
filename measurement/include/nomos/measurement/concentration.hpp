// Concentration: how a total is shared among few or many holders.
//
// Given nonnegative amounts x_1, …, x_n of positive total, s_i = x_i / Σ x are the shares. The
// indices here depend on the shares alone, not on the size of the total, and are symmetric in
// the holders (Hall and Tideman 1967 set out these and further criteria).
//
// - The Herfindahl–Hirschman index H = Σ s_i² (Hirschman 1945, 1964; Herfindahl 1950) lies in
//   [1/n, 1]: 1/n for equal shares, 1 for a single holder. It is often written on shares in per
//   cent, 10⁴ H; it is returned here on fractions. Its normalised form (H − 1/n)/(1 − 1/n) lies
//   in [0, 1]. 1/H is the numbers equivalent: the number of equal holders with the same H
//   (Adelman 1969).
// - The concentration ratio CR_k is the sum of the k largest shares.
// - The entropy index E = −Σ s_i ln s_i (Theil 1967) lies in [0, ln n]; it rises as the
//   total is spread more evenly (it measures diversity, the opposite of concentration).
// - Hannah and Kay's numbers equivalents HK(α) = (Σ s_i^α)^{1/(1−α)}, α > 0, α ≠ 1 (Hannah and
//   Kay 1977), unite them: HK(2) = 1/H, and HK(α) → exp(E) as α → 1, the value returned
//   at α = 1. Each lies in [1, n] and decreases in α.
//
// Accuracy: relative errors of order n u.

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <span>
#include <vector>

#include "inequality.hpp"

namespace nomos::measurement {

// s_i = x_i / Σ x, for x ≥ 0 of positive total.
inline std::vector<double> shares(std::span<const double> x) {
  detail::require(!x.empty(), "at least one holder");
  double total = 0;
  for (double v : x) {
    detail::require(v >= 0 && std::isfinite(v), "nonnegative, finite amounts");
    total += v;
  }
  detail::require(total > 0, "a positive total");
  std::vector<double> s(x.size());
  for (std::size_t i = 0; i < x.size(); i++) s[i] = x[i] / total;
  return s;
}

// H = Σ s², in [1/n, 1].
inline double herfindahl_hirschman(std::span<const double> x) {
  double h = 0;
  for (double s : shares(x)) h += s * s;
  return h;
}

// (H − 1/n)/(1 − 1/n), in [0, 1]; 0 for a single holder's list of one.
inline double herfindahl_hirschman_normalised(std::span<const double> x) {
  const double n = static_cast<double>(x.size());
  if (x.size() == 1) return 0;
  return (herfindahl_hirschman(x) - 1 / n) / (1 - 1 / n);
}

// CR_k: the sum of the k largest shares (all of them when k ≥ n).
inline double concentration_ratio(std::span<const double> x, std::size_t k) {
  auto s = shares(x);
  k = std::min(k, s.size());
  std::partial_sort(s.begin(), s.begin() + static_cast<std::ptrdiff_t>(k), s.end(), std::greater<>());
  double c = 0;
  for (std::size_t i = 0; i < k; i++) c += s[i];
  return c;
}

// E = −Σ s ln s, in [0, ln n].
inline double entropy_index(std::span<const double> x) {
  double e = 0;
  for (double s : shares(x))
    if (s > 0) e -= s * std::log(s);
  return e;
}

// Hannah and Kay's numbers equivalent (Σ s^α)^{1/(1−α)}, α > 0; exp(E) at α = 1.
inline double hannah_kay(std::span<const double> x, double alpha) {
  detail::require(alpha > 0 && std::isfinite(alpha), "an exponent α > 0");
  if (alpha == 1) return std::exp(entropy_index(x));
  double m = 0;
  for (double s : shares(x))
    if (s > 0) m += std::pow(s, alpha);
  return std::pow(m, 1 / (1 - alpha));
}

}  // namespace nomos::measurement
