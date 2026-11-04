// Poverty: the Foster–Greer–Thorbecke family and the Watts index.
//
// A unit is poor when its amount lies strictly below the line z > 0; g_i = (z − x_i)/z is its
// normalised shortfall. The Foster–Greer–Thorbecke family (Foster, Greer and Thorbecke 1984;
// with weights, as surveys use it, Haughton and Khandker 2009, ch. 4) is
//
//   FGT(α) = (1/W) Σ_{x_i < z} w_i g_i^α,   α ≥ 0:
//
// α = 0 is the headcount ratio (the share of the population that is poor), α = 1 the poverty gap
// index (the mean shortfall, the nonpoor counting zero), α = 2 the severity index (the mean
// squared shortfall, which weights the poorest most). For α > 1 the index satisfies the transfer
// axiom (Sen 1976): a transfer from a poor unit to a poorer one lowers it. Every FGT(α) is
// additively decomposable over a partition of the population, FGT = Σ_g (W_g / W) FGT_g (Foster,
// Greer and Thorbecke 1984), which decomposition.hpp computes.
//
// The Watts index (Watts 1968; Zheng 1993) is (1/W) Σ_{x_i < z} w_i ln(z / x_i), the first
// distribution-sensitive poverty index; it needs the poor to have positive amounts.
//
// Caveat. "Strictly below" matters only for α = 0, where a unit exactly at the line is not
// counted; it is the common convention (Haughton and Khandker 2009, ch. 4).

#pragma once

#include <cmath>
#include <span>
#include <vector>

#include "inequality.hpp"

namespace nomos::measurement {

// FGT(α) at the line z > 0, α ≥ 0.
inline double fgt(std::span<const double> x, double z, double alpha, std::span<const double> w = {}) {
  detail::require(z > 0 && std::isfinite(z), "a positive line z");
  detail::require(alpha >= 0 && std::isfinite(alpha), "an exponent α ≥ 0");
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++)
    if (x[i] < z) s += ww[i] * (alpha == 0 ? 1.0 : std::pow((z - x[i]) / z, alpha));
  return s / t.weight;
}

// The headcount ratio FGT(0), the poverty gap index FGT(1), the severity index FGT(2).
inline double headcount_ratio(std::span<const double> x, double z, std::span<const double> w = {}) { return fgt(x, z, 0, w); }
inline double poverty_gap_index(std::span<const double> x, double z, std::span<const double> w = {}) { return fgt(x, z, 1, w); }
inline double poverty_severity(std::span<const double> x, double z, std::span<const double> w = {}) { return fgt(x, z, 2, w); }

// The Watts index (1/W) Σ_{x < z} w ln(z / x).
inline double watts_index(std::span<const double> x, double z, std::span<const double> w = {}) {
  detail::require(z > 0 && std::isfinite(z), "a positive line z");
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  double s = 0;
  for (std::size_t i = 0; i < x.size(); i++)
    if (x[i] < z && ww[i] > 0) {
      detail::require(x[i] > 0, "positive amounts below the line (Watts)");
      s += ww[i] * std::log(z / x[i]);
    }
  return s / t.weight;
}

}  // namespace nomos::measurement
