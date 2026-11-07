// Composite indicators: several factors, each on a scale of its own, brought to one scale and combined by weights
// (OECD and JRC 2008, Handbook on Constructing Composite Indicators: normalisation, §5; weighting and aggregation, §6). A factor
// is scaled between two of its percentiles (so that an outlying unit does not wash the rest into one end, a robust
// min-max), inverted if more of it is worse, and a unit's composite is the weighted arithmetic mean of the factors it has;
// a factor missing for a unit is left out of its mean and the weights renormalised, as a linear aggregation allows.
//
// Caveat. A linear aggregation is compensatory: much of one factor makes up for little of another (Handbook, §6.2). What the
// weights mean is the analyst's to declare, and the composite says nothing a sensitivity analysis to them has not tried.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nomos::measurement {

struct Factor {
  std::vector<double> values;   // one per unit; NaN where the factor is missing
  double weight = 1;            // non-negative
  bool invert = false;          // more of it is worse
};

// A factor scaled to [0, top] between its low and high percentiles (linear interpolation between order statistics).
inline std::vector<double> scaled(const Factor& f, double low = 0.02, double high = 0.98, double top = 1.0) {
  if (!(low >= 0 && low < high && high <= 1)) throw std::invalid_argument("composite: percentiles 0 <= low < high <= 1");
  std::vector<double> known;
  known.reserve(f.values.size());
  for (const double v : f.values) if (std::isfinite(v)) known.push_back(v);
  std::vector<double> out(f.values.size(), std::numeric_limits<double>::quiet_NaN());
  if (known.empty()) return out;
  std::sort(known.begin(), known.end());
  const auto at = [&](double q) {
    const double h = (static_cast<double>(known.size()) - 1) * q;
    const std::size_t i = static_cast<std::size_t>(h);
    return i + 1 < known.size() ? known[i] + (h - static_cast<double>(i)) * (known[i + 1] - known[i]) : known[i];
  };
  const double lo = at(low), hi = at(high), span = hi - lo;
  for (std::size_t k = 0; k < out.size(); k++) {
    if (!std::isfinite(f.values[k])) continue;
    double u = span > 0 ? std::clamp((f.values[k] - lo) / span, 0.0, 1.0) : 0.5;
    if (f.invert) u = 1 - u;
    out[k] = top * u;
  }
  return out;
}

// The composite of the factors, on [0, top]: NaN for a unit none of whose factors is known.
inline std::vector<double> composite(const std::vector<Factor>& factors, double low = 0.02, double high = 0.98, double top = 5.0) {
  if (factors.empty()) throw std::invalid_argument("composite: no factors");
  const std::size_t n = factors.front().values.size();
  std::vector<double> sum(n, 0.0), weight(n, 0.0);
  for (const Factor& f : factors) {
    if (f.values.size() != n) throw std::invalid_argument("composite: factors of different lengths");
    if (f.weight < 0) throw std::invalid_argument("composite: a weight is not negative");
    const auto s = scaled(f, low, high, top);
    for (std::size_t k = 0; k < n; k++) if (std::isfinite(s[k])) { sum[k] += f.weight * s[k]; weight[k] += f.weight; }
  }
  std::vector<double> out(n, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t k = 0; k < n; k++) if (weight[k] > 0) out[k] = sum[k] / weight[k];
  return out;
}

}  // namespace nomos::measurement
