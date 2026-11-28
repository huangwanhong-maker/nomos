// Progress toward a declared target: the trend to target of a series observed at a baseline and at a latest date (decision 0019).
//
// - The compound annual growth rate between a baseline year b and the latest year t with data,
//     r = (x_t / x_b)^{1/(t − b)} − 1,
//   defined for x_b > 0 and x_t ≥ 0.
// - With a target x* for a target year T, the value expected at T if the rate holds, x̂_T = x_t (1 + r)^{T − t}, and the share of the
//   distance from the baseline to the target that the trend would cover,
//     ρ = (x̂_T − x_b) / (x* − x_b),
//   classed by declared floors on ρ (the scheme). ρ is the same whichever way the target lies, so one rule serves a quantity to be
//   raised and one to be lowered; ρ < 0 is movement away from the target. Where x_b = x* the ratio is undefined and the trend is the
//   best class if the latest value is at or beyond the target and the worst otherwise; a latest value at or beyond the target is the
//   best class (the target met), whatever the trend; a baseline already beyond the target with a latest value that is not is the worst
//   class (fallen back from the target), where the ratio's sign would otherwise read the reversal as progress.
// - Without a target, the rate itself is classed by declared floors, its sign turned for a quantity to be lowered.
// - The baseline is the value at a declared base year if observed, else the earliest observed after it; a series with no value at or
//   after the base year, or with fewer than two observed years from it, is not assessed.
//
// The scheme is a declaration. Its default is the one the United Nations Statistics Division states for the global assessment of the
// Sustainable Development Goals (floors 0.95, 0.50, 0.10 and 0 on ρ; 2, 1.25, 0.5 and 0 per cent a year on r; five classes from
// "on track or target met" to "regression"; base year 2015): UN DESA Statistics Division, *Technical Note for Progress Assessment
// 2025*, Section 3.1 (https://unstats.un.org/sdgs/files/report/2025/Technical_Note_for_Progress_Assessment_2025.pdf). The note's
// second and third steps (scores from 0 to 5 and their means by target) are not computed here.
//
// *Caveat.* A compound rate assumes proportional change; a quantity that falls by a fixed amount a year is assessed as slowing. The
// trend is read from two points and ignores everything between them, as the note's method does.

#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::measurement {

enum class Desired { increase, decrease };

struct ProgressScheme {
  std::vector<double> ratio_floors{0.95, 0.50, 0.10, 0.0};      // class k when ρ ≥ floors[k]; below the last, the last class
  std::vector<double> rate_floors{0.02, 0.0125, 0.005, 0.0};    // class k when r ≥ floors[k] (r negated for a decrease)
  std::vector<std::string> names{"on track or target met", "moderate progress, but acceleration needed",
                                 "marginal progress, and significant acceleration needed", "stagnation", "regression"};
  double base_year = 2015;
};

struct Progress {
  double baseline_year = 0, latest_year = 0, baseline = 0, latest = 0;
  double rate = 0;                                                   // compound annual growth rate r
  bool has_target = false;
  double target = 0, target_year = 0;
  double expected = std::numeric_limits<double>::quiet_NaN();       // x̂_T
  double ratio = std::numeric_limits<double>::quiet_NaN();          // ρ, where defined
  bool met = false;                                                  // the latest value at or beyond the target
  std::size_t level = 0;                                             // 0 the best class
  std::string trend;                                                 // its name
  std::string basis;                                                 // "ratio", "rate", "target met", "baseline at target"
};

namespace progress_detail {
inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument(what);
}
inline void check(const ProgressScheme& s) {
  require(s.ratio_floors.size() + 1 == s.names.size() && s.rate_floors.size() + 1 == s.names.size(),
          "a scheme names one class more than it has floors");
  for (std::size_t k = 1; k < s.ratio_floors.size(); k++) require(s.ratio_floors[k] < s.ratio_floors[k - 1], "the floors on the ratio descend");
  for (std::size_t k = 1; k < s.rate_floors.size(); k++) require(s.rate_floors[k] < s.rate_floors[k - 1], "the floors on the rate descend");
}
inline std::size_t classed(double value, const std::vector<double>& floors) {
  for (std::size_t k = 0; k < floors.size(); k++)
    if (value >= floors[k]) return k;
  return floors.size();
}
}  // namespace progress_detail

inline double compound_rate(double baseline, double latest, double years) {
  progress_detail::require(baseline > 0 && latest >= 0 && std::isfinite(baseline) && std::isfinite(latest),
                           "a compound rate needs a positive baseline and a latest value not below zero");
  progress_detail::require(years > 0, "a compound rate needs the latest year after the baseline");
  return std::pow(latest / baseline, 1.0 / years) - 1.0;
}

inline std::size_t level_of_ratio(double ratio, const ProgressScheme& s = {}) {
  progress_detail::check(s);
  return progress_detail::classed(ratio, s.ratio_floors);
}

inline std::size_t level_of_rate(double rate, Desired desired, const ProgressScheme& s = {}) {
  progress_detail::check(s);
  return progress_detail::classed(desired == Desired::increase ? rate : -rate, s.rate_floors);
}

// The trend of a series from (b, x_b) to (t, x_t): toward a target x* at T when one is given, else by its rate.
inline Progress trend(double b, double xb, double t, double xt, Desired desired, const ProgressScheme& s = {}, bool has_target = false,
                      double target = 0, double target_year = 2030) {
  progress_detail::check(s);
  Progress p;
  p.baseline_year = b, p.latest_year = t, p.baseline = xb, p.latest = xt;
  p.rate = compound_rate(xb, xt, t - b);
  p.has_target = has_target;
  if (!has_target) {
    p.level = level_of_rate(p.rate, desired, s);
    p.basis = "rate";
  } else {
    progress_detail::require(std::isfinite(target) && target_year >= t, "a target is a finite value for a year not before the latest");
    p.target = target, p.target_year = target_year;
    p.expected = xt * std::pow(1.0 + p.rate, target_year - t);
    p.met = desired == Desired::increase ? xt >= target : xt <= target;
    const bool base_met = desired == Desired::increase ? xb >= target : xb <= target;
    if (xb == target) {
      p.level = p.met ? 0 : s.names.size() - 1;
      p.basis = "baseline at target";
    } else if (base_met) {  // the target was met at the baseline: met still, or fallen back from it (the ratio's sign would mislead)
      p.level = p.met ? 0 : s.names.size() - 1;
      p.basis = p.met ? "target met" : "fallen back from the target";
    } else {
      p.ratio = (p.expected - xb) / (target - xb);
      p.level = p.met ? 0 : progress_detail::classed(p.ratio, s.ratio_floors);
      p.basis = p.met ? "target met" : "ratio";
    }
  }
  p.trend = s.names[p.level];
  return p;
}

// The baseline year of a series: the base year if observed, else the earliest observed year after it. Years need not be sorted;
// values that are not finite are taken as not observed.
inline std::size_t baseline_index(const std::vector<double>& years, const std::vector<double>& values, double base_year) {
  progress_detail::require(years.size() == values.size(), "a series gives one value for each year");
  std::size_t best = years.size();
  for (std::size_t i = 0; i < years.size(); i++)
    if (std::isfinite(values[i]) && years[i] >= base_year && (best == years.size() || years[i] < years[best])) best = i;
  progress_detail::require(best < years.size(), "the series has no value at or after the base year");
  return best;
}

inline std::size_t latest_index(const std::vector<double>& years, const std::vector<double>& values) {
  std::size_t best = years.size();
  for (std::size_t i = 0; i < years.size(); i++)
    if (std::isfinite(values[i]) && (best == years.size() || years[i] > years[best])) best = i;
  progress_detail::require(best < years.size(), "the series has no value");
  return best;
}

// The trend of an observed series from its baseline (the scheme's base year) to its latest year.
inline Progress trend_of_series(const std::vector<double>& years, const std::vector<double>& values, Desired desired, const ProgressScheme& s = {},
                                bool has_target = false, double target = 0, double target_year = 2030) {
  const std::size_t b = baseline_index(years, values, s.base_year), t = latest_index(years, values);
  progress_detail::require(years[t] > years[b], "the series has fewer than two observed years from its baseline");
  return trend(years[b], values[b], years[t], values[t], desired, s, has_target, target, target_year);
}

}  // namespace nomos::measurement
