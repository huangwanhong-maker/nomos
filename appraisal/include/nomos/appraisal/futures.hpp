// Methods of futures studies that compute: a Delphi panel's rounds, and the space of configurations a morphological field allows
// (decision 0012; the handbook's §10: Delphi, morphological analysis).
//
//   delphi(rounds, o)        a panel's estimates of one quantity, round by round: each round's median, interquartile range (Hyndman
//                            and Fan's type 7), and trimmed mean; consensus when the interquartile range relative to the median falls
//                            to a declared bound; stability when the median moves less than a declared share between rounds (von der
//                            Gracht 2012, Technol. Forecast. Soc. Change 79: 1525–1536, reviews these criteria; neither is a proof of
//                            truth, only of agreement)
//   morphological(field)     Zwicky's morphological box: parameters each with values, a configuration choosing one value of each;
//                            Ritchey's cross-consistency assessment (2011, Wicked Problems – Social Messes, ch. 2) declares pairs of
//                            values that cannot hold together, and the field is reduced to the configurations with no such pair, found
//                            by depth-first search with each pair checked as soon as both values are chosen
//
// Neither method is a forecast: the first records what a panel agrees on, the second what a declared set of exclusions leaves open.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nomos::appraisal {

struct DelphiOptions {
  double consensus = 0.25;   // the interquartile range relative to |median| at which the panel agrees
  double stability = 0.05;   // the relative move of the median between rounds below which it has settled
  double trim = 0.1;         // the share trimmed from each end for the trimmed mean
};

struct DelphiRound {
  double median = 0, q1 = 0, q3 = 0, iqr = 0, trimmed_mean = 0, relative_iqr = 0;
  std::size_t panel = 0;
  bool consensus = false, stable = false;  // stable: compared with the round before
};

namespace detail {
inline double type7(const std::vector<double>& sorted, double p) {
  const double pos = p * static_cast<double>(sorted.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double f = pos - static_cast<double>(i);
  return i + 1 < sorted.size() ? sorted[i] * (1 - f) + sorted[i + 1] * f : sorted.back();
}
}  // namespace detail

inline std::vector<DelphiRound> delphi(const std::vector<std::vector<double>>& rounds, const DelphiOptions& o = {}) {
  std::vector<DelphiRound> out;
  for (const auto& r : rounds) {
    std::vector<double> v;
    for (double x : r)
      if (std::isfinite(x)) v.push_back(x);
    if (v.empty()) throw std::invalid_argument("nomos::appraisal::delphi: a round with at least one estimate");
    std::sort(v.begin(), v.end());
    DelphiRound d;
    d.panel = v.size();
    d.median = detail::type7(v, 0.5), d.q1 = detail::type7(v, 0.25), d.q3 = detail::type7(v, 0.75), d.iqr = d.q3 - d.q1;
    const std::size_t cut = static_cast<std::size_t>(std::floor(o.trim * static_cast<double>(v.size())));
    double s = 0;
    for (std::size_t i = cut; i < v.size() - cut; i++) s += v[i];
    d.trimmed_mean = s / static_cast<double>(v.size() - 2 * cut);
    d.relative_iqr = d.median != 0 ? d.iqr / std::fabs(d.median) : (d.iqr == 0 ? 0.0 : INFINITY);
    d.consensus = d.relative_iqr <= o.consensus;
    if (!out.empty()) {
      const double prev = out.back().median;
      d.stable = prev != 0 ? std::fabs(d.median - prev) / std::fabs(prev) <= o.stability : d.median == 0;
    }
    out.push_back(d);
  }
  return out;
}

struct Field {
  std::vector<std::string> parameters;
  std::vector<std::vector<std::string>> values;  // [parameter][value]
  // Pairs that cannot hold together: ((parameter, value), (parameter, value)).
  std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>>> excluded;
};

struct Morphology {
  std::size_t size = 0;                              // the configurations of the whole box
  std::size_t count = 0;                             // the configurations with no excluded pair
  std::vector<std::vector<std::size_t>> consistent;  // the first `keep` of them (each: a value index per parameter)
  std::size_t visited = 0;                           // the partial configurations the search visited
};

inline Morphology morphological(const Field& f, std::size_t keep = 100000) {
  const std::size_t P = f.parameters.size();
  if (f.values.size() != P || P == 0) throw std::invalid_argument("nomos::appraisal::morphological: values for each parameter");
  std::set<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>>> bad;
  for (auto e : f.excluded) {
    if (e.first.first >= P || e.second.first >= P || e.first.second >= f.values[e.first.first].size() || e.second.second >= f.values[e.second.first].size())
      throw std::invalid_argument("nomos::appraisal::morphological: an exclusion names a value that does not exist");
    if (e.first.first == e.second.first) throw std::invalid_argument("nomos::appraisal::morphological: an exclusion pairs two parameters");
    if (e.first > e.second) std::swap(e.first, e.second);
    bad.insert(e);
  }
  Morphology m;
  m.size = 1;
  for (const auto& v : f.values) m.size *= v.size();
  std::vector<std::size_t> pick;
  auto clash = [&](std::size_t p, std::size_t v) {
    for (std::size_t q = 0; q < pick.size(); q++) {
      std::pair<std::size_t, std::size_t> a{q, pick[q]}, b{p, v};
      if (a > b) std::swap(a, b);
      if (bad.count({a, b})) return true;
    }
    return false;
  };
  auto search = [&](auto&& self) -> void {
    m.visited++;
    if (pick.size() == P) {
      m.count++;
      if (m.consistent.size() < keep) m.consistent.push_back(pick);
      return;
    }
    const std::size_t p = pick.size();
    for (std::size_t v = 0; v < f.values[p].size(); v++) {
      if (clash(p, v)) continue;
      pick.push_back(v);
      self(self);
      pick.pop_back();
    }
  };
  search(search);
  return m;
}

}  // namespace nomos::appraisal
