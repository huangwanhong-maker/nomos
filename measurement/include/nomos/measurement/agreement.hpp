// Agreement among coders: how far several people who classify the same units agree beyond what chance would give, the reliability
// on which coded qualitative data rest (decision 0012; the handbook's §8: interviews and coded material).
//
//   cohen(a, b, categories, weight)   Cohen's κ for two coders, (p_o − p_e)/(1 − p_e) with p_e from the coders' own margins (Cohen
//                                      1960, Educ. Psychol. Meas. 20: 37–46); weighted κ for ordered categories with linear or quadratic
//                                      disagreement weights (Cohen 1968, Psychol. Bull. 70: 213–220)
//   fleiss(counts)                     Fleiss's κ for many coders, each unit classified by the same number of them: counts[u][c] the
//                                      coders who put unit u in category c (Fleiss 1971, Psychol. Bull. 76: 378–382)
//   krippendorff(data, metric)         Krippendorff's α = 1 − D_o/D_e for any number of coders, missing values allowed, by the coincidence
//                                      matrix, with nominal, ordinal, interval or ratio differences (Krippendorff 2011, Computing
//                                      Krippendorff's Alpha-Reliability; Hayes and Krippendorff 2007, Commun. Methods Meas. 1: 77–89);
//                                      units coded by fewer than two coders are left out, as the coefficient requires
//
// Codes are integers 0 .. K − 1 (an unrecorded code is −1); for the ordinal, interval and ratio metrics the integers are the values.

#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::measurement {

enum class Disagreement { nominal, linear, quadratic };

inline double cohen(const std::vector<int>& a, const std::vector<int>& b, std::size_t K, Disagreement w = Disagreement::nominal) {
  if (a.size() != b.size() || a.empty() || K < 2) throw std::invalid_argument("nomos::measurement::cohen: two coders' codes of the same units, two or more categories");
  std::vector<std::vector<double>> p(K, std::vector<double>(K, 0.0));
  double n = 0;
  for (std::size_t u = 0; u < a.size(); u++) {
    if (a[u] < 0 || b[u] < 0) continue;
    if (static_cast<std::size_t>(a[u]) >= K || static_cast<std::size_t>(b[u]) >= K) throw std::invalid_argument("nomos::measurement::cohen: a code outside 0 .. K − 1");
    p[static_cast<std::size_t>(a[u])][static_cast<std::size_t>(b[u])] += 1, n += 1;
  }
  if (n == 0) throw std::invalid_argument("nomos::measurement::cohen: no unit coded by both");
  std::vector<double> ra(K, 0.0), rb(K, 0.0);
  for (std::size_t i = 0; i < K; i++)
    for (std::size_t j = 0; j < K; j++) p[i][j] /= n, ra[i] += p[i][j], rb[j] += p[i][j];
  auto weight = [&](std::size_t i, std::size_t j) {
    const double d = std::fabs(static_cast<double>(i) - static_cast<double>(j)) / static_cast<double>(K - 1);
    return w == Disagreement::nominal ? (i == j ? 0.0 : 1.0) : w == Disagreement::linear ? d : d * d;
  };
  double obs = 0, exp = 0;
  for (std::size_t i = 0; i < K; i++)
    for (std::size_t j = 0; j < K; j++) obs += weight(i, j) * p[i][j], exp += weight(i, j) * ra[i] * rb[j];
  return exp > 0 ? 1 - obs / exp : NAN;
}

inline double fleiss(const std::vector<std::vector<double>>& counts) {
  if (counts.empty()) throw std::invalid_argument("nomos::measurement::fleiss: units");
  const std::size_t K = counts[0].size();
  double raters = -1, N = static_cast<double>(counts.size());
  std::vector<double> pj(K, 0.0);
  double Pbar = 0;
  for (const auto& row : counts) {
    if (row.size() != K) throw std::invalid_argument("nomos::measurement::fleiss: one count per category for every unit");
    double n = 0, sq = 0;
    for (std::size_t j = 0; j < K; j++) n += row[j], sq += row[j] * row[j], pj[j] += row[j];
    if (raters < 0) raters = n;
    if (n != raters || n < 2) throw std::invalid_argument("nomos::measurement::fleiss: every unit coded by the same number of coders, at least two");
    Pbar += (sq - n) / (n * (n - 1));
  }
  Pbar /= N;
  double Pe = 0;
  for (double& v : pj) v /= N * raters, Pe += v * v;
  return Pe < 1 ? (Pbar - Pe) / (1 - Pe) : NAN;
}

enum class Metric { nominal, ordinal, interval, ratio };

// data[coder][unit], −1 where a coder did not code the unit.
inline double krippendorff(const std::vector<std::vector<int>>& data, std::size_t K, Metric metric = Metric::nominal) {
  if (data.size() < 2) throw std::invalid_argument("nomos::measurement::krippendorff: at least two coders");
  const std::size_t U = data[0].size();
  for (const auto& c : data)
    if (c.size() != U) throw std::invalid_argument("nomos::measurement::krippendorff: every coder's row covers the same units");
  // The coincidence matrix: within each unit coded m_u ≥ 2 times, every ordered pair of its values from different coders, weighted 1/(m_u − 1).
  std::vector<std::vector<double>> o(K, std::vector<double>(K, 0.0));
  for (std::size_t u = 0; u < U; u++) {
    std::vector<std::size_t> v;
    for (const auto& c : data)
      if (c[u] >= 0) {
        if (static_cast<std::size_t>(c[u]) >= K) throw std::invalid_argument("nomos::measurement::krippendorff: a code outside 0 .. K − 1");
        v.push_back(static_cast<std::size_t>(c[u]));
      }
    if (v.size() < 2) continue;
    const double w = 1.0 / static_cast<double>(v.size() - 1);
    for (std::size_t i = 0; i < v.size(); i++)
      for (std::size_t j = 0; j < v.size(); j++)
        if (i != j) o[v[i]][v[j]] += w;
  }
  std::vector<double> nc(K, 0.0);
  double n = 0;
  for (std::size_t c = 0; c < K; c++)
    for (std::size_t k = 0; k < K; k++) nc[c] += o[c][k], n += o[c][k];
  if (n < 2) throw std::invalid_argument("nomos::measurement::krippendorff: fewer than two pairable values");
  auto delta2 = [&](std::size_t c, std::size_t k) -> double {
    if (c == k) return 0;
    switch (metric) {
      case Metric::nominal: return 1;
      case Metric::interval: return (static_cast<double>(c) - static_cast<double>(k)) * (static_cast<double>(c) - static_cast<double>(k));
      case Metric::ratio: {
        const double a = static_cast<double>(c), b = static_cast<double>(k);
        return a + b > 0 ? ((a - b) / (a + b)) * ((a - b) / (a + b)) : 0.0;
      }
      case Metric::ordinal: {
        const std::size_t lo = std::min(c, k), hi = std::max(c, k);
        double s = 0;
        for (std::size_t g = lo; g <= hi; g++) s += nc[g];
        s -= (nc[c] + nc[k]) / 2;
        return s * s;
      }
    }
    return 0;
  };
  double Do = 0, De = 0;
  for (std::size_t c = 0; c < K; c++)
    for (std::size_t k = 0; k < K; k++) Do += o[c][k] * delta2(c, k), De += nc[c] * nc[k] * delta2(c, k);
  Do /= n, De /= n * (n - 1);
  return De > 0 ? 1 - Do / De : NAN;
}

}  // namespace nomos::measurement
