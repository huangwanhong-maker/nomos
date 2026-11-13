// Multi-criteria decision analysis: alternatives compared on several criteria, each on a scale of its own, by declared weights;
// and how far the weights can move before the answer changes (decision 0012).
//
// A decision matrix holds alternatives × criteria, each criterion to be maximised (a benefit) or minimised (a cost), with weights
// w ≥ 0 summing to one. What the weights mean is the analyst's to declare; the methods differ in what they presume of them:
//
//   screen(m, limits)          non-compensatory constraints applied first: an alternative that breaks a declared threshold on any
//                              criterion is set aside, whatever its other scores (no score buys back a right; handbook §9.1)
//   weighted_sum(m)            min–max normalisation to [0, 1] and the weighted sum: fully compensatory
//   topsis(m)                  closeness to the ideal and distance from the anti-ideal after vector normalisation (Hwang and Yoon
//                              1981, Multiple Attribute Decision Making, ch. 3; Behzadian et al. 2012, Expert Syst. Appl. 39: 13051)
//   ahp(P)                     weights from a reciprocal matrix of pairwise judgements: its Perron vector (calx), λ_max, the
//                              consistency index CI = (λ_max − n)/(n − 1) and ratio CR = CI/RI with Saaty's random indices (Saaty 1980,
//                              The Analytic Hierarchy Process; RI as Saaty 1987 tabulates); CR above 0.1 is conventionally too
//                              inconsistent to use
//   promethee(m, prefs)        outranking by preference functions of the difference on each criterion (usual, U-shape, V-shape,
//                              linear, Gaussian), the positive, negative and net flows of PROMETHEE II (Brans and Vincke 1985,
//                              Manage. Sci. 31: 647–656): partially compensatory, with indifference and preference thresholds declared
//   weight_stability(m)        for each criterion, the interval of its weight (the others scaled in proportion) over which the
//                              weighted sum's best alternative stays best, found exactly where the lines of the alternatives cross
//   acceptability(m, draws)    SMAA-2's rank acceptability (Lahdelma and Salminen 2001, Oper. Res. 49: 444–454): the share of weights
//                              drawn uniformly on the simplex under which each alternative takes each rank, when no weight is declared
//
// Accuracy: exact arithmetic but for the AHP's eigenvector (inverse iteration, to rounding) and SMAA's Monte Carlo (its standard
// error √(p(1 − p)/draws)).

#pragma once

#include <calx/linalg/dense.hpp>
#include <calx/linalg/nonnegative.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::appraisal {

struct Matrix {
  std::vector<std::string> alternatives, criteria;
  std::vector<std::vector<double>> score;  // [alternative][criterion]
  std::vector<bool> benefit;               // true: more is better
  std::vector<double> weight;              // ≥ 0, summing to one (normalised on use)
};

namespace detail {
inline void require(bool ok, const std::string& what) {
  if (!ok) throw std::invalid_argument("nomos::appraisal: " + what);
}
inline void validate(const Matrix& m) {
  const std::size_t n = m.alternatives.size(), k = m.criteria.size();
  require(n >= 1 && k >= 1 && m.score.size() == n && m.benefit.size() == k, "a matrix of alternatives × criteria, each criterion a benefit or a cost");
  for (const auto& row : m.score) require(row.size() == k, "one score per criterion");
}
inline std::vector<double> weights(const Matrix& m) {
  const std::size_t k = m.criteria.size();
  std::vector<double> w = m.weight.empty() ? std::vector<double>(k, 1.0 / static_cast<double>(k)) : m.weight;
  require(w.size() == k, "one weight per criterion");
  double s = 0;
  for (double v : w) require(v >= 0, "weights ≥ 0"), s += v;
  require(s > 0, "some weight > 0");
  for (double& v : w) v /= s;
  return w;
}
inline std::vector<std::vector<double>> minmax(const Matrix& m) {
  const std::size_t n = m.alternatives.size(), k = m.criteria.size();
  std::vector<std::vector<double>> r(n, std::vector<double>(k, 0.0));
  for (std::size_t j = 0; j < k; j++) {
    double lo = INFINITY, hi = -INFINITY;
    for (std::size_t i = 0; i < n; i++) lo = std::min(lo, m.score[i][j]), hi = std::max(hi, m.score[i][j]);
    for (std::size_t i = 0; i < n; i++) r[i][j] = hi > lo ? (m.benefit[j] ? m.score[i][j] - lo : hi - m.score[i][j]) / (hi - lo) : 1.0;
  }
  return r;
}
}  // namespace detail

// ---- non-compensatory screening ----

struct Limit {
  std::size_t criterion = 0;
  double bound = 0;
  bool at_least = true;  // the score must be ≥ bound (else ≤ bound)
};

struct Screened {
  Matrix kept;
  std::vector<std::string> set_aside;
  std::vector<std::string> why;  // the first limit each breaks
};

inline Screened screen(const Matrix& m, const std::vector<Limit>& limits) {
  detail::validate(m);
  Screened s;
  s.kept = m;
  s.kept.alternatives.clear(), s.kept.score.clear();
  for (std::size_t i = 0; i < m.alternatives.size(); i++) {
    std::string broken;
    for (const auto& l : limits) {
      detail::require(l.criterion < m.criteria.size(), "a limit on a criterion that exists");
      const double v = m.score[i][l.criterion];
      if (l.at_least ? v < l.bound : v > l.bound) {
        broken = m.criteria[l.criterion] + (l.at_least ? " below " : " above ") + std::to_string(l.bound);
        break;
      }
    }
    if (broken.empty()) s.kept.alternatives.push_back(m.alternatives[i]), s.kept.score.push_back(m.score[i]);
    else s.set_aside.push_back(m.alternatives[i]), s.why.push_back(broken);
  }
  return s;
}

// ---- weighted sum and TOPSIS ----

struct Ranking {
  std::vector<double> value;     // per alternative, higher is better
  std::vector<std::size_t> rank; // 1 = best
  std::vector<std::size_t> order;
};

namespace detail {
inline Ranking ranking(std::vector<double> v) {
  Ranking r;
  r.value = std::move(v);
  r.order.resize(r.value.size());
  std::iota(r.order.begin(), r.order.end(), 0);
  std::stable_sort(r.order.begin(), r.order.end(), [&](std::size_t a, std::size_t b) { return r.value[a] > r.value[b]; });
  r.rank.resize(r.value.size());
  for (std::size_t p = 0; p < r.order.size(); p++) r.rank[r.order[p]] = p + 1;
  return r;
}
}  // namespace detail

inline Ranking weighted_sum(const Matrix& m) {
  detail::validate(m);
  const auto w = detail::weights(m);
  const auto r = detail::minmax(m);
  std::vector<double> v(m.alternatives.size(), 0.0);
  for (std::size_t i = 0; i < v.size(); i++)
    for (std::size_t j = 0; j < w.size(); j++) v[i] += w[j] * r[i][j];
  return detail::ranking(v);
}

inline Ranking topsis(const Matrix& m) {
  detail::validate(m);
  const auto w = detail::weights(m);
  const std::size_t n = m.alternatives.size(), k = w.size();
  std::vector<std::vector<double>> v(n, std::vector<double>(k));
  for (std::size_t j = 0; j < k; j++) {
    double norm = 0;
    for (std::size_t i = 0; i < n; i++) norm += m.score[i][j] * m.score[i][j];
    norm = std::sqrt(norm);
    for (std::size_t i = 0; i < n; i++) v[i][j] = norm > 0 ? w[j] * m.score[i][j] / norm : 0.0;
  }
  std::vector<double> best(k), worst(k);
  for (std::size_t j = 0; j < k; j++) {
    double hi = -INFINITY, lo = INFINITY;
    for (std::size_t i = 0; i < n; i++) hi = std::max(hi, v[i][j]), lo = std::min(lo, v[i][j]);
    best[j] = m.benefit[j] ? hi : lo, worst[j] = m.benefit[j] ? lo : hi;
  }
  std::vector<double> c(n);
  for (std::size_t i = 0; i < n; i++) {
    double dp = 0, dm = 0;
    for (std::size_t j = 0; j < k; j++) dp += (v[i][j] - best[j]) * (v[i][j] - best[j]), dm += (v[i][j] - worst[j]) * (v[i][j] - worst[j]);
    dp = std::sqrt(dp), dm = std::sqrt(dm);
    c[i] = dp + dm > 0 ? dm / (dp + dm) : 0.5;
  }
  return detail::ranking(c);
}

// ---- the analytic hierarchy process ----

struct Ahp {
  std::vector<double> weight;  // the Perron vector, summing to one
  double lambda = 0, ci = 0, cr = 0;
  bool consistent = false;     // CR ≤ 0.1
};

inline Ahp ahp(const std::vector<std::vector<double>>& P) {
  const std::size_t n = P.size();
  detail::require(n >= 2, "ahp: at least two items compared");
  calx::linalg::Matrix A(n);
  for (std::size_t i = 0; i < n; i++) {
    detail::require(P[i].size() == n, "ahp: a square matrix of judgements");
    for (std::size_t j = 0; j < n; j++) {
      detail::require(P[i][j] > 0, "ahp: positive judgements");
      detail::require(std::fabs(P[i][j] * P[j][i] - 1) < 1e-9, "ahp: a reciprocal matrix, a_ji = 1/a_ij");
      A(i, j) = P[i][j];
    }
  }
  const auto pf = calx::linalg::perron(A);
  Ahp r;
  r.weight = pf.right;
  r.lambda = pf.lambda;
  r.ci = n > 2 ? (pf.lambda - static_cast<double>(n)) / static_cast<double>(n - 1) : 0.0;
  static const double RI[] = {0, 0, 0, 0.58, 0.90, 1.12, 1.24, 1.32, 1.41, 1.45, 1.49, 1.51, 1.48, 1.56, 1.57, 1.59};  // Saaty 1987
  r.cr = n > 2 && n < 16 ? r.ci / RI[n] : 0.0;
  r.consistent = r.cr <= 0.1;
  return r;
}

// ---- PROMETHEE II ----

struct Preference {
  enum Kind { usual, u_shape, v_shape, linear, gaussian } kind = usual;
  double q = 0, p = 0, s = 0;  // indifference, preference, Gaussian thresholds
  double operator()(double d) const {  // d: the advantage of a over b on the criterion, oriented so that more is better
    if (d <= 0) return 0;
    switch (kind) {
      case usual: return 1;
      case u_shape: return d > q ? 1 : 0;
      case v_shape: return d >= p ? 1 : d / p;
      case linear: return d <= q ? 0 : d >= p ? 1 : (d - q) / (p - q);
      case gaussian: return 1 - std::exp(-d * d / (2 * s * s));
    }
    return 0;
  }
};

struct Flows {
  std::vector<double> positive, negative, net;
  Ranking ranking;
};

inline Flows promethee(const Matrix& m, const std::vector<Preference>& prefs) {
  detail::validate(m);
  const auto w = detail::weights(m);
  const std::size_t n = m.alternatives.size(), k = w.size();
  detail::require(prefs.size() == k, "promethee: one preference function per criterion");
  Flows f;
  f.positive.assign(n, 0.0), f.negative.assign(n, 0.0), f.net.assign(n, 0.0);
  if (n < 2) {
    f.ranking = detail::ranking(f.net);
    return f;
  }
  for (std::size_t a = 0; a < n; a++)
    for (std::size_t b = 0; b < n; b++) {
      if (a == b) continue;
      double pi = 0;
      for (std::size_t j = 0; j < k; j++) {
        const double d = m.benefit[j] ? m.score[a][j] - m.score[b][j] : m.score[b][j] - m.score[a][j];
        pi += w[j] * prefs[j](d);
      }
      f.positive[a] += pi / static_cast<double>(n - 1);
      f.negative[b] += pi / static_cast<double>(n - 1);
    }
  for (std::size_t a = 0; a < n; a++) f.net[a] = f.positive[a] - f.negative[a];
  f.ranking = detail::ranking(f.net);
  return f;
}

// ---- how far the weights can move ----

struct Stability {
  std::vector<double> lower, upper;  // per criterion: the interval of its weight over which the best alternative stays best
  std::size_t best = 0;
};

inline Stability weight_stability(const Matrix& m) {
  detail::validate(m);
  const auto w = detail::weights(m);
  const auto r = detail::minmax(m);
  const std::size_t n = m.alternatives.size(), k = w.size();
  Stability st;
  st.best = weighted_sum(m).order.front();
  for (std::size_t j = 0; j < k; j++) {
    // With w_j = t and the others w_l (1 − t)/(1 − w_j), each alternative's value is linear in t; the best stays best between the
    // nearest crossings with the others on either side of the declared w_j.
    auto value = [&](std::size_t i, double t) {
      double v = t * r[i][j];
      for (std::size_t l = 0; l < k; l++)
        if (l != j) v += (w[j] < 1 ? w[l] * (1 - t) / (1 - w[j]) : (1 - t) / static_cast<double>(k - 1)) * r[i][l];
      return v;
    };
    double lo = 0, hi = 1;
    for (std::size_t i = 0; i < n; i++) {
      if (i == st.best) continue;
      const double d0 = value(st.best, 0) - value(i, 0), d1 = value(st.best, 1) - value(i, 1);
      if (d0 == d1) continue;  // parallel: never cross
      const double t = d0 / (d0 - d1);
      if (t > w[j] + 1e-15 && t < hi) hi = std::min(hi, t);
      if (t < w[j] - 1e-15 && t > lo) lo = std::max(lo, t);
    }
    st.lower.push_back(lo), st.upper.push_back(hi);
  }
  return st;
}

// ---- SMAA-2 rank acceptability ----

inline std::vector<std::vector<double>> acceptability(const Matrix& m, std::size_t draws = 10000, std::uint64_t seed = 1) {
  detail::validate(m);
  const auto r = detail::minmax(m);
  const std::size_t n = m.alternatives.size(), k = m.criteria.size();
  std::vector<std::vector<double>> b(n, std::vector<double>(n, 0.0));  // [alternative][rank − 1]
  calx::random::Stream rng(seed, 0, 7, 0);
  for (std::size_t t = 0; t < draws; t++) {
    // Uniform on the simplex: normalised exponentials (Rubinstein 1982).
    std::vector<double> w(k);
    double s = 0;
    for (double& v : w) v = -std::log(rng.uniform()), s += v;
    std::vector<double> val(n, 0.0);
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t j = 0; j < k; j++) val[i] += w[j] / s * r[i][j];
    const auto rk = detail::ranking(val);
    for (std::size_t i = 0; i < n; i++) b[i][rk.rank[i] - 1] += 1.0 / static_cast<double>(draws);
  }
  return b;
}

}  // namespace nomos::appraisal
