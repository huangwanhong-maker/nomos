// Psychometrics: how well a set of items measures one trait, and where each person stands on it.
//
// Classical test theory. An observed score is a true score plus an error (Lord and Novick 1968, Statistical Theories of
// Mental Test Scores). For k items scored for n persons:
//   - Cronbach's α = k/(k − 1)·(1 − Σ var(item)/var(total)), a lower bound of reliability when the items are essentially
//     τ-equivalent and their errors uncorrelated (Cronbach 1951, Psychometrika 16, 297–334); the same divisor is used for
//     every variance, so it cancels;
//   - each item's correlation with the rest of the scale (the corrected item–total correlation) and α without it;
//   - Spearman and Brown's prophecy, the reliability of a test lengthened m times, m r/(1 + (m − 1) r), and the split-half
//     reliability of odd against even items lengthened back to the whole (Spearman 1910; Brown 1910);
//   - the standard error of measurement, sd·√(1 − reliability).
//
// Item response theory. A person's probability of answering an item correctly depends on the person's trait θ and the
// item's discrimination a and difficulty b: P(θ) = 1/(1 + exp(−a(θ − b))), the two-parameter logistic model (Birnbaum
// 1968); a = 1 for every item is Rasch's model (Rasch 1960). An item's information is a²P(1 − P), largest, a²/4, at θ = b;
// a test's is the sum, and 1/√I(θ) is the standard error of θ there.
//   - fit: marginal maximum likelihood with θ ~ N(0, 1), by the EM algorithm of Bock and Aitkin (1981, Psychometrika 46,
//     443–459) on a grid of Q points equally spaced on [−6, 6] weighted by the normal density (as the mirt package does by
//     default for one dimension: Chalmers 2012, Journal of Statistical Software 48(6)). Each M-step fits each item's
//     logit aθ + d by Newton's method on the expected counts at the nodes; b = −d/a. Under Rasch's model a common
//     discrimination is fitted instead, so that the latent scale stays N(0, 1). The marginal log-likelihood never
//     decreases from one iteration to the next.
//   - persons: the expected a posteriori trait and its posterior standard deviation on the same grid; a person with every
//     item right or every item wrong has a finite EAP and no maximum-likelihood estimate.
//
// Responses are 1 or 0, or −1 where an item was not presented; a missing response contributes nothing to a likelihood.

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::psychometrics {

using Responses = std::vector<std::vector<int>>;  // [person][item]: 1, 0, or −1 not presented

// ---- Classical test theory ---------------------------------------------------------------------------------------------------

namespace detail {
inline double variance(const std::vector<double>& x) {
  const double n = static_cast<double>(x.size());
  const double m = std::accumulate(x.begin(), x.end(), 0.0) / n;
  double s = 0;
  for (const double v : x) s += (v - m) * (v - m);
  return s / (n - 1);
}
inline double correlation(const std::vector<double>& x, const std::vector<double>& y) {
  const double n = static_cast<double>(x.size());
  const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n, my = std::accumulate(y.begin(), y.end(), 0.0) / n;
  double sxy = 0, sxx = 0, syy = 0;
  for (std::size_t i = 0; i < x.size(); i++) sxy += (x[i] - mx) * (y[i] - my), sxx += (x[i] - mx) * (x[i] - mx), syy += (y[i] - my) * (y[i] - my);
  return sxy / std::sqrt(sxx * syy);
}
inline void rectangular(const std::vector<std::vector<double>>& s, std::size_t least_items) {
  if (s.size() < 2) throw std::invalid_argument("nomos::psychometrics: at least two persons");
  const std::size_t k = s[0].size();
  if (k < least_items) throw std::invalid_argument("nomos::psychometrics: at least " + std::to_string(least_items) + " items");
  for (const auto& row : s)
    if (row.size() != k) throw std::invalid_argument("nomos::psychometrics: every person has a score on every item");
}
inline std::vector<double> column(const std::vector<std::vector<double>>& s, std::size_t j) {
  std::vector<double> c;
  for (const auto& row : s) c.push_back(row[j]);
  return c;
}
}  // namespace detail

// Scores [person][item], every item scored for every person.
inline double cronbach_alpha(const std::vector<std::vector<double>>& scores) {
  detail::rectangular(scores, 2);
  const std::size_t k = scores[0].size();
  std::vector<double> total(scores.size(), 0);
  double items = 0;
  for (std::size_t j = 0; j < k; j++) {
    const auto c = detail::column(scores, j);
    items += detail::variance(c);
    for (std::size_t i = 0; i < c.size(); i++) total[i] += c[i];
  }
  const double kk = static_cast<double>(k);
  return kk / (kk - 1) * (1 - items / detail::variance(total));
}

struct ItemAnalysis {
  std::vector<double> rest_correlation;  // each item with the sum of the others
  std::vector<double> alpha_without;     // α of the scale without the item (NaN with two items)
};

inline ItemAnalysis item_analysis(const std::vector<std::vector<double>>& scores) {
  detail::rectangular(scores, 2);
  const std::size_t k = scores[0].size();
  ItemAnalysis a;
  for (std::size_t j = 0; j < k; j++) {
    std::vector<double> rest(scores.size(), 0);
    std::vector<std::vector<double>> without;
    for (std::size_t i = 0; i < scores.size(); i++) {
      std::vector<double> row;
      for (std::size_t l = 0; l < k; l++)
        if (l != j) rest[i] += scores[i][l], row.push_back(scores[i][l]);
      without.push_back(std::move(row));
    }
    a.rest_correlation.push_back(detail::correlation(detail::column(scores, j), rest));
    a.alpha_without.push_back(k > 2 ? cronbach_alpha(without) : std::numeric_limits<double>::quiet_NaN());
  }
  return a;
}

inline double spearman_brown(double reliability, double lengthened) {
  if (!(lengthened > 0)) throw std::invalid_argument("nomos::psychometrics::spearman_brown: a positive factor of length");
  return lengthened * reliability / (1 + (lengthened - 1) * reliability);
}

// Odd items against even ones (the first item is odd), their correlation lengthened back to the whole test.
inline double split_half(const std::vector<std::vector<double>>& scores) {
  detail::rectangular(scores, 2);
  std::vector<double> odd(scores.size(), 0), even(scores.size(), 0);
  for (std::size_t i = 0; i < scores.size(); i++)
    for (std::size_t j = 0; j < scores[i].size(); j++) (j % 2 == 0 ? odd : even)[i] += scores[i][j];
  return spearman_brown(detail::correlation(odd, even), 2);
}

inline double standard_error_of_measurement(double sd, double reliability) {
  if (!(reliability >= 0 && reliability <= 1)) throw std::invalid_argument("nomos::psychometrics: a reliability in [0, 1]");
  return sd * std::sqrt(1 - reliability);
}

// ---- Item response theory ----------------------------------------------------------------------------------------------------------

struct Item {
  double a = 1, b = 0;  // discrimination, difficulty
};

inline double probability(const Item& it, double theta) { return 1 / (1 + std::exp(-it.a * (theta - it.b))); }

inline double information(const Item& it, double theta) {
  const double p = probability(it, theta);
  return it.a * it.a * p * (1 - p);
}

inline double information(const std::vector<Item>& items, double theta) {
  double s = 0;
  for (const auto& it : items) s += information(it, theta);
  return s;
}

struct Grid {
  std::vector<double> nodes, weights;  // weights sum to one: the N(0, 1) prior on the grid
};

inline Grid grid(std::size_t points = 61, double limit = 6) {
  if (points < 3 || !(limit > 0)) throw std::invalid_argument("nomos::psychometrics::grid: at least three points on a positive range");
  Grid g;
  double total = 0;
  for (std::size_t q = 0; q < points; q++) {
    const double t = -limit + 2 * limit * static_cast<double>(q) / static_cast<double>(points - 1);
    g.nodes.push_back(t);
    g.weights.push_back(std::exp(-t * t / 2));
    total += g.weights.back();
  }
  for (auto& w : g.weights) w /= total;
  return g;
}

enum class Model { rasch, two_parameter };

struct Fit {
  std::vector<Item> items;
  double log_likelihood = 0;
  std::vector<double> trace;  // the marginal log-likelihood at each iteration
  std::size_t iterations = 0;
  bool converged = false;
};

namespace detail {
inline void validate(const Responses& r) {
  if (r.empty() || r[0].empty()) throw std::invalid_argument("nomos::psychometrics: no responses");
  for (const auto& row : r) {
    if (row.size() != r[0].size()) throw std::invalid_argument("nomos::psychometrics: every person has a response, or −1, for every item");
    for (const int x : row)
      if (x != 0 && x != 1 && x != -1) throw std::invalid_argument("nomos::psychometrics: responses are 1, 0 or −1");
  }
}
// log P(responses of person | θ_q) for every node.
inline std::vector<double> log_likelihoods(const std::vector<int>& row, const std::vector<Item>& items, const Grid& g) {
  std::vector<double> l(g.nodes.size(), 0);
  for (std::size_t q = 0; q < g.nodes.size(); q++)
    for (std::size_t j = 0; j < items.size(); j++) {
      if (row[j] < 0) continue;
      const double z = items[j].a * (g.nodes[q] - items[j].b);
      // log σ(z) and log(1 − σ(z)) without overflow.
      l[q] += row[j] ? -std::log1p(std::exp(-z)) : -std::log1p(std::exp(z));
    }
  return l;
}
// The posterior over the nodes, and the person's contribution log Σ_q w_q L_q to the marginal log-likelihood.
inline double posterior(const std::vector<int>& row, const std::vector<Item>& items, const Grid& g, std::vector<double>& post) {
  const auto l = log_likelihoods(row, items, g);
  const double top = *std::max_element(l.begin(), l.end());
  double z = 0;
  post.assign(l.size(), 0);
  for (std::size_t q = 0; q < l.size(); q++) post[q] = g.weights[q] * std::exp(l[q] - top), z += post[q];
  for (auto& p : post) p /= z;
  return top + std::log(z);
}
}  // namespace detail

inline double marginal_log_likelihood(const Responses& r, const std::vector<Item>& items, const Grid& g = grid()) {
  detail::validate(r);
  std::vector<double> post;
  double ll = 0;
  for (const auto& row : r) ll += detail::posterior(row, items, g, post);
  return ll;
}

inline Fit fit(const Responses& r, Model model = Model::two_parameter, const Grid& g = grid(), std::size_t most = 500, double tolerance = 1e-8) {
  detail::validate(r);
  const std::size_t n = r.size(), k = r[0].size(), Q = g.nodes.size();
  Fit f;
  f.items.assign(k, Item{1, 0});
  for (std::size_t j = 0; j < k; j++) {  // start from the proportion correct
    double right = 0, seen = 0;
    for (const auto& row : r)
      if (row[j] >= 0) right += row[j], seen++;
    const double p = std::clamp((right + 0.5) / (seen + 1), 0.01, 0.99);
    f.items[j].b = -std::log(p / (1 - p));
  }
  std::vector<double> post;
  double previous = -std::numeric_limits<double>::infinity();
  for (f.iterations = 0; f.iterations < most; f.iterations++) {
    // E-step: expected numbers presented (n_jq) and right (r_jq) at each node.
    std::vector<std::vector<double>> seen(k, std::vector<double>(Q, 0)), right(k, std::vector<double>(Q, 0));
    double ll = 0;
    for (std::size_t i = 0; i < n; i++) {
      ll += detail::posterior(r[i], f.items, g, post);
      for (std::size_t j = 0; j < k; j++) {
        if (r[i][j] < 0) continue;
        for (std::size_t q = 0; q < Q; q++) {
          seen[j][q] += post[q];
          if (r[i][j]) right[j][q] += post[q];
        }
      }
    }
    f.trace.push_back(ll);
    if (ll - previous < tolerance * (1 + std::abs(ll))) {
      f.converged = true;
      f.log_likelihood = ll;
      break;
    }
    previous = ll;
    // M-step: logit aθ + d for each item by Newton's method on the expected counts; Rasch: a common a, each d.
    if (model == Model::two_parameter) {
      for (std::size_t j = 0; j < k; j++) {
        double a = f.items[j].a, d = -f.items[j].a * f.items[j].b;
        for (int step = 0; step < 25; step++) {
          double ga = 0, gd = 0, haa = 0, had = 0, hdd = 0;
          for (std::size_t q = 0; q < Q; q++) {
            const double t = g.nodes[q], p = 1 / (1 + std::exp(-(a * t + d))), e = right[j][q] - seen[j][q] * p, v = seen[j][q] * p * (1 - p);
            ga += e * t, gd += e, haa += v * t * t, had += v * t, hdd += v;
          }
          const double det = haa * hdd - had * had;
          if (!(det > 0)) break;
          const double da = (hdd * ga - had * gd) / det, dd = (haa * gd - had * ga) / det;
          a += da, d += dd;
          if (std::abs(da) + std::abs(dd) < 1e-10) break;
        }
        f.items[j] = {a, -d / a};
      }
    } else {
      double a = f.items[0].a;
      std::vector<double> d(k);
      for (std::size_t j = 0; j < k; j++) d[j] = -f.items[j].a * f.items[j].b;
      for (int step = 0; step < 50; step++) {
        // Newton on (a, d_1, …, d_k): the Hessian is an arrow, solved by eliminating the d's.
        double ga = 0, haa = 0, change = 0;
        std::vector<double> gd(k, 0), hdd(k, 0), had(k, 0);
        for (std::size_t j = 0; j < k; j++)
          for (std::size_t q = 0; q < Q; q++) {
            const double t = g.nodes[q], p = 1 / (1 + std::exp(-(a * t + d[j]))), e = right[j][q] - seen[j][q] * p, v = seen[j][q] * p * (1 - p);
            ga += e * t, gd[j] += e, haa += v * t * t, had[j] += v * t, hdd[j] += v;
          }
        double schur = haa, rhs = ga;
        for (std::size_t j = 0; j < k; j++) schur -= had[j] * had[j] / hdd[j], rhs -= had[j] * gd[j] / hdd[j];
        const double da = rhs / schur;
        a += da;
        change += std::abs(da);
        for (std::size_t j = 0; j < k; j++) {
          const double dd = (gd[j] - had[j] * da) / hdd[j];
          d[j] += dd;
          change += std::abs(dd);
        }
        if (change < 1e-10) break;
      }
      for (std::size_t j = 0; j < k; j++) f.items[j] = {a, -d[j] / a};
    }
  }
  if (!f.converged) f.log_likelihood = marginal_log_likelihood(r, f.items, g);
  return f;
}

struct Score {
  double eap = 0, sd = 0;  // the posterior mean of θ and its standard deviation
};

inline std::vector<Score> score(const Responses& r, const std::vector<Item>& items, const Grid& g = grid()) {
  detail::validate(r);
  if (r[0].size() != items.size()) throw std::invalid_argument("nomos::psychometrics::score: one item for each response");
  std::vector<Score> out;
  std::vector<double> post;
  for (const auto& row : r) {
    detail::posterior(row, items, g, post);
    double m = 0, s2 = 0;
    for (std::size_t q = 0; q < post.size(); q++) m += post[q] * g.nodes[q];
    for (std::size_t q = 0; q < post.size(); q++) s2 += post[q] * (g.nodes[q] - m) * (g.nodes[q] - m);
    out.push_back({m, std::sqrt(s2)});
  }
  return out;
}

}  // namespace nomos::psychometrics
