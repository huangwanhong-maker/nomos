// Decompositions by groups: inequality into a part within groups and a part between them, and
// poverty into the contributions of groups.
//
// Partition the units into groups g, with population shares v_g = W_g / W, means μ_g and shares
// of the total s_g = v_g μ_g / μ. Every index of the generalised entropy family decomposes
// additively (Bourguignon 1979 for Theil's indices; Shorrocks 1980 for the whole family, which he
// shows is the only additively decomposable class of relative indices):
//
//   GE(α) = Σ_g v_g^{1−α} s_g^α GE_g(α) + GE_B(α),
//
// where GE_g is the index within group g and GE_B the index of the distribution in which every
// unit holds its group's mean (Shorrocks 1980; Cowell 2011, ch. 3). For Theil's T (α = 1) the
// weights are the shares of the total s_g; for Theil's L (α = 0), the population shares v_g, so
// that only for α = 0 do the weights sum to one whatever the groups. The Gini coefficient does
// not decompose so (its residual depends on how the groups overlap), and is not here.
//
// The Foster–Greer–Thorbecke indices decompose as FGT(α) = Σ_g v_g FGT_g(α) (Foster, Greer and
// Thorbecke 1984), the contribution of group g being v_g FGT_g(α).

#pragma once

#include <cmath>
#include <span>
#include <vector>

#include "inequality.hpp"
#include "poverty.hpp"

namespace nomos::measurement {

// The decomposition of GE(α) over groups 0 … G − 1 (G = 1 + the largest label).
struct EntropyDecomposition {
  double alpha = 0;
  double total = 0, within = 0, between = 0;  // total = within + between
  std::vector<double> population_share;       // v_g
  std::vector<double> amount_share;           // s_g
  std::vector<double> mean;                   // μ_g
  std::vector<double> index;                  // GE_g(α); 0 for an empty group
  std::vector<double> weight;                 // v_g^{1−α} s_g^α
};

namespace detail {
struct Groups {
  std::size_t count = 0;
  std::vector<std::vector<double>> x, w;  // the amounts and weights of each group
  std::vector<double> weight, sum;         // W_g, Σ w x in g
};
inline Groups split(std::span<const double> x, std::span<const std::size_t> group, const std::vector<double>& w) {
  require(group.size() == x.size(), "a group for every unit");
  Groups g;
  for (std::size_t k : group) g.count = std::max(g.count, k + 1);
  g.x.resize(g.count);
  g.w.resize(g.count);
  g.weight.assign(g.count, 0.0);
  g.sum.assign(g.count, 0.0);
  for (std::size_t i = 0; i < x.size(); i++) {
    g.x[group[i]].push_back(x[i]);
    g.w[group[i]].push_back(w[i]);
    g.weight[group[i]] += w[i];
    g.sum[group[i]] += w[i] * x[i];
  }
  return g;
}
}  // namespace detail

inline EntropyDecomposition decompose_generalised_entropy(std::span<const double> x, std::span<const std::size_t> group, double alpha,
                                                          std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  detail::require(t.mean > 0, "a positive mean");
  const auto g = detail::split(x, group, ww);
  EntropyDecomposition d;
  d.alpha = alpha;
  d.total = generalised_entropy(x, alpha, ww);
  d.population_share.assign(g.count, 0.0);
  d.amount_share.assign(g.count, 0.0);
  d.mean.assign(g.count, 0.0);
  d.index.assign(g.count, 0.0);
  d.weight.assign(g.count, 0.0);
  std::vector<double> means, masses;  // the smoothed distribution: each nonempty group at its mean
  for (std::size_t k = 0; k < g.count; k++) {
    if (g.weight[k] == 0) continue;
    d.population_share[k] = g.weight[k] / t.weight;
    d.mean[k] = g.sum[k] / g.weight[k];
    d.amount_share[k] = d.population_share[k] * d.mean[k] / t.mean;
    detail::require(d.mean[k] > 0, "a positive mean in every group");
    d.index[k] = generalised_entropy(g.x[k], alpha, g.w[k]);
    const double v = d.population_share[k], s = d.amount_share[k];
    d.weight[k] = std::pow(v, 1 - alpha) * std::pow(s, alpha);
    d.within += d.weight[k] * d.index[k];
    means.push_back(d.mean[k]);
    masses.push_back(g.weight[k]);
  }
  d.between = generalised_entropy(means, alpha, masses);
  return d;
}

// The decomposition of FGT(α) at the line z over groups: FGT = Σ_g v_g FGT_g.
struct PovertyDecomposition {
  double total = 0;
  std::vector<double> population_share;  // v_g
  std::vector<double> index;             // FGT_g(α); 0 for an empty group
  std::vector<double> contribution;      // v_g FGT_g(α), summing to the total
};

inline PovertyDecomposition decompose_fgt(std::span<const double> x, std::span<const std::size_t> group, double z, double alpha,
                                          std::span<const double> w = {}) {
  const auto ww = detail::weights(x, w);
  const auto t = detail::totals(x, ww);
  const auto g = detail::split(x, group, ww);
  PovertyDecomposition d;
  d.total = fgt(x, z, alpha, ww);
  d.population_share.assign(g.count, 0.0);
  d.index.assign(g.count, 0.0);
  d.contribution.assign(g.count, 0.0);
  for (std::size_t k = 0; k < g.count; k++) {
    if (g.weight[k] == 0) continue;
    d.population_share[k] = g.weight[k] / t.weight;
    d.index[k] = fgt(g.x[k], z, alpha, g.w[k]);
    d.contribution[k] = d.population_share[k] * d.index[k];
  }
  return d;
}

}  // namespace nomos::measurement
