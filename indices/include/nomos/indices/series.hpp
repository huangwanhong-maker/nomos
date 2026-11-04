// Series of index numbers over periods: direct and chained.
//
// Prices p[t][i] and quantities q[t][i] of the same N items in periods t = 0, …, T − 1 give two
// series for a bilateral formula P (bilateral.hpp):
//
//   direct    I_t = P(p^b, p^t, q^b, q^t), every period against the base b;
//   chained   I_t = Π_{s=1}^{t} P(p^{s−1}, p^s, q^{s−1}, q^s), the product of the links between
//             consecutive periods, divided by its value at the base (ILO et al. 2004, ch. 15, on
//             the chain principle; Szulc 1983).
//
// The two agree for every data set exactly when the formula is transitive, P(0, 2) = P(0, 1)
// P(1, 2) (Fisher's circularity test), which no formula using quantities satisfies in general
// (Eichhorn and Voeller 1976); Jevons' elementary index does. When prices and quantities return
// to where they started, a direct index returns to one (the identity test) but a chained one
// need not: the chain drift of Szulc (1983), largest for Laspeyres and Paasche and small for the
// superlative indices when changes are smooth (ILO et al. 2004, ch. 19).

#pragma once

#include <vector>

#include "bilateral.hpp"

namespace nomos::indices {

// p[t][i] or q[t][i]: one vector of the N items for each period.
using Periods = std::vector<std::vector<double>>;

// The links P(p^{t−1}, p^t, q^{t−1}, q^t), t = 1, …, T − 1.
inline std::vector<double> chain_links(Formula f, const Periods& p, const Periods& q) {
  detail::require(!p.empty() && p.size() == q.size(), "prices and quantities for the same periods");
  std::vector<double> links;
  for (std::size_t t = 1; t < p.size(); t++) links.push_back(price_index(f, p[t - 1], p[t], q[t - 1], q[t]));
  return links;
}

// I_t = P(p^b, p^t, q^b, q^t).
inline std::vector<double> direct_series(Formula f, const Periods& p, const Periods& q, std::size_t base = 0) {
  detail::require(!p.empty() && p.size() == q.size() && base < p.size(), "prices and quantities for the same periods, and a base among them");
  std::vector<double> I(p.size());
  for (std::size_t t = 0; t < p.size(); t++) I[t] = t == base ? 1.0 : price_index(f, p[base], p[t], q[base], q[t]);
  return I;
}

// I_t = Π links, divided by its value at the base.
inline std::vector<double> chained_series(Formula f, const Periods& p, const Periods& q, std::size_t base = 0) {
  detail::require(base < p.size(), "a base among the periods");
  const auto links = chain_links(f, p, q);
  std::vector<double> I(p.size(), 1.0);
  for (std::size_t t = 1; t < p.size(); t++) I[t] = I[t - 1] * links[t - 1];
  const double at_base = I[base];
  for (double& v : I) v /= at_base;
  return I;
}

}  // namespace nomos::indices
