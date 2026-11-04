// Bilateral index numbers: how much prices (or quantities) changed between two situations.
//
// Two situations 0 and 1 (periods, places) have prices p⁰, p¹ and quantities q⁰, q¹ of the same N
// items. A price index is a function P(p⁰, p¹, q⁰, q¹) > 0, written in this order throughout, as
// in the axiomatic literature (Diewert 1992; ILO et al. 2004, ch. 16). With p·q = Σ_i p_i q_i,
// expenditure shares s_i^t = p_i^t q_i^t / p^t·q^t and the value ratio V = p¹·q¹ / p⁰·q⁰:
//
//   Laspeyres            P_L = p¹·q⁰ / p⁰·q⁰                        (the basket of situation 0)
//   Paasche              P_P = p¹·q¹ / p⁰·q¹                        (the basket of situation 1)
//   Fisher               P_F = (P_L P_P)^{1/2}                      (Fisher 1922)
//   Törnqvist            ln P_T = Σ_i ½(s_i⁰ + s_i¹) ln(p_i¹ / p_i⁰)   (Törnqvist 1936)
//   Walsh                P_W = Σ p_i¹ (q_i⁰ q_i¹)^{1/2} / Σ p_i⁰ (q_i⁰ q_i¹)^{1/2}   (Walsh 1901)
//   Marshall–Edgeworth   P_ME = p¹·(q⁰ + q¹) / p⁰·(q⁰ + q¹)
//
// (ILO et al. 2004, ch. 15, which derives each and relates them.) The corresponding quantity index Q is the same formula with the roles of prices and quantities
// exchanged, Q(p⁰, p¹, q⁰, q¹) = P(q⁰, q¹, p⁰, p¹). Fisher's index alone among these satisfies
// the factor reversal test P Q = V (Fisher 1922); Fisher's, Törnqvist's and Walsh's are
// superlative (exact for flexible functional forms; Diewert 1976) and satisfy time reversal,
// P(p⁰, p¹, q⁰, q¹) P(p¹, p⁰, q¹, q⁰) = 1; Laspeyres and Paasche do not, their product with the
// reversed index being P_L / P_P.
//
// Indices with a basket from elsewhere:
//   Lowe     P_Lo = p¹·q^b / p⁰·q^b, the basket q^b of a reference period b (Lowe 1823; ILO et
//            al. 2004, ch. 15), which is how most consumer price indices are compiled;
//   Young    P_Y = Σ_i s_i^b (p_i¹ / p_i⁰), with the expenditure shares s^b of a period b
//            (Young 1812; ILO et al. 2004, ch. 15).
//
// Elementary indices, of prices alone (ILO et al. 2004, ch. 20), for n quotes of one item:
//   Jevons   (Π p_i¹ / p_i⁰)^{1/n}           Dutot   Σ p_i¹ / Σ p_i⁰
//   Carli    (1/n) Σ p_i¹ / p_i⁰             harmonic  (1/n Σ p_i⁰ / p_i¹)⁻¹
//   CSWD     (Carli · harmonic)^{1/2}  (Carruthers, Sellwood and Ward 1980; Dalén 1992).
// Jevons is transitive and satisfies time reversal; Carli fails time reversal with an upward
// bias, Carli(0, 1) Carli(1, 0) ≥ 1 (by the inequality of arithmetic and harmonic means); Dutot
// is not invariant to the units of the items (fails commensurability).
//
// Accuracy: plain sums, relative error of order N u. Every price must be positive; quantities
// nonnegative, with positive values p·q where they are divided by.

#pragma once

#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace nomos::indices {

namespace detail {

inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::indices: ") + what);
}

inline double dot(std::span<const double> a, std::span<const double> b) {
  double s = 0;
  for (std::size_t i = 0; i < a.size(); i++) s += a[i] * b[i];
  return s;
}

inline void items(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  const std::size_t n = p0.size();
  require(n > 0 && p1.size() == n && q0.size() == n && q1.size() == n, "prices and quantities of the same items");
  for (std::size_t i = 0; i < n; i++) {
    require(p0[i] > 0 && p1[i] > 0 && std::isfinite(p0[i]) && std::isfinite(p1[i]), "positive, finite prices");
    require(q0[i] >= 0 && q1[i] >= 0 && std::isfinite(q0[i]) && std::isfinite(q1[i]), "nonnegative, finite quantities");
  }
}

}  // namespace detail

enum class Formula { laspeyres, paasche, fisher, tornqvist, walsh, marshall_edgeworth };

// P_L = p¹·q⁰ / p⁰·q⁰.
inline double laspeyres(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  return detail::dot(p1, q0) / detail::dot(p0, q0);
}

// P_P = p¹·q¹ / p⁰·q¹.
inline double paasche(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  return detail::dot(p1, q1) / detail::dot(p0, q1);
}

// P_F = (P_L P_P)^{1/2}.
inline double fisher(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  return std::sqrt(laspeyres(p0, p1, q0, q1) * paasche(p0, p1, q0, q1));
}

// ln P_T = Σ ½(s⁰ + s¹) ln(p¹/p⁰).
inline double tornqvist(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  const double v0 = detail::dot(p0, q0), v1 = detail::dot(p1, q1);
  double l = 0;
  for (std::size_t i = 0; i < p0.size(); i++) l += 0.5 * (p0[i] * q0[i] / v0 + p1[i] * q1[i] / v1) * std::log(p1[i] / p0[i]);
  return std::exp(l);
}

// P_W = Σ p¹ √(q⁰q¹) / Σ p⁰ √(q⁰q¹).
inline double walsh(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  double a = 0, b = 0;
  for (std::size_t i = 0; i < p0.size(); i++) {
    const double m = std::sqrt(q0[i] * q1[i]);
    a += p1[i] * m;
    b += p0[i] * m;
  }
  return a / b;
}

// P_ME = p¹·(q⁰ + q¹) / p⁰·(q⁰ + q¹).
inline double marshall_edgeworth(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  double a = 0, b = 0;
  for (std::size_t i = 0; i < p0.size(); i++) {
    a += p1[i] * (q0[i] + q1[i]);
    b += p0[i] * (q0[i] + q1[i]);
  }
  return a / b;
}

// The price index of a formula.
inline double price_index(Formula f, std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  switch (f) {
    case Formula::laspeyres: return laspeyres(p0, p1, q0, q1);
    case Formula::paasche: return paasche(p0, p1, q0, q1);
    case Formula::fisher: return fisher(p0, p1, q0, q1);
    case Formula::tornqvist: return tornqvist(p0, p1, q0, q1);
    case Formula::walsh: return walsh(p0, p1, q0, q1);
    case Formula::marshall_edgeworth: return marshall_edgeworth(p0, p1, q0, q1);
  }
  throw std::invalid_argument("nomos::indices::price_index: an unknown formula");
}

// The quantity index of a formula: P with prices and quantities exchanged. Quantities must then
// be positive where they enter as "prices".
inline double quantity_index(Formula f, std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  return price_index(f, q0, q1, p0, p1);
}

// V = p¹·q¹ / p⁰·q⁰.
inline double value_ratio(std::span<const double> p0, std::span<const double> p1, std::span<const double> q0, std::span<const double> q1) {
  detail::items(p0, p1, q0, q1);
  return detail::dot(p1, q1) / detail::dot(p0, q0);
}

// The Lowe index p¹·q^b / p⁰·q^b, of a basket q^b from a reference period.
inline double lowe(std::span<const double> p0, std::span<const double> p1, std::span<const double> qb) {
  detail::items(p0, p1, qb, qb);
  return detail::dot(p1, qb) / detail::dot(p0, qb);
}

// The Young index Σ s^b (p¹/p⁰) / Σ s^b, of expenditure shares s^b from a reference period.
inline double young(std::span<const double> p0, std::span<const double> p1, std::span<const double> sb) {
  detail::items(p0, p1, sb, sb);
  double a = 0, b = 0;
  for (std::size_t i = 0; i < p0.size(); i++) {
    a += sb[i] * p1[i] / p0[i];
    b += sb[i];
  }
  return a / b;
}

// Elementary indices of n price quotes.
namespace detail {
inline void quotes(std::span<const double> p0, std::span<const double> p1) {
  require(!p0.empty() && p0.size() == p1.size(), "quotes of the same items");
  for (std::size_t i = 0; i < p0.size(); i++) require(p0[i] > 0 && p1[i] > 0, "positive prices");
}
}  // namespace detail

inline double jevons(std::span<const double> p0, std::span<const double> p1) {
  detail::quotes(p0, p1);
  double l = 0;
  for (std::size_t i = 0; i < p0.size(); i++) l += std::log(p1[i] / p0[i]);
  return std::exp(l / static_cast<double>(p0.size()));
}

inline double dutot(std::span<const double> p0, std::span<const double> p1) {
  detail::quotes(p0, p1);
  double a = 0, b = 0;
  for (std::size_t i = 0; i < p0.size(); i++) {
    a += p1[i];
    b += p0[i];
  }
  return a / b;
}

inline double carli(std::span<const double> p0, std::span<const double> p1) {
  detail::quotes(p0, p1);
  double a = 0;
  for (std::size_t i = 0; i < p0.size(); i++) a += p1[i] / p0[i];
  return a / static_cast<double>(p0.size());
}

inline double harmonic(std::span<const double> p0, std::span<const double> p1) {
  detail::quotes(p0, p1);
  double a = 0;
  for (std::size_t i = 0; i < p0.size(); i++) a += p0[i] / p1[i];
  return static_cast<double>(p0.size()) / a;
}

// Carruthers, Sellwood, Ward and Dalén: the geometric mean of Carli and harmonic.
inline double cswd(std::span<const double> p0, std::span<const double> p1) { return std::sqrt(carli(p0, p1) * harmonic(p0, p1)); }

}  // namespace nomos::indices
