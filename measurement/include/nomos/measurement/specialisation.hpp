// Specialisation: how the structure of a part departs from the structure of the whole, and how
// far flows in two directions match.
//
// A table x of nonnegative amounts (a calx::linalg::Rectangular) has rows for units and columns
// for categories. Write x_r· and x_·c for the row and column sums and x_·· for the total.
//
// - The location quotient of unit r in category c (Isard 1960), which is also Balassa's
//   index of revealed comparative advantage (Balassa 1965), is
//     LQ_rc = (x_rc / x_r·) / (x_·c / x_··):
//   the share of category c in unit r over its share in the whole. LQ > 1 marks specialisation.
//   In each category its mean over the units, weighted by their sizes x_r· / x_··, is one; in
//   each unit its mean over the categories, weighted by x_·c / x_··, is one. Undefined (NaN)
//   where a row or column sums to zero. Its symmetric form (LQ − 1)/(LQ + 1), in [−1, 1),
//   treats specialisation and its absence alike (Laursen 2015).
// - The index of dissimilarity of two distributions a and b over the same categories, D =
//   ½ Σ_c |a_c / Σa − b_c / Σb| (Duncan and Duncan 1955), is the share of either that would have
//   to move between categories for the two to have the same structure; in [0, 1].
// - The Grubel–Lloyd index of two flows a_c and b_c in opposite directions in each category
//   (Grubel and Lloyd 1975) is GL_c = 1 − |a_c − b_c| / (a_c + b_c), in [0, 1]: 1 when the
//   flows match, 0 when one way only. The aggregate is the mean weighted by a_c + b_c,
//     GL = Σ (a_c + b_c − |a_c − b_c|) / Σ (a_c + b_c),
//   and, since an overall imbalance Σa ≠ Σb forces some categories to be one-way, Grubel and
//   Lloyd's adjusted aggregate divides by what the imbalance leaves possible:
//     C = Σ (a_c + b_c − |a_c − b_c|) / (Σ (a_c + b_c) − |Σa − Σb|)
//   (Grubel and Lloyd 1975). C = 1 when every category is as balanced as the totals allow.
//   (Aquino 1978 proposed a different correction, which is not here.)

#pragma once

#include <calx/linalg/rectangular.hpp>

#include <cmath>
#include <limits>
#include <span>
#include <vector>

#include "inequality.hpp"

namespace nomos::measurement {

// LQ_rc for a table of units × categories; NaN where a row or column sum is zero.
inline calx::linalg::Rectangular location_quotients(const calx::linalg::Rectangular& x) {
  detail::require(x.rows > 0 && x.cols > 0, "a nonempty table");
  std::vector<double> row(x.rows, 0.0), col(x.cols, 0.0);
  double total = 0;
  for (std::size_t r = 0; r < x.rows; r++)
    for (std::size_t c = 0; c < x.cols; c++) {
      const double v = x(r, c);
      detail::require(v >= 0 && std::isfinite(v), "nonnegative, finite amounts");
      row[r] += v;
      col[c] += v;
      total += v;
    }
  calx::linalg::Rectangular q(x.rows, x.cols);
  for (std::size_t r = 0; r < x.rows; r++)
    for (std::size_t c = 0; c < x.cols; c++)
      q(r, c) = row[r] > 0 && col[c] > 0 ? (x(r, c) / row[r]) / (col[c] / total) : std::numeric_limits<double>::quiet_NaN();
  return q;
}

// Balassa's index of revealed comparative advantage: the location quotient under its other name.
inline calx::linalg::Rectangular balassa_index(const calx::linalg::Rectangular& x) { return location_quotients(x); }

// (LQ − 1)/(LQ + 1), elementwise.
inline calx::linalg::Rectangular symmetric_location_quotients(calx::linalg::Rectangular lq) {
  for (double& v : lq.a) v = (v - 1) / (v + 1);
  return lq;
}

// D = ½ Σ |a/Σa − b/Σb|.
inline double dissimilarity(std::span<const double> a, std::span<const double> b) {
  detail::require(a.size() == b.size() && !a.empty(), "two distributions over the same categories");
  double sa = 0, sb = 0;
  for (std::size_t i = 0; i < a.size(); i++) {
    detail::require(a[i] >= 0 && b[i] >= 0, "nonnegative amounts");
    sa += a[i];
    sb += b[i];
  }
  detail::require(sa > 0 && sb > 0, "positive totals");
  double d = 0;
  for (std::size_t i = 0; i < a.size(); i++) d += std::fabs(a[i] / sa - b[i] / sb);
  return d / 2;
}

// GL = 1 − |a − b|/(a + b) for one category (NaN when a + b = 0).
inline double grubel_lloyd(double a, double b) {
  detail::require(a >= 0 && b >= 0, "nonnegative flows");
  return a + b > 0 ? 1 - std::fabs(a - b) / (a + b) : std::numeric_limits<double>::quiet_NaN();
}

namespace detail {
struct TwoWay {
  double gross = 0, net = 0, sa = 0, sb = 0;  // Σ(a + b), Σ|a − b|, Σa, Σb
};
inline TwoWay two_way(std::span<const double> a, std::span<const double> b) {
  require(a.size() == b.size() && !a.empty(), "flows in the same categories");
  TwoWay t;
  for (std::size_t i = 0; i < a.size(); i++) {
    require(a[i] >= 0 && b[i] >= 0, "nonnegative flows");
    t.gross += a[i] + b[i];
    t.net += std::fabs(a[i] - b[i]);
    t.sa += a[i];
    t.sb += b[i];
  }
  require(t.gross > 0, "some flow");
  return t;
}
}  // namespace detail

// The aggregate Grubel–Lloyd index Σ(a + b − |a − b|)/Σ(a + b).
inline double grubel_lloyd(std::span<const double> a, std::span<const double> b) {
  const auto t = detail::two_way(a, b);
  return (t.gross - t.net) / t.gross;
}

// Grubel and Lloyd's aggregate adjusted for the overall imbalance:
// Σ(a + b − |a − b|)/(Σ(a + b) − |Σa − Σb|); NaN when all flow is one way.
inline double grubel_lloyd_adjusted(std::span<const double> a, std::span<const double> b) {
  const auto t = detail::two_way(a, b);
  const double room = t.gross - std::fabs(t.sa - t.sb);
  return room > 0 ? (t.gross - t.net) / room : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace nomos::measurement
