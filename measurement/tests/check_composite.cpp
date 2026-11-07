// Checks of composite indicators (composite.hpp).
//
// - A factor scaled between its percentiles: linearly between them, clipped beyond, an outlier leaving the rest in place; inverted, its complement.
// - Percentiles interpolated between order statistics, as the Handbook's robust min-max.
// - The composite is the weighted mean of the scaled factors; a factor missing for a unit is left out of its mean, its weights renormalised; a unit none
//   of whose factors is known stays missing.
// - A weight of nought leaves a factor out; doubling every weight changes nothing; a linear aggregation compensates: much of one offsets little of another.
// - Refusals: mismatched lengths, negative weights, percentiles out of order.

#include <nomos/measurement/composite.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nm = nomos::measurement;

int main() {
  calx::verification::Suite s;
  const double no = std::numeric_limits<double>::quiet_NaN();
  {
    const nm::Factor f{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 1, false};
    const auto u = nm::scaled(f, 0.1, 0.9, 1.0);
    s.check("the lower percentile goes to nought", u[1], 1e-12);
    s.check("the upper to the top", u[9] - 1, 1e-12);
    s.check("linearly between them", u[5] - 0.5, 1e-12);
    s.check("clipped below", u[0], 0);
    s.check("and above", u[10] - 1, 0);
    nm::Factor g = f;
    g.invert = true;
    const auto v = nm::scaled(g, 0.1, 0.9, 1.0);
    s.check("inverted, its complement", v[3] + u[3] - 1, 1e-12);
    nm::Factor h = f;
    h.values[10] = 1e9;
    s.check("an outlier leaves the rest where they were", nm::scaled(h, 0.1, 0.9, 1.0)[5] - u[5], 1e-12);
    s.check("the top of the scale is asked for", nm::scaled(f, 0.1, 0.9, 5.0)[5] - 2.5, 1e-12);
    const nm::Factor q{{0, 10}, 1, false};
    s.check("percentiles between order statistics", nm::scaled(q, 0.25, 0.75, 1.0)[0], 0);
    const nm::Factor same{{3, 3, 3}, 1, false};
    s.check("a factor with no spread sits in the middle", nm::scaled(same)[1] - 0.5, 0);
  }
  {
    const nm::Factor a{{0, 1, 2, no}, 1, false}, b{{4, 2, 0, no}, 3, false};
    const auto c = nm::composite({a, b}, 0, 1, 5);
    s.check("the weighted mean of the scaled factors", c[0] - (0 * 1 + 5 * 3) / 4.0, 1e-12);
    s.check("a linear aggregation compensates", c[1] - 2.5, 1e-12);
    s.require("a unit none of whose factors is known stays missing", std::isnan(c[3]));
    const nm::Factor d{{0, 2, no}, 1, false}, e{{0, 2, 2}, 1, false};
    const auto m = nm::composite({d, e}, 0, 1, 5);
    s.check("a factor missing for a unit is left out of its mean", m[2] - 5, 1e-12);
    const nm::Factor z{{9, 1, 5}, 0, false}, y{{0, 1, 2}, 1, false};
    s.check("a weight of nought leaves a factor out", nm::composite({z, y}, 0, 1, 5)[0], 0);
    nm::Factor a2 = a, b2 = b;
    a2.weight = 2, b2.weight = 6;
    s.check("doubling every weight changes nothing", nm::composite({a2, b2}, 0, 1, 5)[1] - c[1], 1e-12);
  }
  {
    bool refused = false;
    try { nm::composite({nm::Factor{{1, 2}, 1, false}, nm::Factor{{1}, 1, false}}); } catch (const std::invalid_argument&) { refused = true; }
    s.require("factors of different lengths are refused", refused);
    refused = false;
    try { nm::composite({nm::Factor{{1, 2}, -1, false}}); } catch (const std::invalid_argument&) { refused = true; }
    s.require("a negative weight is refused", refused);
    refused = false;
    try { nm::scaled(nm::Factor{{1, 2}, 1, false}, 0.9, 0.1); } catch (const std::invalid_argument&) { refused = true; }
    s.require("percentiles out of order are refused", refused);
  }
  return s.done();
}
