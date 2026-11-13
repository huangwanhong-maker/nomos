// Checks of nomos::measurement's agreement among coders (agreement.hpp): Cohen's κ, unweighted and with linear and quadratic weights,
// Fleiss's κ, and Krippendorff's α on his own example (four coders, twelve units, missing values) at the four levels of measurement,
// against statsmodels 0.14 and the krippendorff package 0.9 (appraisal/tests/fixtures/make_appraisal_references.py, 2026-11-13);
// and perfect agreement, κ = α = 1.

#include <nomos/measurement/agreement.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <vector>

using namespace nomos::measurement;

int main() {
  calx::verification::Suite s;
  const std::vector<int> a{0, 1, 2, 1, 0, 2, 2, 1, 0, 1, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1};
  const std::vector<int> b{0, 1, 2, 2, 0, 1, 2, 1, 1, 1, 2, 0, 0, 1, 2, 0, 1, 2, 1, 2};
  s.check("Cohen's κ against statsmodels: 0.5454545…", std::fabs(cohen(a, b, 3) - 0.5454545454545454), 1e-14);
  s.check("Cohen's linearly weighted κ against statsmodels: 0.6428571…", std::fabs(cohen(a, b, 3, Disagreement::linear) - 0.6428571428571429), 1e-14);
  s.check("Cohen's quadratically weighted κ against statsmodels: 0.75", std::fabs(cohen(a, b, 3, Disagreement::quadratic) - 0.75), 1e-14);
  s.check("Cohen's κ of a coder with herself is 1", std::fabs(cohen(a, a, 3) - 1), 1e-15);
  const std::vector<std::vector<double>> counts{{5, 0, 0, 0}, {2, 2, 1, 0}, {0, 4, 1, 0}, {1, 1, 1, 2}, {0, 0, 5, 0},
                                                {3, 2, 0, 0}, {0, 0, 0, 5}, {1, 3, 1, 0}, {0, 1, 4, 0}, {2, 0, 0, 3}};
  s.check("Fleiss's κ against statsmodels: 0.4105037…", std::fabs(fleiss(counts) - 0.4105037513397641), 1e-14);
  const int n = -1;
  const std::vector<std::vector<int>> data{{1, 2, 3, 3, 2, 1, 4, 1, 2, n, n, n}, {1, 2, 3, 3, 2, 2, 4, 1, 2, 5, n, 3},
                                           {n, 3, 3, 3, 2, 3, 4, 2, 2, 5, 1, n}, {1, 2, 3, 3, 2, 4, 4, 1, 2, 5, 1, n}};
  s.check("Krippendorff's α, nominal, on his example: 0.7434 (the krippendorff package)", std::fabs(krippendorff(data, 6, Metric::nominal) - 0.743421052631579), 1e-12);
  s.check("Krippendorff's α, ordinal: 0.8154", std::fabs(krippendorff(data, 6, Metric::ordinal) - 0.8153875037548814), 1e-12);
  s.check("Krippendorff's α, interval: 0.8491", std::fabs(krippendorff(data, 6, Metric::interval) - 0.8491071428571428), 1e-12);
  s.check("Krippendorff's α, ratio: 0.7974", std::fabs(krippendorff(data, 6, Metric::ratio) - 0.7974027747116121), 1e-12);
  const std::vector<std::vector<int>> same{{0, 1, 2, 1}, {0, 1, 2, 1}, {0, 1, 2, n}};
  s.check("Krippendorff's α of coders who always agree is 1", std::fabs(krippendorff(same, 3) - 1), 1e-15);
  return s.done();
}
