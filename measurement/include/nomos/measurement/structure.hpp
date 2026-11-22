// The structure of regional economies: how a region's employment changed against its reference (shift-share), how varied it is and
// how related its parts are (variety), which activities lie near one another (relatedness), and how complex its set of activities is
// (economic complexity) (decision 0018).
//
// - Shift-share (Dunn 1960; Esteban-Marquillas 1972). A region's employment e_i by industry at two dates, and its reference's (the
//   nation) E_i, with growth rates G (the reference's total), G_i (the reference's industry i). The change Δe_i is
//     national share e_i G + industrial mix e_i (G_i − G) + competitive Δe_i − e_i G_i,
//   an identity (the competitive term written so that an industry new to the region has it all). Esteban-Marquillas split the
//   competitive term with the homothetic employment e*_i = e E_i / E: competitive effect e*_i (g_i − G_i) and allocation effect
//   (e_i − e*_i)(g_i − G_i), g_i the region's growth in industry i, defined where every e_i is positive.
// - Variety (Frenken, Van Oort and Verburg 2007). Shares p_i over fine classes nested in groups g with shares P_g: the entropy
//   H = Σ p_i log(1/p_i) is unrelated variety UV = Σ P_g log(1/P_g) plus related variety RV = Σ P_g H_g, H_g the entropy within g
//   (the chain rule of entropy; exact). In natural logarithms unless a base is given (Frenken et al. use base 2).
// - Relatedness (Hidalgo, Klinger, Barabási and Hausmann 2007). The binary incidence M_rp = [RCA_rp ≥ threshold] (RCA the location
//   quotient, nomos's `specialisation`), the diversity k_r = Σ_p M_rp and the ubiquity k_p = Σ_r M_rp. The proximity
//   φ_pq = Σ_r M_rp M_rq / max(k_p, k_q), the lesser of the two conditional probabilities; the density of region r around p,
//   ω_rp = Σ_{q ≠ p} φ_pq M_rq / Σ_{q ≠ p} φ_pq (p itself left out: a convention, stated).
// - Economic complexity (Hidalgo and Hausmann 2009). The index of regions is the second eigenvector of the walk D⁻¹ M U⁻¹ Mᵀ, which
//   Mealy, Farmer and Teytelboym (2019) show to be the relaxation of a normalised cut and nothing more: tyche's `bipartite_ordering`.
//   Its raw form has zero dot product with diversity (the walk's stationary law is proportional to diversity); the index of
//   activities is proportional to the mean index of the regions that hold them, and the index of regions to the mean index of their
//   activities, exactly. Reported standardised (mean 0, variance 1, divisor n), the regions' signed to rise with diversity (tyche's convention) and the
//   activities' taken as the mean regional index of the regions that hold them, so that their signs agree. An incidence
//   that falls apart into blocks is refused, since its second eigenvector would only name a block.
// - Fitness and complexity (Tacchella, Cristelli, Caldarelli, Gabrielli and Pietronero 2012): F_r ∝ Σ_p M_rp Q_p,
//   Q_p ∝ 1 / Σ_r M_rp / F_r, each normalised to mean one after every step. *Caveat.* On matrices that are not nested the iteration
//   may drive some fitnesses toward zero (Servedio et al. 2018); the number of steps is a declared choice and is reported.
//
// Nothing here names a region or an industry. The measures predict, and do not identify causes (Hidalgo et al. 2018 review the
// out-of-sample evidence; Mealy et al. 2019 their meaning). Sources: E. S. Dunn, "A statistical and analytical technique for regional
// analysis", Pap. Proc. Reg. Sci. Assoc. 6 (1960) 97–112; J. M. Esteban-Marquillas, "A reinterpretation of shift-share analysis",
// Reg. Urban Econ. 2 (1972) 249–255; K. Frenken, F. Van Oort and T. Verburg, "Related variety, unrelated variety and regional economic
// growth", Reg. Stud. 41 (2007) 685–697; C. A. Hidalgo et al., "The product space conditions the development of nations", Science 317
// (2007) 482–487; C. A. Hidalgo and R. Hausmann, "The building blocks of economic complexity", PNAS 106 (2009) 10570–10575; P. Mealy,
// J. D. Farmer and A. Teytelboym, "Interpreting economic complexity", Sci. Adv. 5 (2019) eaau1705; A. Tacchella et al., "A new metrics
// for countries' fitness and products' complexity", Sci. Rep. 2 (2012) 723.

#pragma once

#include <calx/linalg/rectangular.hpp>
#include <tyche/networks/latent.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "inequality.hpp"
#include "specialisation.hpp"

namespace nomos::measurement {

// ---- Shift-share -----------------------------------------------------------------------------------------------------------------

struct ShiftShare {
  std::vector<double> national, mix, competitive;           // by industry
  std::vector<double> competitive_effect, allocation;       // Esteban-Marquillas's split (NaN where undefined)
  double change = 0, national_total = 0, mix_total = 0, competitive_total = 0, competitive_effect_total = 0, allocation_total = 0;
  bool homothetic = false;                                  // whether Esteban-Marquillas's split is defined (every e_i positive)
};

inline ShiftShare shift_share(const std::vector<double>& e0, const std::vector<double>& e1, const std::vector<double>& E0, const std::vector<double>& E1) {
  const std::size_t n = e0.size();
  detail::require(n > 0 && e1.size() == n && E0.size() == n && E1.size() == n, "employment by industry at two dates for the region and its reference");
  double r0 = 0, R0 = 0, R1 = 0;
  for (std::size_t i = 0; i < n; i++) {
    detail::require(e0[i] >= 0 && e1[i] >= 0 && E0[i] >= 0 && E1[i] >= 0, "employment is not negative");
    detail::require(E0[i] > 0 || e0[i] == 0, "an industry the region holds is held by its reference");
    r0 += e0[i], R0 += E0[i], R1 += E1[i];
  }
  detail::require(R0 > 0, "a reference with employment");
  const double G = R1 / R0 - 1;
  ShiftShare S;
  S.national.resize(n), S.mix.resize(n), S.competitive.resize(n), S.competitive_effect.assign(n, std::numeric_limits<double>::quiet_NaN());
  S.allocation.assign(n, std::numeric_limits<double>::quiet_NaN());
  S.homothetic = r0 > 0;
  for (std::size_t i = 0; i < n; i++) S.homothetic = S.homothetic && (e0[i] > 0 || E0[i] == 0);
  for (std::size_t i = 0; i < n; i++) {
    const double Gi = E0[i] > 0 ? E1[i] / E0[i] - 1 : 0;
    S.national[i] = e0[i] * G, S.mix[i] = e0[i] * (Gi - G), S.competitive[i] = (e1[i] - e0[i]) - e0[i] * Gi;
    S.change += e1[i] - e0[i], S.national_total += S.national[i], S.mix_total += S.mix[i], S.competitive_total += S.competitive[i];
    if (S.homothetic && E0[i] > 0) {
      const double gi = e1[i] / e0[i] - 1, star = r0 * E0[i] / R0;
      S.competitive_effect[i] = star * (gi - Gi), S.allocation[i] = (e0[i] - star) * (gi - Gi);
      S.competitive_effect_total += S.competitive_effect[i], S.allocation_total += S.allocation[i];
    } else if (S.homothetic) {
      S.competitive_effect[i] = 0, S.allocation[i] = 0;
    }
  }
  return S;
}

// ---- Variety ---------------------------------------------------------------------------------------------------------------------

struct Variety {
  double entropy = 0, unrelated = 0, related = 0;  // H = UV + RV
  std::map<std::string, double> within;            // H_g by group
};

// Amounts over fine classes and the group of each; logarithms to `base` (e by default).
inline Variety variety(const std::vector<double>& amounts, const std::vector<std::string>& group, double base = std::exp(1.0)) {
  detail::require(!amounts.empty() && group.size() == amounts.size(), "a group for every class");
  detail::require(base > 1, "a logarithm's base above one");
  double T = 0;
  for (double a : amounts) detail::require(a >= 0 && std::isfinite(a), "amounts are not negative"), T += a;
  detail::require(T > 0, "a positive total");
  const double lb = std::log(base);
  std::map<std::string, double> P;
  for (std::size_t i = 0; i < amounts.size(); i++) P[group[i]] += amounts[i] / T;
  Variety V;
  for (std::size_t i = 0; i < amounts.size(); i++) {
    const double p = amounts[i] / T;
    if (p > 0) V.entropy += p * std::log(1 / p) / lb, V.within[group[i]] += (p / P[group[i]]) * std::log(P[group[i]] / p) / lb;
  }
  for (const auto& [g, q] : P)
    if (q > 0) V.unrelated += q * std::log(1 / q) / lb, V.related += q * V.within[g];
  return V;
}

// ---- Relatedness and complexity --------------------------------------------------------------------------------------------------

struct Incidence {
  calx::linalg::Rectangular M;            // regions × activities, 0 or 1
  std::vector<double> diversity, ubiquity;
};

// M_rp = [RCA_rp ≥ threshold] from a table of amounts (regions × activities); NaN quotients count as absent.
inline Incidence incidence(const calx::linalg::Rectangular& x, double threshold = 1) {
  const auto q = location_quotients(x);
  Incidence I{calx::linalg::Rectangular(x.rows, x.cols), std::vector<double>(x.rows, 0), std::vector<double>(x.cols, 0)};
  for (std::size_t r = 0; r < x.rows; r++)
    for (std::size_t p = 0; p < x.cols; p++) {
      const double m = q(r, p) >= threshold ? 1.0 : 0.0;
      I.M(r, p) = m, I.diversity[r] += m, I.ubiquity[p] += m;
    }
  return I;
}

// An incidence given as a binary matrix (regions × activities).
inline Incidence from_binary(const calx::linalg::Rectangular& M) {
  Incidence I{M, std::vector<double>(M.rows, 0), std::vector<double>(M.cols, 0)};
  for (std::size_t r = 0; r < M.rows; r++)
    for (std::size_t p = 0; p < M.cols; p++) {
      detail::require(M(r, p) == 0 || M(r, p) == 1, "an incidence of zeros and ones");
      I.diversity[r] += M(r, p), I.ubiquity[p] += M(r, p);
    }
  return I;
}

// φ_pq = Σ_r M_rp M_rq / max(k_p, k_q); zero where an activity is held nowhere.
inline calx::linalg::Rectangular proximity(const Incidence& I) {
  const std::size_t P = I.M.cols;
  calx::linalg::Rectangular phi(P, P);
  for (std::size_t p = 0; p < P; p++)
    for (std::size_t q = 0; q < P; q++) {
      const double k = std::max(I.ubiquity[p], I.ubiquity[q]);
      double s = 0;
      for (std::size_t r = 0; r < I.M.rows; r++) s += I.M(r, p) * I.M(r, q);
      phi(p, q) = k > 0 ? s / k : 0;
    }
  return phi;
}

// ω_rp = Σ_{q ≠ p} φ_pq M_rq / Σ_{q ≠ p} φ_pq; NaN where p is near nothing.
inline calx::linalg::Rectangular density(const Incidence& I, const calx::linalg::Rectangular& phi) {
  const std::size_t R = I.M.rows, P = I.M.cols;
  calx::linalg::Rectangular w(R, P);
  for (std::size_t p = 0; p < P; p++) {
    double den = 0;
    for (std::size_t q = 0; q < P; q++)
      if (q != p) den += phi(p, q);
    for (std::size_t r = 0; r < R; r++) {
      double num = 0;
      for (std::size_t q = 0; q < P; q++)
        if (q != p) num += phi(p, q) * I.M(r, q);
      w(r, p) = den > 0 ? num / den : std::numeric_limits<double>::quiet_NaN();
    }
  }
  return w;
}

struct Complexity {
  std::vector<double> regions, activities;          // standardised indices (ECI, PCI)
  std::vector<double> regions_raw, activities_raw;  // the walks' second eigenvectors as tyche orders them
};

namespace detail {
inline std::vector<double> standardise(std::vector<double> v) {
  double m = 0, s = 0;
  for (double x : v) m += x / v.size();
  for (double x : v) s += (x - m) * (x - m) / v.size();
  s = std::sqrt(s);
  for (double& x : v) x = s > 0 ? (x - m) / s : 0;
  return v;
}
}  // namespace detail

// The complexity indices of regions and of activities, on the incidence's regions and activities that hold anything.
inline Complexity complexity(const Incidence& I) {
  detail::require(I.M.rows > 2 && I.M.cols > 2, "at least three regions and three activities");
  for (double d : I.diversity) detail::require(d > 0, "every region holds an activity (drop the empty rows)");
  for (double u : I.ubiquity) detail::require(u > 0, "every activity is held somewhere (drop the empty columns)");
  // A walk on an incidence that falls apart into blocks has the eigenvalue one more than once, and its "second" eigenvector only names
  // a block: refused.
  {
    const std::size_t R = I.M.rows, P = I.M.cols;
    std::vector<bool> seen_r(R, false), seen_p(P, false);
    std::vector<std::size_t> stack{0};
    seen_r[0] = true;
    while (!stack.empty()) {
      const std::size_t r = stack.back();
      stack.pop_back();
      for (std::size_t p = 0; p < P; p++)
        if (I.M(r, p) > 0 && !seen_p[p]) {
          seen_p[p] = true;
          for (std::size_t q = 0; q < R; q++)
            if (I.M(q, p) > 0 && !seen_r[q]) seen_r[q] = true, stack.push_back(q);
        }
    }
    for (std::size_t r = 0; r < R; r++) detail::require(seen_r[r], "the incidence falls apart into separate blocks; a complexity index would only name them");
  }
  Complexity C;
  C.regions_raw = tyche::networks::bipartite_ordering(I.M);
  // The activities' eigenvector of U⁻¹ Mᵀ D⁻¹ M is U⁻¹ Mᵀ applied to the regions' (Hidalgo and Hausmann 2009): each activity's mean
  // index over the regions that hold it, which fixes its sign with the regions'.
  C.activities_raw.assign(I.M.cols, 0.0);
  for (std::size_t p = 0; p < I.M.cols; p++) {
    for (std::size_t r = 0; r < I.M.rows; r++) C.activities_raw[p] += I.M(r, p) * C.regions_raw[r];
    C.activities_raw[p] /= I.ubiquity[p];
  }
  C.regions = detail::standardise(C.regions_raw), C.activities = detail::standardise(C.activities_raw);
  return C;
}

struct Fitness {
  std::vector<double> fitness, complexity;  // each with mean one
  int steps = 0;
  double change = 0;                        // the largest relative change in the last step
};

inline Fitness fitness_complexity(const Incidence& I, int steps = 200, double tolerance = 1e-12) {
  const std::size_t R = I.M.rows, P = I.M.cols;
  for (double d : I.diversity) detail::require(d > 0, "every region holds an activity");
  for (double u : I.ubiquity) detail::require(u > 0, "every activity is held somewhere");
  Fitness out{std::vector<double>(R, 1.0), std::vector<double>(P, 1.0), 0, 0};
  for (int s = 0; s < steps; s++) {
    std::vector<double> F(R, 0), Q(P, 0);
    for (std::size_t r = 0; r < R; r++)
      for (std::size_t p = 0; p < P; p++) F[r] += I.M(r, p) * out.complexity[p];
    for (std::size_t p = 0; p < P; p++) {
      double d = 0;
      for (std::size_t r = 0; r < R; r++)
        if (I.M(r, p) > 0) d += 1 / out.fitness[r];
      Q[p] = 1 / d;
    }
    double mf = 0, mq = 0;
    for (double v : F) mf += v / R;
    for (double v : Q) mq += v / P;
    double change = 0;
    for (std::size_t r = 0; r < R; r++) F[r] /= mf, change = std::max(change, std::fabs(F[r] - out.fitness[r]) / std::max(out.fitness[r], 1e-300));
    for (std::size_t p = 0; p < P; p++) Q[p] /= mq, change = std::max(change, std::fabs(Q[p] - out.complexity[p]) / std::max(out.complexity[p], 1e-300));
    out.fitness = F, out.complexity = Q, out.steps = s + 1, out.change = change;
    if (change < tolerance) break;
  }
  return out;
}

}  // namespace nomos::measurement
