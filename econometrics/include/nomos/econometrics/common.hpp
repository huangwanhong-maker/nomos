// What nomos's econometrics shares with tyche: the absorption of fixed effects (tyche::panel) and the pieces of robust covariances
// (tyche::regression's sandwich), named here as the estimators that stay in nomos use them. Panels, instruments and differences in
// differences moved in peira's stage 2 (ADR 0082; nomos's decision 0007): panels to tyche::panel, instruments and differences in
// differences to peira::quasi.

#pragma once

#include <tyche/panel/absorb.hpp>
#include <tyche/regression/linear.hpp>
#include <tyche/regression/sandwich.hpp>

#include <calx/linalg/rectangular.hpp>

#include <stdexcept>
#include <string>

namespace nomos::econometrics {

namespace tr = tyche::regression;
using calx::linalg::Matrix;
using calx::linalg::Rectangular;
using tyche::panel::AbsorbOptions;
using tyche::panel::demean;
using tyche::panel::Factor;
using tyche::panel::factor;
using tyche::panel::recover_effects;

namespace detail {
inline void require(bool ok, const char* what) {
  if (!ok) throw std::invalid_argument(std::string("nomos::econometrics: ") + what);
}
using tyche::regression::join_columns;
using tyche::regression::sandwich;
using tyche::regression::sandwich_meat;
using tyche::regression::solve_spd;
}  // namespace detail

}  // namespace nomos::econometrics
