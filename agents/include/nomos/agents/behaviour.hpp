// Behaviours of three origins behind one interface, kept apart (decision 0016).
//
// An agent's behaviour maps what it observes to what it does. Where the behaviour came from decides what a run built on
// it can establish, so it is declared with the behaviour and never lost:
//   rule      written down: a heuristic or a theory's decision rule, with its source (a paper, a declaration);
//   learned   fitted to a record or trained by reinforcement, with the data and the method of the fit;
//   proposed  put forward by a language model (or anyone not yet checked), with the model and the prompt.
// A run made only of written rules is computed from its declarations; one with learned behaviour is an estimate; one in
// which proposed behaviour decided anything is a proposal (the standings of heuresis's assessment: verified, computed,
// estimated, proposed). Learned policies show what an optimiser or a record would do, not what people do; a language
// model's behaviour beyond what has been validated is not evidence of human behaviour (the survey of evolutionary
// allocation, Part VI §6.5).
//
// Each behaviour declares what it observes and which acts it may take. Asked to decide, it is shown only the observations
// it declared, and an act it did not declare is refused, naming it: what a behaviour can see and do is part of the model
// (ODD+D's decision section, Müller et al. 2013, "Describing human decisions in agent-based models — ODD + D",
// Environ. Model. Softw. 48, 37–48).
//
// A tally counts the decisions and acts of each origin over a run, so that the share of a result that model-generated
// behaviour decided is reported beside it.
#pragma once

#include <calx/linalg/qr.hpp>
#include <calx/random/stream.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nomos::agents {

enum class Origin { rule, learned, proposed };

inline const char* name(Origin o) {
  switch (o) {
    case Origin::rule: return "rule";
    case Origin::learned: return "learned";
    default: return "proposed";
  }
}

inline Origin origin_of(const std::string& s) {
  for (Origin o : {Origin::rule, Origin::learned, Origin::proposed})
    if (s == name(o)) return o;
  throw std::invalid_argument("nomos::agents: there is no origin '" + s + "' (there are rule, learned and proposed)");
}

// The standing a run's results can have when the behaviour of this origin decided in it.
inline const char* standing(Origin o) {
  switch (o) {
    case Origin::rule: return "computed";
    case Origin::learned: return "estimated";
    default: return "proposed";
  }
}

struct Declared {
  std::string name;
  Origin origin = Origin::rule;
  std::vector<std::string> observes;  // the observations it may read
  std::vector<std::string> acts;      // the acts it may take
  std::string source;                 // a citation; the data and method of a fit; the model and prompt of a proposal
};

using Observation = std::map<std::string, double>;

// One act: its name, the agent or option it is directed to, and an amount (what each means is the model's).
struct Action {
  std::string act;
  std::size_t target = 0;
  double amount = 0;
};

using Decide = std::function<std::vector<Action>(const Observation&, calx::random::Stream&)>;

struct Behaviour {
  Declared declared;
  Decide decide;
};

struct Refused : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Decisions and acts by origin over a run.
struct Tally {
  std::array<long, 3> decisions{}, acts{};
  void add(Origin o, std::size_t n) {
    decisions[static_cast<std::size_t>(o)]++;
    acts[static_cast<std::size_t>(o)] += static_cast<long>(n);
  }
  long total() const { return decisions[0] + decisions[1] + decisions[2]; }
  // The share of decisions taken by learned or proposed behaviour.
  double model_share() const { return total() > 0 ? static_cast<double>(decisions[1] + decisions[2]) / static_cast<double>(total()) : 0; }
  // The weakest origin that decided anything gives the run its standing.
  const char* standing() const {
    if (decisions[2] > 0) return nomos::agents::standing(Origin::proposed);
    if (decisions[1] > 0) return nomos::agents::standing(Origin::learned);
    return nomos::agents::standing(Origin::rule);
  }
};

// Asks a behaviour to decide, showing it only what it declared it observes, and refusing an act it did not declare.
inline std::vector<Action> ask(const Behaviour& b, const Observation& world, calx::random::Stream& rng, Tally* tally = nullptr) {
  if (!b.decide) throw Refused("the behaviour '" + b.declared.name + "' has no way to decide");
  Observation seen;
  for (const auto& o : b.declared.observes) {
    const auto it = world.find(o);
    if (it == world.end()) throw Refused("the behaviour '" + b.declared.name + "' observes '" + o + "', which the model does not show");
    seen.emplace(o, it->second);
  }
  auto acts = b.decide(seen, rng);
  for (const auto& a : acts)
    if (std::find(b.declared.acts.begin(), b.declared.acts.end(), a.act) == b.declared.acts.end())
      throw Refused("the behaviour '" + b.declared.name + "' (" + name(b.declared.origin) + ") took the act '" + a.act + "', which it did not declare");
  if (tally) tally->add(b.declared.origin, acts.size());
  return acts;
}

// A policy learned as a linear function of the observations, its weights fitted elsewhere (`source` says how): the act's
// amount is max(0, w₀ + Σ_k w_k x_k), directed to `target`.
inline Behaviour linear_policy(std::string name, std::vector<std::string> observes, std::vector<double> weights, std::string act, std::size_t target, std::string source) {
  if (weights.size() != observes.size() + 1) throw std::invalid_argument("nomos::agents: a linear policy has an intercept and one weight per observation");
  Behaviour b;
  b.declared = {std::move(name), Origin::learned, observes, {act}, std::move(source)};
  b.decide = [observes, weights, act, target](const Observation& x, calx::random::Stream&) {
    double v = weights[0];
    for (std::size_t k = 0; k < observes.size(); k++) v += weights[k + 1] * x.at(observes[k]);
    return std::vector<Action>{{act, target, std::max(0.0, v)}};
  };
  return b;
}

// Least-squares weights for a linear policy from a record of observations and the amounts chosen, by calx's QR with
// column pivoting (the basic solution when the record does not determine every weight); a ridge λ ≥ 0 on the weights
// (not the intercept) is taken as √λ rows appended.
inline std::vector<double> fit_linear(const std::vector<std::vector<double>>& X, const std::vector<double>& y, double ridge = 0) {
  if (X.empty() || X.size() != y.size()) throw std::invalid_argument("nomos::agents: one amount for each observation");
  const std::size_t p = X[0].size() + 1, extra = ridge > 0 ? p - 1 : 0;
  calx::linalg::Rectangular A(X.size() + extra, p);
  std::vector<double> b(X.size() + extra, 0.0);
  for (std::size_t r = 0; r < X.size(); r++) {
    if (X[r].size() + 1 != p) throw std::invalid_argument("nomos::agents: every observation has the same length");
    A(r, 0) = 1;
    for (std::size_t k = 1; k < p; k++) A(r, k) = X[r][k - 1];
    b[r] = y[r];
  }
  for (std::size_t k = 1; k <= extra; k++) A(X.size() + k - 1, k) = std::sqrt(ridge);
  return calx::linalg::least_squares(A, b);
}

}  // namespace nomos::agents
