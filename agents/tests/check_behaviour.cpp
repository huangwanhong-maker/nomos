// Checks of nomos/agents/behaviour: a behaviour is shown only what it declared it observes; an act it did not declare is
// refused, and so is an observation the model does not show; the tally gives a run the standing of the weakest origin
// that decided in it and the share model-generated behaviour decided; a learned linear policy is fitted exactly from a
// record it can reproduce, clamped at zero, and the fit's basic solution still reproduces a record with a redundant column.

#include <nomos/agents/behaviour.hpp>

#include <calx/verification/check.hpp>

#include <cmath>
#include <string>

using namespace nomos::agents;

int main() {
  calx::verification::Suite s;
  calx::random::Stream rng(7, 1, 0, 0);
  {
    Behaviour b;
    b.declared = {"spend what was received", Origin::rule, {"received"}, {"pay"}, "a declared rule"};
    std::size_t seen = 0;
    bool only = true;
    b.decide = [&](const Observation& o, calx::random::Stream&) {
      seen = o.size();
      only = o.count("received") == 1 && o.count("held") == 0;
      return std::vector<Action>{{"pay", 3, o.at("received")}};
    };
    Tally t;
    const auto a = ask(b, {{"received", 5}, {"held", 40}}, rng, &t);
    s.require("a behaviour is shown only the observations it declared", seen == 1 && only);
    s.require("... and its act passes", a.size() == 1 && a[0].act == "pay" && a[0].target == 3 && a[0].amount == 5);
    s.require("a run of written rules alone is computed", std::string(t.standing()) == "computed" && t.model_share() == 0);
    bool refused = false;
    try {
      ask(b, {{"held", 40}}, rng);
    } catch (const Refused& e) {
      refused = std::string(e.what()).find("'received'") != std::string::npos;
    }
    s.require("an observation the model does not show is refused, naming it", refused);
  }
  {
    Behaviour b;
    b.declared = {"a model's suggestion", Origin::proposed, {"held"}, {"pay"}, "a language model, prompt kept"};
    b.decide = [](const Observation&, calx::random::Stream&) { return std::vector<Action>{{"borrow", 0, 10}}; };
    bool refused = false;
    try {
      ask(b, {{"held", 1}}, rng);
    } catch (const Refused& e) {
      refused = std::string(e.what()).find("'borrow'") != std::string::npos && std::string(e.what()).find("proposed") != std::string::npos;
    }
    s.require("an act a behaviour did not declare is refused, naming the act and the origin", refused);
  }
  {
    Tally t;
    t.add(Origin::rule, 1), t.add(Origin::rule, 2), t.add(Origin::learned, 1);
    s.require("one learned decision makes the run an estimate", std::string(t.standing()) == "estimated");
    s.check("... and model-generated behaviour decided a third of the decisions", t.model_share() - 1.0 / 3, 1e-15);
    t.add(Origin::proposed, 0);
    s.require("one proposed decision makes it a proposal", std::string(t.standing()) == "proposed");
    s.require("the acts of each origin are counted apart", t.acts[0] == 3 && t.acts[1] == 1 && t.acts[2] == 0);
  }
  {
    // A record made by amount = 2 + 0.5 x − 0.25 y: the fit recovers it exactly.
    std::vector<std::vector<double>> X;
    std::vector<double> y;
    for (int i = 0; i < 12; i++) {
      const double a = i, b = (i * 7) % 5;
      X.push_back({a, b});
      y.push_back(2 + 0.5 * a - 0.25 * b);
    }
    const auto w = fit_linear(X, y);
    s.check("a linear policy fitted from a record it can reproduce: intercept", w[0] - 2, 1e-12);
    s.check("... and weights", std::fabs(w[1] - 0.5) + std::fabs(w[2] + 0.25), 1e-12);
    const auto p = linear_policy("fitted spending", {"x", "y"}, w, "pay", 1, "least squares on twelve records");
    s.require("a learned policy declares its origin and its fit", p.declared.origin == Origin::learned && !p.declared.source.empty());
    const auto a = ask(p, {{"x", 1}, {"y", 20}}, rng);
    s.check("... and is clamped at zero where the line falls below it", a[0].amount, 0.0);
    // A redundant column (twice the first): the basic solution still reproduces the record.
    std::vector<std::vector<double>> Z;
    for (const auto& r : X) Z.push_back({r[0], r[1], 2 * r[0]});
    const auto v = fit_linear(Z, y);
    double worst = 0;
    for (std::size_t r = 0; r < Z.size(); r++) worst = std::max(worst, std::fabs(v[0] + v[1] * Z[r][0] + v[2] * Z[r][1] + v[3] * Z[r][2] - y[r]));
    s.check("a record with a redundant column: the basic solution reproduces it", worst, 1e-10);
    const auto r = fit_linear(X, y, 1e6);
    s.require("a large ridge shrinks the weights, not the intercept", std::fabs(r[1]) < 0.01 && std::fabs(r[2]) < 0.01);
  }
  s.require("origins read back from their names", origin_of("rule") == Origin::rule && origin_of("learned") == Origin::learned && origin_of("proposed") == Origin::proposed);
  return s.done();
}
