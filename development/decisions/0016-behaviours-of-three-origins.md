# 0016. Behaviours of Three Origins, Kept Apart

- **Status:** Accepted, provisionally (2026-11-21).
- **Date:** 2026-11-21
- **Decided by:** Claude, for the infrastructure's agent economies whose every payment is posted (its ADR 0139, Stage 7 of evolutionary allocation), after the survey of agent simulation of economic circulation (its report *Evolutionary heterogeneous resource allocation*, Part VI §6.1 and §6.5).

## 1. Context

An agent's behaviour may be a written rule, a policy learned from a record or by reinforcement, or a proposal of a language model. What a simulation can establish depends on which: a run of written rules is computed from its declarations; learned behaviour shows what an optimiser or a record would do, not what people do; a model's behaviour beyond what has been validated is not evidence of human behaviour. `agents` held schedules, draws, logs and choices, and nothing that kept the origin of a behaviour with it.

## 2. Decision

`agents/behaviour.hpp`:

- `Origin{rule, learned, proposed}` and the standing each gives a run (computed, estimated, proposed: the standings of heuresis's assessment);
- `Declared{name, origin, observes, acts, source}`: what a behaviour may observe and do, and where it came from (a citation, the data and method of a fit, the model and prompt of a proposal);
- `Behaviour{declared, decide}`, `Action{act, target, amount}`, `Observation` (named numbers);
- `ask`: the behaviour is shown only the observations it declared, and an act it did not declare is refused, naming the act and the origin; an observation the model does not show is refused, naming it (ODD+D's decision section, Müller et al. 2013);
- `Tally`: decisions and acts by origin, the share decided by learned or proposed behaviour, and the run's standing, the weakest origin that decided in it;
- `linear_policy` and `fit_linear`: a learned policy as a linear function of observations, clamped at zero, its weights fitted by calx's QR with column pivoting (with a ridge on the weights), its fit recorded in its source.

## 3. Alternatives Considered

- **Origin as a comment.** Lost at the first composition; a tally cannot be made of it.
- **The economics domain keeping it.** Nothing in it names an economy: the origin of a behaviour matters to any agent model.

## 4. Consequences

`check_nomos_behaviour` (16 checks): only declared observations shown; an undeclared act and an unshown observation refused by name; standings and shares of the tally; a linear policy recovered exactly from a record it can reproduce, clamped at zero, the basic solution with a redundant column, a ridge shrinking the weights and not the intercept. The economics domain's `ledger` runs its agents through `ask` and reports the tally with every reading.
