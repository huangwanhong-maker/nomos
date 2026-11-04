# 0005. Kirman's Herding

- **Status:** Accepted, provisionally (two views; more views and agents on a network planned)
- **Date:** 2026-11-05
- **Decided by:**
  - **The request.** The user asked to visualise "such algorithms in economics state perception and decision making", among them the agent-based models of decision 0004.
  - **What Claude did.** It wrote the model, its check and this record; the infrastructure's Economy Atlas draws it ([ADR 0074](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0074-perception-and-decision-in-the-economy-atlas.md)).

## 1. Context

`agents` held one model of the social sciences, Schelling's sorting (decision 0004). The economics of herding needed a second model, whose stationary law is known exactly: Kirman's (1993) ants, the recruitment of opinion in a market. N agents each hold one of two views. An agent changes its view alone at rate ε, or is recruited by another of the other view at rate δ per such agent, a share of N − 1. The number k holding view 1 is then a birth–death chain whose stationary law is a product of the ratios of its rates. Below ε = δ/(N − 1) that law is bimodal: although no agent prefers either view, the market spends most of its time near one or the other and switches between them suddenly.

## 2. Decision

**`agents/choice.hpp`, `Kirman`:**
- a model for the schedules of `engine.hpp`: the rate of agent i is ε + δ · (others of the other view)/(N − 1); a change of view changes every rate, so every agent is affected;
- `ones(x)`, the number holding view 1;
- `stationary(N)`, the exact law of that number, by detailed balance of the birth–death chain;
- `describe()`, its ODD+D description.

## 3. Alternatives Considered

- **The chain of k alone, without agents.** Rejected: the law of k is the check; the model is the agents, so that it runs on the same schedules as every other model, with draws keyed by agent and time.
- **Checking the simulated law of k directly.** Rejected after its check: which end the market holds mixes slowly, a switch every few thousand units of time, so a run of 200,000 units leaves the two ends unequal by chance. The shape of the herd, min(k, N − k), mixes fast, and is checked against the folded exact law.

## 4. Consequences

- **Checks.** `check_nomos_agents` gains two (16 in all), with N = 20, ε = 0.008 and δ = 0.3 (below δ/(N − 1) = 0.0158):
  - the exact law is bimodal, its ends the likeliest;
  - the herd's shape, sampled once a unit of time for 200,000 units, lies within total variation 0.03 of the folded exact law (0.013).
- **Its first use.** The Economy Atlas's Decisions view draws a market of 40 agents and its law over a run of 200,000 units against the exact one.
- **Planned.** More than two views; recruitment along a network rather than by everyone.
