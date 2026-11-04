# 0004. Agent-Based Models

- **Status:** Accepted, provisionally (discrete-time and Gillespie schedules; calibration, social force and agents that plan are planned)
- **Date:** 2026-11-04
- **Decided by:**
  - **The request.** The user asked for agent-based models and POMDPs "for modern urban planning", reusable and kept apart from the applications. Asked where they belong, they chose nomos's planned `agents` for the first and tyche's new `decision` for the second (tyche decision 0009).
  - **What Claude did.** It wrote the module and this record. The design is the infrastructure's [ADR 0071](https://github.com/huangwanhong-maker/The-Science-Infrastructures/blob/main/development/decisions/0071-agents-and-decisions-for-urban-planning.md).

## 1. Context

**The plan.** Decision 0003 planned `agents`: "agent-based models: agents, schedules and interactions on calx's random streams". The infrastructure's urban-science domain now needs it for households choosing dwellings, walkers on the street network, and owners and occupiers changing a district's uses. history, economics and the planned domains of the social sciences need it too.

**Why it is nomos's.** An agent-based model names no society: an agent is numbered, a rule is a kernel, and a space says who is next to whom. It animates the agents whose static problems nomos already holds: their games, their matchings and their choices.

## 2. Decision

**A new module, `agents`**, header-only, on calx's `random` and `scheduler`.

1. **Schedules** (`engine.hpp`). The schedule is part of the model, never hidden:
   - synchronous: a product of kernels, order-free, run on calx's pool when given one;
   - random sequential: a composition, in an order drawn each tick;
   - Gillespie's direct method in continuous time, on a Fenwick tree of rates, recomputing only the rates a firing affects.
2. **Draws by identity.** Agent i's draws at tick t come from the stream (seed, "agnt", i, t); the order of a tick, and Gillespie's draws, from their own streams. A synchronous run is therefore the same on any number of threads. Its event log is kept per agent and joined in the agents' order.
3. **The event log.** When, which agent, which rule, what it touched. A domain reads it into its record of events, so that a simulated history is analysed as an observed one is.
4. **Ensembles.** Replicates by seed, with means and standard errors.
5. **The ODD+D description** (Grimm et al. 2006, 2010, 2020; Müller et al. 2013) as a structure each model fills.
6. **Spaces** (`space.hpp`). Lattices (Moore and von Neumann, bounded or toroidal), graphs from edges, and the plane binned for radius queries.
7. **Choices** (`choice.hpp`).
   - The multinomial logit (McFadden 1974), computed stably.
   - Schelling's model of segregation (Schelling 1971), as the canonical model and check, and as the base of the infrastructure's model of residential choice.

## 3. Alternatives Considered

- **A borrowed framework** (Mesa, NetLogo, Repast, MASON, Agents.jl). Rejected for the engine. The infrastructure computes in C++ on its own streams, so that a run is reproducible by identity, and its checks hold the engine to exact laws. The frameworks stay references for models and for their conventions.
- **A class hierarchy of agents.** Rejected. A model is a type with `update` (discrete time) or `rate` and `fire` (continuous time), so that its state can be any structure: a vector of kinds, a district's cells, a network of walkers.
- **Synchronous updating as the default.** Rejected: no schedule is the default. Each model states its own, since synchronous updating can create artefacts the asynchronous model lacks (Huberman and Glance 1993).

## 4. Consequences

- **Checks.** `check_nomos_agents` (14):
  - the rate tree against a linear search;
  - Gillespie's SIS of 5 agents against physis's master equation, exactly (Pearson's χ² 1.7 on 5 degrees of freedom);
  - Kurtz's mean field on 1000 agents (0.4643 ± 0.0040 against the logistic 0.4659);
  - the voter model's martingale (0.239 ± 0.007 against 0.25);
  - the same synchronous states and log on one thread and on four;
  - Schelling's sorting: 0.50 to 0.74 like neighbours at a threshold of ⅓, and 0.87 at ½;
  - logit choice;
  - the plane's queries.
- **calx** gained `random/categorical.hpp` for the agents' choices and tyche's particles.
- **The tests depend on physis's `statistical_mechanics`**, through the sibling the root CMake already adds; the module itself does not.
- **Planned.**
  - Agents that plan, through tyche's `decision`.
  - Calibration by indirect inference and ABC, tyche's to write, with the engine as their simulator.
  - Social force for pedestrians.
  - Births and deaths with stable identities beyond an index.
