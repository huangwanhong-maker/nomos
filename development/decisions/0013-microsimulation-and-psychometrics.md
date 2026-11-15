# 0013. Microsimulation of Taxes and Benefits, and Psychometrics

- **Status:** Accepted, provisionally (2026-11-15).
- **Date:** 2026-11-15
- **Decided by:** the user, who placed quantitative analysis among the layers beneath policy making (the infrastructure's ADR 0116, its fourth strand); and Claude.

## 1. Context

Policy making asks what a reform of taxes and benefits would do to households, and social analysis measures traits that no one observes directly. Neither was in nomos, whose README planned `psychometrics`.

## 2. Decision

- `microsimulation`: static tax-benefit microsimulation of declared rules, a reform compared with its baseline, and incentives found by applying the rules again.
- `psychometrics`: classical test theory, and item response theory by marginal maximum likelihood.

General factor and latent-class models went to tyche (its decision 0025). nomos keeps the scales and items of the social sciences.

## 3. Consequences

`check_nomos_microsimulation` (16), `check_nomos_psychometrics` (17).
