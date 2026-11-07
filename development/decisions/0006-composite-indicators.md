# 0006. Composite Indicators

- **Status:** Accepted, provisionally (2026-11-07).
- **Decided by:** Claude, for the Urban Science Workbench's *A composite of place* (the infrastructure's ADR 0078), after GeoE3's method; none of GeoE3's code (GPL-3.0) is used.

## 1. Decision

`measurement/composite.hpp` holds factors on scales of their own brought to one scale and combined by weights: each scaled between two percentiles (a robust min-max), inverted if more of it is worse, and their weighted arithmetic mean, a factor missing for a unit left out of its mean and the weights renormalised. After the OECD and JRC's *Handbook on Constructing Composite Indicators* (2008), §5 and §6.

It belongs here rather than in tyche since composites are a measure of the social sciences, as inequality is.

## 2. Caveat

A linear aggregation is compensatory: much of one factor makes up for little of another. What the weights mean is the analyst's to declare. The workbench's tool says that the index is a description by one's own weights, and not a measure.

## 3. Verification

`check_nomos_composite`: closed forms for the scaling, the percentiles, the mean and missing values, invariance to doubling every weight, and refusals.
