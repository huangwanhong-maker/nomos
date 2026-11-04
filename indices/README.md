# indices: Index Numbers

*Status: draft (2026-10-25), [decision 0003](../development/decisions/0003-the-first-modules.md).*

## 1. Concepts

An index number compares the prices (or quantities) of the same N items in two situations, periods or places. A price index is a function P(p⁰, p¹, q⁰, q¹) of the prices and quantities of both, written in this order as the axiomatic literature writes it (Diewert 1992; ILO et al. 2004, ch. 16).

- **Fixed baskets.** Laspeyres values the basket of situation 0, Paasche that of situation 1; Lowe the basket of a reference period, which is how most consumer price indices are compiled, and Young the expenditure shares of one.
- **Symmetric and superlative indices.** Fisher's ideal index is the geometric mean of Laspeyres and Paasche (Fisher 1922); Törnqvist's averages the logarithmic price changes with the mean of the two periods' shares (Törnqvist 1936); Walsh's values the geometric mean of the baskets (Walsh 1901); Marshall–Edgeworth's their arithmetic mean. Fisher, Törnqvist and Walsh are superlative: exact for flexible aggregators, and close to one another (Diewert 1976, 1978).
- **Quantity indices** exchange the roles of prices and quantities.
- **Elementary indices** compare price quotes without quantities (ILO et al. 2004, ch. 20): Jevons (the geometric mean of price relatives), Dutot (the ratio of mean prices), Carli (their arithmetic mean), the harmonic mean, and CSWD, the geometric mean of Carli and harmonic (Carruthers, Sellwood and Ward 1980; Dalén 1992).
- **The axiomatic approach** judges a formula by the tests it passes. Identity, proportionality, commensurability and the mean value test are passed by all six formulas with quantities; time reversal by Fisher, Törnqvist, Walsh and Marshall–Edgeworth; factor reversal (P·Q = V, the value ratio) by Fisher alone. Among elementary indices Jevons is transitive, Carli fails time reversal upward, Dutot fails commensurability.
- **Series.** A direct series compares every period with the base; a chained series multiplies the links between consecutive periods. They agree when the formula is transitive (circular), which no formula with quantities is in general; when prices and quantities return to where they started a chained index may not return to one: chain drift (Szulc 1983).

## 2. Content

| Header | Content |
|---|---|
| `bilateral.hpp` | `enum class Formula {laspeyres, paasche, fisher, tornqvist, walsh, marshall_edgeworth}`; `laspeyres`, `paasche`, `fisher`, `tornqvist`, `walsh`, `marshall_edgeworth`, each `(p0, p1, q0, q1)`; `price_index(f, p0, p1, q0, q1)`, `quantity_index(f, p0, p1, q0, q1)`, `value_ratio(p0, p1, q0, q1)`; `lowe(p0, p1, qb)`, `young(p0, p1, sb)`; elementary `jevons(p0, p1)`, `dutot`, `carli`, `harmonic`, `cswd` |
| `series.hpp` | `Periods` (one vector of items per period); `chain_links(f, p, q)`, `direct_series(f, p, q, base)`, `chained_series(f, p, q, base)` |

Prices must be positive, quantities nonnegative; violations throw `std::invalid_argument`.

## 3. Checks

`tests/check_indices.cpp` (35 checks):

- the artificial data set of the CPI Manual (ILO et al. 2004, ch. 19; six commodities, five periods): fixed-base Laspeyres (1, 1.4200, 1.3450, 1.3550, 1.4400), Paasche (1, 1.3824, 1.2031, 1.0209, 0.7968), Fisher (1, 1.4011, 1.2721, 1.1762, 1.0712) and Törnqvist (1, 1.4052, 1.2890, 1.2268, 1.2477), and chained Laspeyres (1, 1.4200, 1.3646, 1.3351, 1.3306), to the four printed decimals; chained Fisher and Törnqvist within 0.01 of each other and within 3·10⁻⁴ in period 5; rebasing;
- on random data: identity, proportionality, invariance to proportional changes in current quantities (failed by Marshall–Edgeworth), commensurability (failed by Dutot), the mean value test; time reversal (and its failure by Laspeyres, whose product with the reversed index is P_L/P_P); factor reversal (Fisher passes, Törnqvist fails); Laspeyres' quantity index; Lowe and Young with period 0's basket and shares are Laspeyres;
- elementary: Jevons transitive and reversible, Carli biased upward, harmonic the reciprocal of reversed Carli, CSWD reversible, harmonic ≤ Jevons ≤ Carli;
- bouncing prices: every direct index returns to one, chained Laspeyres drifts, chained Fisher less, chained Jevons not at all.

## 4. Dependencies

None in the headers. In the checks, calx `random` and `verification`.

## 5. References

- Carruthers, A. G., Sellwood, D. J. and Ward, P. W. (1980), Recent developments in the retail prices index, *The Statistician* 29, 1–32.
- Dalén, J. (1992), Computing elementary aggregates in the Swedish consumer price index, *J. Off. Stat.* 8, 129–147.
- Diewert, W. E. (1976), Exact and superlative index numbers, *J. Econometrics* 4, 115–145; (1978), Superlative index numbers and consistency in aggregation, *Econometrica* 46, 883–900; (1992), Fisher ideal output, input, and productivity indexes revisited, *J. Prod. Anal.* 3, 211–248.
- Eichhorn, W. and Voeller, J. (1976), *Theory of the Price Index*, Springer.
- Fisher, I. (1922), *The Making of Index Numbers*, Houghton Mifflin.
- ILO, IMF, OECD, UNECE, Eurostat and World Bank (2004), *Consumer Price Index Manual: Theory and Practice*, International Labour Office, Geneva.
- Lowe, J. (1823), *The Present State of England in Regard to Agriculture, Trade and Finance*, 2nd ed., London.
- Szulc, B. J. (1983), Linking price index numbers, in W. E. Diewert and C. Montmarquette (eds), *Price Level Measurement*, Statistics Canada.
- Törnqvist, L. (1936), The Bank of Finland's consumption price index, *Bank of Finland Monthly Bulletin* 10, 1–8.
- Walsh, C. M. (1901), *The Measurement of General Exchange Value*, Macmillan.
- Young, A. (1812), *An Inquiry into the Progressive Value of Money in England*, London.

## 6. Open Questions

1. **Multilateral indices** for comparisons across many places (GEKS, Geary–Khamis, CCD) and their transitivity.
2. **Hedonic and unit-value indices**, when a domain brings item-level data.
