# RESULT — rank 6: Lattès / isogeny layer maps in a PKM tower factor base

## Measurement
Setup per PROTOCOL.md: p in {211,521,1031}, s = 16/24/32, m=2, (k,d) in {(3,2),(2,3),(2,4)}; families generic, kummer, isogeny (Velu chains, rational 2-/3-isogenies; for d=4 two chained 2-isogenies per layer), lattes (d=4, x([4]) of an auxiliary curve). 3 seeds each, 90 systems, all solved by both engines (no timeouts). Full table: `analysis.txt` (raw: `raw/run1.json`, `raw/run2.json`, `raw/cert.txt`).

Correctness: S_3 vs brute force p=101, 397,953 triples, 0 false pos / 0 false neg; all Velu maps and Lattès map verified to send curve points to curve points; every system's planted decomposition is a root of all equations and survives Singular std (90/90).

Condensed (medians over 3 seeds; msolve F4 max degree / Singular std max degree):

| p, (k,d) | generic | kummer | isogeny | lattes |
|---|---|---|---|---|
| 211, (2,4) | 16 / 53 | 19 / 50 | 16 / 52 | 16 / 52 |
| 521, (2,4) | 24 / 67 | 27 / 66 | 24 / 68 | 24 / 68 |
| 1031, (2,4) | 32 / 81 | 35 / 82 | 32 / 84 | 32 / 84 |
| 211, (2,3) | 18 / 41 | 18 / 43 | 18 / 46 | - |
| 521, (2,3) | 26 / 57 | 26 / 59 | 26 / 62 | - |
| 1031, (2,3) | 34 / 70 | 34 / 75 | 34 / 78 | - |
| 211, (3,2) | 18 / 37 | 18 / 40 | 18 / 41 | - |
| 521, (3,2) | 26 / 56 | 26 / 54 | 26 / 58 | - |
| 1031, (3,2) | 34 / 69 | 34 / 72 | 34 / 74 | - |

Bézout = 4((d+1)^k s)^2 (0.26M-3M); vars 6 (k=2) or 8 (k=3). vdim 2-84, no systematic family effect visible. Pooled decomposition count / random-base control: generic 1.14, kummer 1.09, isogeny 1.05, lattes 1.01 (all ~1).

Twice-run diff: 90/90 systems identical (hash, both degrees, vdim, |F|, decomposition counts): 0 mismatches.

## Falsifier outcome
- msolve: structured = generic in 12/12 comparisons (degree = s + small constant, i.e. set by the top product polynomial). Falsifier holds (neutral). Kummer is +1..+3 in the d=4 rows.
- Singular: structured NOT lower. isogeny/lattes higher than generic in 10/12 comparisons (+3..+9, outside generic seed range), lower by 1 in 2/12 (within seed spread).

## Interpretation (separate)
Verdict: NEUTRAL on msolve, slight RAISES on Singular. No lowering on either engine; the arithmetic structure of Lattès/isogeny layers buys no degree and decomposition counts are random-like. Caveat: structured bases tend to have larger |F| and Singular degree rises mildly with it, so the Singular "raise" is partly confounded by |F| (fixed s, not fixed |F|); I did not disentangle it. The d_solv is dominated by the degree-s base polynomial, which is the same for all families, so this design has limited sensitivity.

Tier: T3 (single measured study, raw data, twice-run identical; engine comparison but not independently replicated). m=3 not run (T5).
Limitations: toy p <= 1031, 3 seeds, s fixed rather than |F| fixed, auxiliary curves independent of target, only k<=3.
Next step: fix |F| (not s) and use a larger-than-product base condition (e.g. x_k = 0 with more layers) to remove the s-dominated floor.

Files: PROTOCOL.md, tower.py, cert.py, one.py, driver.py, analyze.py, analysis.txt, raw/, fig_rank6_lattes_isogeny.png/.svg.
