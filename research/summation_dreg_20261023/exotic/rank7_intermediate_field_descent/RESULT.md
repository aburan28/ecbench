# RESULT — rank 7: intermediate-field (GF(2^k)) descent vs F_2 descent
Date 2026-10-23. Code: `r7.py` (all configs), `prime_n.py`, `analyze.py`. Raw: `raw/run1.json`, `raw/run2.json`,
`raw/prime{1,2}.json`, `raw/analysis.txt`. Figure: `fig_rank7_intermediate_field_descent.{png,svg}`.

## Verdict
**Neutral / harness check passed (T3).** The harness reproduces the Diem/Gaudry bounded-degree behaviour for the
GF(q)-descent: descended equation degrees ([4,2] for S_3, [12,9,9] for S_4), Bezout (8, 972) and vdim without field
equations (8, 384) are identical for every k tested. Falsifier NOT triggered. The "F_2-descent degree grows with k*m"
expectation holds only partly (see table): it holds for Singular (both metrics), not for msolve F4 (flat 3 / 6).
No claim about prime n beyond the structural fact below.

## Measurement
Systems: A0 = GF(q)-descent, Singular, ring `(2^k,a)`, no field equations; A1 = A0 + x_i^q+x_i; B = F_2-descent
(k*m boolean unknowns + c^2=c) with the factor base the embedded subfield (basis alpha^j), msolve F4 and Singular std.
msolve cannot do GF(2^k), k>1 (prime fields only), so A has no msolve column.
Tower: GF(2^n)=GF(q)[t]/g, g = min poly of the harness generator z over GF(q) (computed with `harness/gf2.py`
arithmetic; GF(q) uses the harness IRRED[k], which equals Singular's Conway poly for k=2,3,4, checked).
Singular degree is given as `[k:d]` (harness convention) / `label` (largest degree label in the std protocol); the two
disagree in size, so both are shown.

| n | k | m | curve | #planted-class sols (brute force = A1 vdim = B vdim) | A0: eq degs / Bezout / vdim | A0 Singular [k:d]/label | B vars | B msolve | B Singular [k:d]/label |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 2 | 2 | a6 generic | 4 | [4,2] / 8 / 8 | 2/6 | 4 | 3 | 2/5 |
| 6 | 3 | 2 | Koblitz | 4 | [4,2] / 8 / 8 | 2/6 | 6 | 3 | 3/5 |
| 6 | 3 | 2 | generic | 2 | [4,2] / 8 / 8 | 2/6 | 6 | 3 | 3/5 |
| 8 | 4 | 2 | Koblitz | 4 | [4,2] / 8 / 8 | 2/6 | 8 | 3 | 3/7 |
| 8 | 4 | 2 | generic | 2 | [4,2] / 8 / 8 | 2/6 | 8 | 3 | 3/7 |
| 6 | 2 | 3 | generic | 6 | [12,9,9] / 972 / 384 | 2/28 | 6 | 6 | 3/9 |
| 9 | 3 | 3 | generic | 6 | [12,9,9] / 972 / 384 | 2/27 | 9 | 6 | 4/17 |
| 12 | 4 | 3 | generic | 6 | [12,9,9] / 972 / 384 | 2/27 | 12 | 6 | 4/33 |

Generic = a2=1, a6 = z^2+z+1 (outside the subfield). Koblitz m=3 and (4,2,2) Koblitz: no planted decomposition exists
(see Finding 2); recorded as SKIP in the raw JSON.
Correctness gate: for all 8 systems, number of ordered x-tuples with a signed point-sum = +-R by brute-force group law
== zeros of the harness S_{m+1} formula == vdim(A1) == vdim(B); planted tuple is a root of A0, A1, B (0 FP / 0 FN in counts).
Twice-run diff (system hashes, degrees, vdim, counts; run1 vs run2 with different hash seeds): **identical** (after
fixing one nondeterminism in my own export, set-iteration order of monomials; the first attempt differed only in B_hash).

Prime n (`raw/prime1.json`): n=11,13 have divisors {1,n}; #{x: x^(2^d)=x} = 2 (d=1) and 2^n (d=n): the only subfield is
F_2 itself. n=17,131,233,283 prime as well. The only factor-base lever is an arbitrary l-dimensional F_2-subspace
(B-type): l=3,m=2 at n=11 and 13: 6 vars, msolve 3, Singular [k:d] 3, vdim 2, planted recovered.

## Findings (measurement)
1. A0 degrees/Bezout/vdim independent of k (m=2: 8; m=3: 972 / 384), also independent of curve (Koblitz vs generic).
2. Structural: a curve defined over GF(q) (e.g. Koblitz) has E(GF(q)) closed under addition; for odd m the points with
   x in GF(q) have y in GF(q) only, so every sum has x_R in GF(q): descended equations become dependent (no
   isolated problem). For m=2 the twist points (y in GF(q^2)) make Koblitz work. So the Diem regime for m=3 needs a
   curve NOT defined over GF(q).
3. B: msolve flat (3 for m=2, 6 for m=3); Singular [k:d] 2→3→3 (m=2), 3→4→4 (m=3); Singular label 5,5,7 and 9,17,33.
   vdim equals solution count (tiny) everywhere; A0's vdim (8, 384) counts all closure solutions, A1's the rational ones.

## Interpretation / limitations
- Harness reproduces the theorem's bounded-degree statement for the algebraic data (equation degrees, Bezout, vdim).
  A0 Singular solving degree (label 6, 27-28) is likewise k-flat, near the Macaulay-type bound.
- The F_2 side's growth depends on engine and metric; with msolve it does not grow in this range. Honest reading: the
  descent-to-F_2 penalty is visible in Singular only; k <= 4 is too small to separate trends from field-equation/
  non-homogeneity artifacts of Singular's label metric.
- GF(q) A0 is Singular-only (no msolve). A0 degrees are for the system without field equations (the Diem regime);
  with field equations x^q+x the degree would be forced up to ~q.
- Prime n: no intermediate field, so this lever is structurally unavailable; this experiment gives no information
  on how hard real targets (131/233/283) are. Toy sizes only.
- Tier T3 (single measured study, raw data, twice-run identical, count gate T1-style on toy curves).
- Next step: larger k (5,6,7) at m=2 with msolve-capable prime-field analogue, and an engine-independent degree
  (e.g. Macaulay-matrix first-fall) for the GF(q) side.
