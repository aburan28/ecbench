# RESULT -- rank 2: order-6 automorphism quotient on j=0 prime curves (m=2)

## Verdict
**LOWERS on both engines (tier T3), but the mechanism is the factor-base-equation degree (x^d -> u^(d/3)), not a Bezout/S3 gain.**
Falsifier (PROTOCOL.md, unchanged) NOT triggered: "quotient >= plain on both engines" held in only 8/58 instances (need >= 2/3).
Quotient < plain: msolve 46/58, Singular 50/58. Bezout (B_ci) ratio >= 2x in 26/58 and exactly 1.0x for all d >= 36.
m=3 (S_3 chain) not run (T5).

## Measurements (58 instances = 6 primes p in {61,97,193,1009,4093,65521}, 16 values of d (3|d, d|p-1), 2 curves/targets each; run1 == run2)
Degrees are means over instances with the same d (full per-instance data: raw/summary_table.json, run1/results.json).
Vars = 2 and max eq degree = d (plain) / 12 (quotient, d>=36) for both; Singular and msolve 'degree' defined in PROTOCOL Amendments 1-2.

| d | n | msolve plain | msolve quot | Singular plain | Singular quot | B_ci plain/quot |
|---|---|---|---|---|---|---|
| 6 | 10 | 7 | 12 | 14.2 | 14.4 | 24/4 = 6x |
| 9 | 2 | 10 | 12 | 20 | 17 | 36/9 = 4x |
| 12 | 10 | 13 | 12 | 25.8 | 18.2 | 3x |
| 18 | 4 | 19 | 12 | 38.5 | 22 | 2x |
| 24 | 4 | 25 | 14 | 49 | 24.8 | 1.5x |
| 36 | 4 | 37 | 18 | 74 | 33.8 | 1.0x |
| 48 | 4 | 49 | 22 | 93.8 | 38 | 1.0x |
| 96 | 4 | 97 | 38 | 181.5 | 66 | 1.0x |
| 126 | 2 | 127 | 48 | 253 | 92 | 1.0x |
(other d: 30,33,45,60,63,66,90 same trend.) Crossover: at d=6 (and msolve d=9) the quotient is NOT lower (its fixed degree-12 orbit polynomial dominates); from d=12 up it is lower, ratio -> ~2.0x (msolve) / ~2.7x (Singular) at d=96-126.

Bezout: B_all (product of the three degrees) plain/quot = 3x for every d>=... (4d^2 vs 4d^2/3: (d*d*4)/((d/3)^2*12)=3). B_ci (min pair, honest complete-intersection bound) = 6x at d=6, falling to 1x at d>=36 because the S3 hypersurface in u has degree 12 = 3*4, so the quotient never reduces the S3-vs-fieldeq intersection count.
Solutions / vdim (Singular; radical vdim == brute-force count in all 116 systems): total plain 620 vs quotient 534 across instances; one quotient solution lifts to 1 (not 9) plain solution in general, since the target x_R is fixed (the 27-factor target-orbit variant equals Q^3 up to a unit: no new information). So cost per decomposition improves by about the degree ratio only (~1.16 plain solutions per quotient solution).

## Correctness gates (all passed)
- Symbolic: S3(wx1,wx2,wx3)=w*S3 (residue 0); orbit product Q is a polynomial in u1=x1^3,u2=x2^3, degree 6 in each, total 12, 49 terms (raw/Q_u_poly.json). 27-factor variant = unit*Q^3 (numeric, 5 primes).
- S3 vs real group law: 204 curves (p in {7..61}), 0 false positives / 0 false negatives, semi-invariance holds (raw/cert_S3.json).
- Every instance: quotient solutions lifted by 9 rotations equal the plain solution set; planted (x1,x2) recovered in both; Q=0 iff one of the 9 rotated S3 vanishes; Groebner radical vdim = brute-force count.

## Twice-run diff
run1/results.json vs run2/results.json: byte-identical (`cmp`), covering system hashes, degrees from both engines, vdim, radical vdim, solution counts. Deterministic.

## Interpretation (separate from measurement)
The drop is dominated by replacing x^d-c with u^(d/3)-c: degree reached tracks ~d+1 (plain msolve) vs ~d/3+10 (quotient). That is a factor-base-size effect, tied to V being 3x smaller in u-space; with Bezout (B_ci) unchanged at large d it is not evidence that the S3 ideal itself became easier. Fair control not run: plain system on a factor base with the same count of free values (d/3) -- would then show whether quotient beats a plain smaller base; likely roughly neutral since the plain degree also tracks base size. Small d: quotient is worse (orbit polynomial degree 12). Singular degree here is max pair-level (bracket values uninformative, Amendment 1); parser splits concatenated levels (Amendment 2).

## Limitations / next step
Single study (T3), m=2 only, overdetermined 2-variable systems (3 equations), base a coset x^d=c only. Next: control with matched base size; m=3 chain with orbit quotient on the intermediate variable; unknown-base-size scaling at d up to ~1000.

## Files
PROTOCOL.md, RESULT.md, sym.py, cert.py, run.py, analyze.py, raw/{Q_u_poly,sym_gate,cert_S3,summary_table}.json, raw/sym_log.txt, run1/, run2/ (results.json, log.txt), fig_rank2_j0_order6_quotient.png/.svg

---
## Control section (added after coordinator request; PROTOCOL Amendment 3 written before the run)
Control = plain system on the index-3 sub-coset W = {x^(d/3)=c'} (same size d/3 as the quotient base), target re-planted in W (same p, b, x0). 58 instances, run1 == run2 (`cmp run1/control.json run2/control.json` identical); control radical vdim == brute-force count in all 58. Data: run?/control.json, analyze_ctrl.py.

Means by d (msolve / Singular degree; decompositions of R; msolve max matrix rows*cols = cost proxy):
| d | deg plain | deg quot | deg ctrl | N plain (=quot lifted) | N ctrl | maxmat plain / quot / ctrl |
|---|---|---|---|---|---|---|
| 6 | 7 / 14.2 | 12 / 14.4 | 3 / 6.6 | 2.4 | 1.4 | 100 / 2254 / 24 |
| 12 | 13 / 25.8 | 12 / 18.2 | 5 / 10.4 | 3.5 | 1.6 | 120 / 1536 / 72 |
| 36 | 37 / 74 | 18 / 33.8 | 13 / 26 | 2.8 | 2.0 | 120 / 4635 / 120 |
| 96 | 97 / 181.5 | 38 / 66 | 33 / 65 | 61.2 | 7.0 | 120 / 68780 / 120 |
| 126 | 127 / 253 | 48 / 92 | 43 / 86 | 8.0 | 2.0 | 120 / 55699 / 120 |
Totals: decompositions plain/quotient-lifted 620 vs control 132 (control finds 0.21x). Control degree < quotient in 58/58 (msolve) and 57/58 (Singular); quotient/control degree ratio averages 2.05 (msolve) / 1.5 (Singular), but at d >= 36 the quotient is only ~10-25% above the control (e.g. 48 vs 43, 92 vs 86), at d <= 12 it is 2-4x above.

Decision rule outcome: for d >= ~18-36 the control's degree is within 25% of the quotient's (both engines) while finding ~0.2x (about 1/5, between 1/3 and 1/9) as many decompositions => by the pre-registered rule the quotient is a genuine per-decomposition gain **in degree terms**: it covers the full base V at roughly the degree of the 3x-smaller base. For d <= 12 the rule gives the opposite (quotient much worse). Counts per system are tiny (1-60, mostly ~2), so the "x" ratio is noisy.
Cost-proxy caveat (contradicts a cost-gain reading): the quotient's dense degree-12 orbit polynomial makes msolve matrices 15-500x larger than plain/control (plain and control stay at ~100-144 cells since these 3-equation 2-variable systems are tiny); decompositions per matrix cell is worse for the quotient than both plain and control at nearly all d (e.g. d=96: plain 0.51, quotient 0.0009... see analyze_ctrl.py output units x1e3: 510 / 0.89 / 58). So "degree" and "matrix-size cost" disagree; degree here is largely a pair-level (lcm) effect of x^d-c, not actual work.

## Revised bottom line
The quotient's lower d_solv is mostly the base-size effect (control with the same base size is even lower), but it retains ~4.7x more decompositions at only ~10-25% higher degree for d >= 36, so per decomposition (by degree) it is a real gain over the orbit-fold-like small-base control; by the matrix-size proxy it is a loss. Net: lowers by degree per decomposition (T3, small counts), not shown to lower real cost; the algebraic lever is not demonstrated to beat the orbit fold in work terms. Figure not regenerated (table above).
