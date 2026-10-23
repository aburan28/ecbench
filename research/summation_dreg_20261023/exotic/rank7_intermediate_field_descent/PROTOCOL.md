# PROTOCOL — rank 7: descend to an intermediate field GF(2^k) vs. all the way to F_2
Written 2026-10-23 BEFORE any run.

## Hypothesis (harness check against Diem / Gaudry)
For n = k*m, E over GF(2^n)=GF(q^m), q=2^k, factor base F = {P : x(P) in GF(q)}.
The decomposition equation S_{m+1}(x_1..x_m, x_R)=0 Weil-descends over GF(q) to m
equations in m unknowns whose total degrees are bounded by deg S_{m+1} (independent of k),
so the Bezout number of the descended system is bounded independent of k.
Descending the same equation to F_2 (k*m boolean unknowns, field eqs) gives a system whose
Groebner solving degree should grow with k*m.
For prime n there is no proper subfield, so the GF(q)-descent does not exist.

## Systems compared (same curve, same planted decomposition, same R)
A. GF(q)-descent: unknowns x_i in GF(2^k), Singular `ring r=(2^k,a),(x1..xm,t),dp`; GF(2^n) built as
   GF(q)[t]/(g), g = min poly of the harness generator z over GF(q) (computed with harness gf2 arithmetic,
   so the tower is the harness field). Equations = the m t-coordinates of S_{m+1} (reduced mod g).
   A0: no field equations (generic, Diem regime; degrees/Bezout independent of q).
   A1: plus x_i^q + x_i (rational-point count).
B. F_2-descent: x_i = sum_j c_ij * alpha^j (alpha in GF(2^n) a root of IRRED[k]; basis of the embedded
   subfield), boolean_weil.var_elt, weil_equations, + c^2+c. Degree from msolve F4 AND Singular std.
msolve cannot do GF(2^k), k>1 (prime fields only), so A is Singular-only; B has both engines.
S_3 for m=2, S_4 (resultant form, harness S4_unsplit) for m=3.
Configs (n,k,m): m=2: (4,2),(6,3),(8,4); m=3: (6,2),(9,3),(12,4). Curves a2=1,a6=1 (Koblitz) and a2=1,a6=generic
element (non-Koblitz) for the small configs.

## Correctness gate
For each config: brute-force the set of ordered x-tuples in V^m with some point choice (+-) summing to +-R via
the group law (bin_add); require: #solutions of A1 (vdim) == #solutions of B (vdim) == brute-force count,
and the planted tuple is a root of A0, A1, B. (0 false pos / 0 false neg in count.)

## Falsifier (sharpened from survey_E §5 rank 7)
Lever "lowers degree" is FALSIFIED (harness not reproducing Diem) if, in A0 with m fixed, the Singular
solving degree or the descended-equation total degrees grow with k, OR if the descended equation degrees exceed
deg S_{m+1}. Secondary (expectation): B's degree grows with k*m while A0's stays <= a function of m only.
Outcomes: A0 bounded & B grows => harness reproduces the theorem (consistent; says nothing new about prime n).
Otherwise: harness/encoding issue or theorem-regime misunderstood => inconclusive.
Prime n in {11,13}: show no proper subfield (divisors) and that B with an l-dim subspace is the only option.

## Determinism
Each measurement run twice (run1/run2); diff system hash, degrees, vdim. Timings not compared.
Resources: nice -n 10, timeout 300, prlimit --as=2.5e9, msolve -t 1.
