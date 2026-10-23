# PROTOCOL (written 2026-10-23 before any run) -- rank 2: order-6 automorphism quotient, j=0

Hypothesis: for E: y^2=x^3+b over F_p (p=1 mod 3, w a primitive cube root of 1), S_3 (a=0) is
semi-invariant, S3(wx1,wx2,wx3)=w*S3. Hence the orbit product
Q(x1,x2;xR)=prod_{k1,k2 in Z/3} S3(w^k1 x1, w^k2 x2, xR) is a polynomial in u1=x1^3,u2=x2^3.
Factor base V={x : x^d=c}, 3|d, d|p-1, c = x0^d for a point x0 (so V closed under x->wx,
and u=x^3 ranges over all d/3 roots of u^(d/3)=c, each with 3 preimages in V).

Systems (m=2, target x_R of R=P1+P2 with P1,P2 planted in V; both systems always solvable):
- PLAIN:    x1^d-c, x2^d-c, S3(x1,x2,xR)           (vars 2)
- QUOTIENT: u1^(d/3)-c, u2^(d/3)-c, Q(u1,u2;xR)    (vars 2)
- target-orbit variant (xR -> w^k xR): by semi-invariance S3(w^a x1,w^b x2,w^k xR)=w^k S3(w^(a-k)x1,w^(b-k)x2,xR),
  so the 27-factor product equals Q^3 up to unit; it adds nothing (expected; verified symbolically).
m=3 (S_3 chain) is NOT run (time; recorded as T5).

Correctness gates (before degrees): (i) semi-invariance and u-polynomiality symbolic;
(ii) S3 vs real group law on small curves (0 FP/0 FN on rational-y points);
(iii) for each instance: solution set of each system from Groebner (Singular, vdim + brute-force
enumerated V^2) equals brute force; quotient solutions expanded by w-cube-roots (9 lifts) and
filtered by S3 must equal the plain solution set exactly; planted (x1,x2) recovered.

Measurements per (p,d,b,R): variables, max equation degree, Bezout.
 B_all = product of the 3 equation degrees; B_ci = min over pairs of equations of degree product
 (complete-intersection bound on a 2-var system; this is the Bezout used by the falsifier).
 vdim (Singular), solutions count; solving degree = max degree in msolve F4 round table (-t 1 -v 2 -g 2)
 and max [k:d] in Singular std protocol (dp ordering, option(prot)). Timeout -> recorded "not reached".
Costs per decomposition: report vdim and N_plain/N_quot (solutions) so cost per plain decomposition found
= (solve quotient + 9-lift check) vs solve plain.
Every measurement run twice (run1/run2); diff on system hash, degrees, vdim, solution counts.

FALSIFIER (fixed, per survey_E rank 2, sharpened): over the grid of instances, the lever is FALSIFIED
(verdict neutral/raises) if quotient d_solv >= plain d_solv on BOTH engines in a majority (>=2/3)
of instances AND the B_ci ratio plain/quotient is < 2x in a majority of instances.
"Lowers" requires quotient d_solv < plain d_solv on both engines in a majority of instances.
Mixed engine outcomes -> report per engine; if the engines disagree on direction, verdict inconclusive-by-engine.

## Amendment 1 (2026-10-23, before the grid run; falsifier text above unchanged)
Singular `[k:d]` bracket values are uninformative in this build (a single constant `[4294967295:2]` per run), so
the Singular degree is the maximum pair-level number printed in the std protocol (the bare integers, e.g.
`4(2)s12ss14...`); the bracket maximum is recorded separately as `bracket_deg`. The same Singular rule is used for
both systems. Orbit-variant check done numerically (27-factor product = unit * Q^3) plus symbolic for Q.
Amendment 2 (same day, before final run): Singular prints consecutive pair levels with no separator (e.g. `147148`);
the parser splits such digit runs into two non-decreasing levels (<=8 apart). Parser is identical for both systems.

## Amendment 3 (2026-10-23, written BEFORE running the control; coordinator request)
CONTROL: plain system x1^(d/3)-c', x2^(d/3)-c', S3(x1,x2,xR') on the index-3 sub-coset W={x: x^(d/3)=c'}, c'=x0^(d/3)
(x0 = the planted x1 of the same (p,d,b,seed) instance; W is a subset of V, same size d/3 as the quotient's u-base).
Target xR' re-planted as x0 + P2, P2 with x in W (same p, curve, c; target differs from the original because the
original x2 is in W only w.p. 1/3). Both original systems re-measured for a matrix-size proxy.
Reported per instance: d_solv (msolve, Singular), vdim, distinct decompositions (plain: N on V^2; quotient: lifted count,
equal to plain's; control: N on W^2), decompositions per system, and decompositions per cost proxy
(proxy = max msolve F4 matrix rows*cols over rounds; labelled proxy).
DECISION RULE: if control d_solv ~ quotient d_solv (within 25%, both engines) while control finds ~1/3-1/9 as many
decompositions (mean N_ctrl <= 0.5 * mean N_plain) => quotient is a genuine per-decomposition gain (searches full V at the
reduced-base degree). If control matches on degree AND decompositions (N_ctrl >= 0.5 N_plain) => only the orbit fold /
base-size effect. Otherwise inconclusive. Run twice.
