# PROTOCOL: rank 3 - solving in invariants of the signed-permutation group (written before running)

## Hypothesis
For E: y^2 = x^3 + a x^2 + c^2 x over F_{p^n} (rational 2-torsion T=(0,0), c in F_p^*, a in F_{p^n} generic,
p odd), the coordinate Y=(x-c)/(x+c) satisfies Y(P+T)=-Y(P), Y(-P)=Y(P). The Gaudry-Diem decomposition
polynomial S_{m+1}(Y_1..Y_m, Y_R) (R fixed) is invariant under D_m (even sign changes) x S_m. Rewriting the
Weil-descended system in D_m-invariants may lower the Groebner-engine solving degree d_solv relative to
S_m-invariants, not only the Bezout number. (Odd sign changes map R -> R+T, i.e. Y_R -> -Y_R, so for fixed R the
honest group is D_m; the full B_m acts only on the product S(Y,Y_R)*S(Y,-Y_R).)
Odd p only: in char 2 the Z/2 factor is modular (char | |Z/2|), invariant theory differs; this test does not cover it.

## Setting
n = m in {2,3}; p odd small; E over F_{p^n} with a in F_{p^n} \ F_p random (so factor-base points are not a
subgroup), b=c^2, c in F_p^*. Factor base: points with x (equivalently Y) in F_p. Planted instance: choose
P_1..P_m in E(F_{p^n}) with x_i in F_p, R = sum P_i. Descent: Y_i in F_p, so expand the F_{p^n}-coefficients of S
in the basis 1,t,..,t^{n-1} -> n equations over F_p in m unknowns.

## Systems compared (same ideal, different generators/ring)
(a) plain: unknowns Y_1..Y_m.
(b) S_m-invariants: e_1..e_m of Y.
(c1) D_m invariants: u_k=e_k(Y^2) (k<m) and pi=Y_1..Y_m, with u_m=pi^2 substituted (m variables).
(c1') same with u_m kept and equation pi^2-u_m added (m+1 variables).
(c2) B_m invariants of the product S(Y,Y_R)S(Y,-Y_R): unknowns u_k=e_k(Y^2), k=1..m (solves for R or R+T).
Measured: msolve F4 (-t 1 -v 2, degree = max first column of round table) and Singular std (dp, option(prot),
max [k:d]); vdim (Singular), number of variables, equation degrees (Bezout = product of equation total degrees),
planted decomposition recovered (planted point lies in V(GB): std(G + maximal ideal at planted) != <1>) in each.
Correctness first: S_3 and S_4 in Y checked against brute-force sign-sum over E(F_p) (0 false pos / 0 false neg),
and S(Y_planted; Y_R)=0 over F_{p^n}.
Each deterministic measurement run twice (run1/run2), diff of system hash, degrees, vdim.

## Falsifier (fixed)
FALSIFIED (neutral) if (c1) d_solv == (b) d_solv on BOTH engines in every instance (the extra Z/2 buys nothing
in degree, only in Bezout/vdim). "Lowers" requires (c1) < (b) on both engines in all/most instances; if only one
engine drops, report "engine-dependent". "Raises" if (c1) > (b). Also reported: per-solution expansion factor
(|S_m|=m!, |D_m|=2^{m-1} m!, |B_m|=2^m m! ordered Y-tuples per invariant solution, generic).
Tier ceiling: T3 (single measured study, raw data) unless certificates are independently replicated.
