# RESULT: rank 3 - solving in invariants of the signed-permutation group

Odd p only (11..101): in characteristic 2 the Z/2 factor is modular, not covered. Toy instances, academic measurement.

## Measurement
Setup: E: y^2 = x^3 + a x^2 + c^2 x over F_{p^n}, a in F_{p^n}\F_p, c in F_p^*, Y=(x-c)/(x+c); factor base x in F_p;
n=m. Semaev S_{m+1} derived in Y (`semaev_Y.py`; S3 degree 2/var, S4 degree 4/var), Weil descent = coefficient
expansion (Y_i in F_p) -> n equations over F_p. Instances: m=2: p in {11,13,17,31,101} x 2 seeds (10); m=3:
p in {11,13,17,19,23} x 3 seeds (15). Planted R = sum of m factor-base points; Y_R in F_p instances rejected (descent degenerates).

Correctness (`cert_semaevY.py`, `raw/cert_semaevY.json`): 11 curve/field cases, 0 false positives / 0 false
negatives vs brute-force sign sums; every system (a,b,c1,c1',c2) vanishes at the planted point; planted point lies in
the variety of the GB on both engines (Singular std(G+m_planted) != 1; msolve GB evaluated at planted): all True.

Degree definitions: msolve = max first column of the F4 round table (blank = no round, then input degree).
Singular = least `degBound` d at which std already has the full vdim (the `[k:d]` protocol number is a constant 2 here
and the digit stream is ambiguous, so the harness convention was NOT usable); d_eff = max(d, input degree).
Invariant variables have degree 1 (standard grading).

Systems: a plain Y; b e_k(Y); c1 (u_k=e_k(Y^2) k<m, pi=prod Y, u_m=pi^2 substituted); c1' (u_1..u_m, pi plus
equation pi^2-u_m); c2 (B_m: product S(Y,Y_R)S(Y,-Y_R) in e_k(Y^2)).

m=n=3 (15 instances; ranges):
| system | vars | max eq deg | Bezout | vdim | Singular d_eff | msolve d |
|---|---|---|---|---|---|---|
| a Y | 3 | 12 | 1200-1728 | 384 | 27-30 | 17 |
| b e_k(Y) (S_3) | 3 | 4 | 64 | 64 | 6-8 | 8-10 |
| c1 D_3 (u,pi) | 3 | 4 | 36-64 | 16 | 6-7 | 5-6 |
| c1' D_3 + relation | 4 | 2 | 16 | 16 | 3-4 | 4-5 |
| c2 B_3 product | 3 | 4 | 64 | 64 | 10 | 11 |

m=n=2 (10 instances): a: vdim 8, Sing 5-7, msolve 4-6; b: vdim 4, Sing 2 (input degree), msolve 2; c1: vdim 2, Sing 2, msolve 2;
c1': vdim 2, Sing 2, msolve 1; c2: vdim 4, Sing 3, msolve 4.

Per-instance comparison with b (m=3): c1 < b: msolve 15/15; Singular 12/15 lower, 3/15 higher (7 vs 6).
c1' < b: 15/15 on both engines. (m=2: c1 = b on both; c1' = b on Singular, lower on msolve, trivially.)

Expansion factors (vdim ratios, exact): a->b = m! (8->4, 384->64: 2, 6); b->c1 = 2^(m-1) (2, 4); a->c1 = |D_m| = 2^(m-1) m! (4, 24).
One invariant solution <-> |D_m| ordered Y-tuples (24 at m=3); each Y also stands for the pair +-P.
c2 is NOT a faithful B_m quotient: vdim 64 = b's, not 16: the Weil-descended components of a product force
a choice of factor per Frobenius conjugate (2^n branches), producing spurious cross-conjugate solutions; also degree 10-11, worse than b.

Twice-run diff (run1 vs run2, system hashes, both-engine degrees, vdim, planted flags, 25 instances x 5 systems): identical.

## Interpretation (separate)
- Falsifier (c1 d_solv = b d_solv on both engines): NOT triggered at m=3. At m=2 it holds (tie, degenerate: 2 quadrics).
- Verdict, m=3: D_m invariants LOWER d_solv: msolve 15/15 (b 8-10 -> c1 5-6 -> c1' 4-5); Singular is engine-dependent for the
  minimal variant c1 (12 lower / 3 higher, 6-8 -> 6-7) but clear for c1' (6-8 -> 3-4, 15/15). Lowering is not only Bezout/vdim (64 -> 16).
- Mechanism: in Z=Y^2 the generators have degree 2 per variable, so c1' is a complete intersection of 4 quadrics (Bezout 16 = vdim);
  d_solv tracks the Macaulay-type bound (b: three quartics, bound 10 = msolve max; c1': 4 quadrics, bound 5 = msolve max).
  So the gain is a generator-degree effect from the presentation, matching the Bezout count; c1' needs an extra variable.
- B_m is unavailable for fixed R (odd flips send R to R+T); the product trick loses.
Tier: T3 (single study, raw data, twice-run identical, no independent replication). m<=3, p<=101 only; extrapolation to m=4 (T4/T5) untested.
Limitations: grading choice (variables unit degree; natural weights 2k, m would change the comparison); Singular degBound metric is mine;
n=m only; odd p only.
Next step: m=4 over F_{p^4} (D_4, Bezout 4096->512) with c1', both engines, to see if d_solv follows the quadric bound.

Files: PROTOCOL.md, semaev_Y.py, lib.py, systems.py, run_solvers.py, driver.py, cert_semaevY.py, analyze.py,
raw/{run1,run2}.json, raw/cert_semaevY.json, raw/summary_table.txt, fig_rank3_reflection_group_invariants.{png,svg}.
