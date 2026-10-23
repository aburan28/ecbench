# RESULT - rank 4: Artin-Schreier trace elimination (char 2)

## Verdict (measurement)
* Singular std: residual d_solv < FULL in 22/27 cells (the 5 ties are m=2, already at the floor 2), but residual == CONTROL in 27/27. Lowers vs FULL, with zero Artin-Schreier-specific content: NEUTRAL (generic variable-fixing / hybrid artefact).
* msolve F4: residual < FULL in 19/27, == FULL in 8/27; residual == CONTROL in 25/27 (AS<CTL once, AS>CTL once; +-1 noise). Same reading: NEUTRAL-by-hybrid.
* Per model: unsplit (S4_unsplit): Singular 3->2, 4->3; msolve 6->3 (l=2), 7->6 / 6->6 (l=3). Split chain: Singular 5->4 (l=2), 6->5 (l=3); msolve 3-4 -> 3-4 (mostly tie or -1). m=2 single S3: Singular 3->2 (l=3) else tie.
* Literal survey_E falsifier (residual d_solv = full d_solv) is NOT met (residual is lower), but the sharpened one (AS < CONTROL somewhere) is: the lower degree is exactly what fixing ANY one linear form per unknown gives. The trace form is not special. The degree-driving part is not the trace block.
* The reduction removes n bits of the aux only to a fixed trace bit (n -> n-1) and one bit per factor-base element; the hard core (degree 5-6 Singular on the split chain, 6 on msolve unsplit at the input-degree floor) remains.

## Reduction (see PROTOCOL.md for the exact definition)
Constant F2-linear change of unknowns to a trace-adapted basis (Tr(b_1)=1, others 0; aux e_0 trace 1), then enumerate (fix) all trace bits, solve each of 2^k branches on the trace-zero coordinates; lift x_i = t_i b_1 + sum r_j b_j. Residual degree = max over branches. CONTROL fixes a non-trace coordinate per unknown instead.
Note: pure change of basis cannot change the ideal; only substitution (branching) can, hence the reduction is effectively a hybrid guess of one F2-form per unknown, and the control isolates what is AS-specific.

## Correctness (before degree) - all 27 cells x 2 runs
* Numeric solution set SF of each model is enumerated by brute force; true point decompositions BF: planted decomposition in SF in every cell; BF subset SF for m=2 and unsplit; for the split chain every NON-degenerate BF (P1+-P2 != O) is in SF (BF misses: only degenerate partial-sum-O tuples, which the S3 chain cannot express - known property, not caused by the reduction). Unsplit/split also contain spurious twist/non-F_q tuples (counts in raw JSON).
* FULL Singular vdim == numeric solution count; for every branch of both arms the branch vdim == the numeric count of solutions with that trace-bit pattern; sum of branch vdims == FULL vdim (0 lost, 0 invented): `ok=True` in all 27 cells. Planted assignment satisfies the system (asserted at build) and lies in the planted branch.
* Singular/msolve timeouts: none occurred in the 27 cells.

## Table (run1; seed-1 harness subspace, bumped only if Tr|_V = 0). ff = Macaulay first-fall for FULL ('-' = none up to D)
| n | l | m | model | vars FULL | vars/branch | #branches | max in-deg | vdim | Sing FULL | Sing AS | Sing CTL | msolve FULL | msolve AS | msolve CTL | FULL ff (<=D) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 2 | 2 | single | 4 | 2 | 4 | 2 | 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2(<= 5) |
| 8 | 2 | 3 | split | 14 | 10 | 16 | 3 | 7 | 5 | 4 | 4 | 4 | 3 | 3 | 2(<= 5) |
| 8 | 2 | 3 | unsplit | 6 | 3 | 8 | 6 | 10 | 3 | 2 | 2 | 6 | 3 | 3 | 5(<= 5) |
| 8 | 3 | 2 | single | 6 | 4 | 4 | 2 | 2 | 3 | 2 | 2 | 3 | 2 | 2 | 2(<= 5) |
| 8 | 3 | 3 | split | 17 | 13 | 16 | 3 | 6 | 6 | 5 | 5 | 4 | 4 | 4 | 2(<= 4) |
| 8 | 3 | 3 | unsplit | 9 | 6 | 8 | 6 | 6 | 4 | 3 | 3 | 7 | 6 | 6 | -(<= 5) |
| 9 | 2 | 2 | single | 4 | 2 | 4 | 2 | 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2(<= 5) |
| 9 | 2 | 3 | split | 15 | 11 | 16 | 3 | 3 | 5 | 4 | 4 | 3 | 3 | 3 | 2(<= 4) |
| 9 | 2 | 3 | unsplit | 6 | 3 | 8 | 6 | 10 | 3 | 2 | 2 | 6 | 3 | 3 | 4(<= 5) |
| 9 | 3 | 2 | single | 6 | 4 | 4 | 2 | 2 | 3 | 2 | 2 | 3 | 2 | 3 | 2(<= 5) |
| 9 | 3 | 3 | split | 18 | 14 | 16 | 3 | 7 | 6 | 5 | 5 | 4 | 3 | 3 | 2(<= 4) |
| 9 | 3 | 3 | unsplit | 9 | 6 | 8 | 6 | 22 | 4 | 3 | 3 | 7 | 6 | 6 | -(<= 5) |
| 10 | 2 | 2 | single | 4 | 2 | 4 | 2 | 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2(<= 5) |
| 10 | 2 | 3 | split | 16 | 12 | 16 | 3 | 7 | 5 | 4 | 4 | 3 | 3 | 3 | 2(<= 4) |
| 10 | 2 | 3 | unsplit | 6 | 3 | 8 | 6 | 10 | 3 | 2 | 2 | 6 | 3 | 3 | 4(<= 5) |
| 10 | 3 | 2 | single | 6 | 4 | 4 | 2 | 1 | 3 | 2 | 2 | 2 | 2 | 2 | 2(<= 5) |
| 10 | 3 | 3 | split | 19 | 15 | 16 | 3 | 3 | 6 | 5 | 5 | 4 | 4 | 3 | 2(<= 3) |
| 10 | 3 | 3 | unsplit | 9 | 6 | 8 | 6 | 3 | 4 | 3 | 3 | 7 | 6 | 6 | -(<= 5) |
| 11 | 2 | 2 | single | 4 | 2 | 4 | 2 | 2 | 2 | 2 | 2 | 3 | 2 | 2 | 2(<= 5) |
| 11 | 2 | 3 | split | 17 | 13 | 16 | 3 | 3 | 6 | 5 | 5 | 4 | 3 | 3 | 2(<= 4) |
| 11 | 2 | 3 | unsplit | 6 | 3 | 8 | 6 | 3 | 3 | 2 | 2 | 6 | 3 | 3 | 5(<= 5) |
| 11 | 3 | 2 | single | 6 | 4 | 4 | 2 | 2 | 3 | 2 | 2 | 3 | 3 | 3 | 2(<= 5) |
| 11 | 3 | 3 | split | 20 | 16 | 16 | 3 | 7 | 6 | 5 | 5 | 4 | 3 | 3 | 2(<= 3) |
| 11 | 3 | 3 | unsplit | 9 | 6 | 8 | 6 | 22 | 4 | 3 | 3 | 6 | 6 | 6 | -(<= 5) |
| 12 | 2 | 2 | single | 4 | 2 | 4 | 2 | 1 | 2 | 2 | 2 | 2 | 2 | 2 | 2(<= 5) |
| 12 | 2 | 3 | split | 18 | 14 | 16 | 3 | 1 | 6 | 5 | 5 | 4 | 3 | 3 | 2(<= 4) |
| 12 | 2 | 3 | unsplit | 6 | 3 | 8 | 6 | 1 | 3 | 2 | 2 | 6 | 3 | 3 | 5(<= 5) |

Full columns (original-coordinate FULL, planted-branch degrees, planted-branch first-fall) are in raw/summary.json. First-fall of FULL: 2 for single S3 and for the split chain (flat, the Kosters-Yeo feature), 4-5 or >D for unsplit; the AS planted branch first-fall is mixed (2-3 on split, 2-5 or none up to D on unsplit; see summary.json), no systematic first-fall drop; AS vs CONTROL first-fall was not compared.

## Twice-run diff
run1 vs run2: 27/27 cells identical on the deterministic content (system SHA, degrees both engines, vdim, branch data, first-fall, checks): diff empty (aggregate.py; timings excluded).

## Tier
T3 (single measured study with raw data, run twice, correctness cross-checked against brute force). Not T1/T2.

## Limitations
n=12, l=3 (3 cells) not run (90-min budget on a loaded host); one seed per cell; per-branch max is the metric (a user would also pay 2^k branches); AS arm is a hybrid, not a closed-form polynomial elimination: the true AS solution A = (beta/alpha)*HalfTrace(alpha*gamma/beta^2) is rational in the factor-base coordinates, so a polynomial closed-form elimination of the aux trace/solution would blow up degree (Diem-type), and Tr(1/x) validity constraints are non-linear, so they cannot be eliminated linearly - this is the structural reason the lever is hollow. Engines disagree strongly on FULL unsplit (Singular 3-4 vs msolve 6-7), as in REPORT section on engine dependence.

## Next step
None recommended for the AS lever. If pursued: add the validity constraint Tr(x+a2+a6/x^2)=0 as extra equations (extra consequence generators) and test only that against FULL, with the same control.

## Files
PROTOCOL.md, as_lib.py, run_cells.py, aggregate.py, fig.py, raw/run1, raw/run2 (cell JSON), raw/*.log, raw/summary.json, raw/aggregate_out.txt, fig_rank4_artin_schreier.png/.svg.
