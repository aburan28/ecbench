# PROTOCOL - rank 4: Artin-Schreier (AS) trace elimination before Groebner (char 2)

Written before any measurement. Academic toy-scale measurement only.

## Hypothesis
The low first-fall degree of Weil-descended binary S_3 systems comes from the trace /
Artin-Schreier structure (x on E iff Tr(x + a2 + a6/x^2) = 0; S_3 is F2-AFFINE in its third
argument: S3(x1,x2,A) = alpha*A^2 + beta*A + gamma, solvable iff Tr(alpha*gamma/beta^2)=0).
Splitting every unknown field element into the part the absolute trace sees (one F2-linear
form Tr(.)) and the trace-zero part, and eliminating the trace block up front, may lower the
degree a Groebner engine needs on the residual (H1), or only remove the easy part (H0, null).

## Design of the reduction (exact)
Curve y^2+xy=x^3+x^2+1 (a2=a6=1) over GF(2^n), n in {8..12}; V = F2-subspace of dim l in
{2,3} (seed-1 subspace of the harness), m in {2,3}.
* Trace-adapted coordinates (invertible constant F2-linear change of unknowns, ideal unchanged):
  - factor-base x_i in V: basis b_1..b_l of the same V with Tr(b_1)=1, Tr(b_j)=0 (j>=2), so
    Tr(x_i) = c_{i,1} (the "trace bit"). If Tr|_V == 0 the block is empty (reduction vacuous; next seed is used).
  - free chain aux A (split model, n bits): basis e_0 (Tr=1), e_1..e_{n-1} (Tr=0); Tr(A)=a_0.
* RESIDUAL (AS arm): fix every trace bit t in F2^k (k = #unknown field elements: m for the
  unsplit model, m+(m-2) = 2m-2... i.e. x_1..x_m and A_1 for the split m=3 model), substitute,
  and solve each of the 2^k branches over the remaining trace-zero coordinates (l-1 bits per x_i,
  n-1 bits per aux). The trace block is eliminated in closed form (the enumerated linear form);
  the solution set is the union over branches; lift back x_i = t_i b_1 + sum r_j b_j, A = a_0 e_0 + sum a_k e_k.
  Residual degree = MAX over branches (worst case the solver must face) and also the planted branch.
* CONTROL arm (same cost structure, no trace): fix a non-trace linear form per unknown element
  (last V-coordinate c_{i,l}, aux bit a_{n-1}). If CONTROL == AS everywhere, any change is a plain
  variable-count (hybrid) effect and not AS-specific.
* Why this differs from the program's "trace-zero pre-filter" (0/29,624 hits on odd-order IC targets):
  that is a FILTER on candidate points/relations (reject those whose trace is wrong), leaving the
  polynomial system untouched, and it never fires. Here the system itself is transformed
  (substitution, fewer unknowns, trace-block removed from the aux too) and it does not depend on
  hitting any target property; it is a change of the ideal handed to the engine.

## Systems compared
FULL (trace-adapted coordinates; also FULL in original harness coordinates as reference) vs
RESIDUAL-AS vs CONTROL, for models: m=2 (single S3 with xR const), m=3 unsplit (`S4_unsplit`),
m=3 split (S3 chain with one free aux). n in {8,9,10,11,12} x l in {2,3}.

## Correctness checks (before degree)
Per cell: (1) brute force over V^m (and aux in F_q) of the numeric S3/S4 system gives solution set
SF; brute-force point decompositions BF (true point sums with signs) satisfy BF subset SF;
(2) union over AS branches of branch solutions == SF exactly and sum of branch vdim == FULL vdim ==
numeric solution count (0 lost / 0 invented); (3) planted decomposition satisfies the planted
branch system and that branch has vdim>=1; (4) the numeric count equals Singular vdim.

## Measurements
d_solv from msolve F4 (`-t 1 -v 2`, max degree column) AND Singular std (max [k:d]); vdim;
nvars, neqs, max input degree; first-fall degree via `pdp.macaulay_rank_seq` where affordable.
All solver calls: `nice -n 10 timeout 300 prlimit --as=2500000000`. Timeout = recorded ">= not reached".
Run twice (run1/run2); diff system hash, degrees, vdim.

## Falsifier (fixed; from survey_E rank 4, sharpened)
For a cell and engine: RESIDUAL_AS max-branch d_solv == FULL d_solv. Lever is FALSIFIED on an engine if this
holds in every cell (or all but where FULL is already at the input-degree floor); LOWERS if AS
< FULL on both engines in most cells AND AS < CONTROL somewhere (AS-specific); NEUTRAL-by-hybrid if
AS < FULL but AS == CONTROL (generic variable-fixing artefact, not Artin-Schreier); RAISES if AS > FULL.
Interpretation of outcomes:
* equal: trace block was the easy part; hard core unchanged (H0).
* lower & > control: trace block carried degree (H1).
* lower & == control: no AS content, only fewer variables.
No edits to this falsifier after results; amendments will be dated below.

## Amendments
(none)
