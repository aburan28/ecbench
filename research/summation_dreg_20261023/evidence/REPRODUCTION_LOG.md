# Reproduction log (run once, then twice)

Host: 4 vCPU, 15 GB, Linux 6.18.44. Tools: msolve 0.6.5, Singular 4.3.3.0,
ecbench (pinned crypto `0dc9ef4`). crypto repo at `af663ad7c`.

## A. Independent harness — my own systems, deterministic invariants

Canonical grid (20 systems: chain m=2/m=3, unsplit S4, split) run **twice**:
- `canonical_run1.json` vs `canonical_run2.json`: **FULL FILES BIT-IDENTICAL** (0 mismatches).
- Reproducible invariants: system sha256, first-fall degree, Singular solving
  degree, vdim, dim, Macaulay rank profile. (Timings excluded by construction.)

Correctness certificates (independent of repo):
- Prime S3: 0 false pos / 0 false neg over 5 non-singular curves (p=13..29).
- Binary S3: 0/0 over 6 curves incl. Koblitz a2∈{0,1},a6=1 (n=3,4,5).
- Binary S4 (resultant of two S3 quadratics): 0/0 over 4 curves (n=4,5,6).

## B. Repo's PKM-tower msolve degrees — reproduced TWICE, exact

`msolve -t 1 -v 2 -f research/ic_next10_20261014/03-pkm-degree/systems/<name>.ms`

| system | recorded D / rounds / rows×cols | run1 | run2 |
|---|---|---|---|
| m2_t6_x464833983 (N=12) | 5 / 14 / 2803×2064 | 5/14/2803×2064 | 5/14/2803×2064 |
| m3_t3_x388393539 (N=9)  | 6 / 17 / 3827×2997 | 6/17/3827×2997 | 6/17/3827×2997 |
| m4_t2_x92234565 (N=8)   | 6 / 22 / 2009×1906 | 6/22/2009×1906 | 6/22/2009×1906 |
| m3_t4_x334226296 (N=12) | 6 / 19 / 12900×12287 | 6/19/12900×12287 | 6/19/12900×12287 |

run1 = run2 = recorded, every field, every system. (m4_t3_x259303215, recorded
D=7 at 594 s, not re-run twice for budget; single run matches dims 26025×21665.)

## C. ecbench replays — reproducible vs non-reproducible, each run TWICE

- PASS: `ecbench verify pkm_tower_ecbench_20261013/sessions/p20 --replay-all`
  run1 & run2 identical: ok=true, 192 records, 192 verified, 144 replays, 0 problems.
- FAIL: `ecbench verify ecbench_rootfind_engine_20261012/sessions/claim-n13 --replay-all`
  run1 & run2 identical: ok=false, 224 records, 224 verified, 168 replays, **96 problems**.
  All 96 failures are IC (F4) arms; all 72 rho rows match. The failure is
  deterministic on this host but the mismatch is cross-code-path: the recorded
  Apple-arm64 word-XOR counts differ from this x86 host's BMI2 pext/pdep counts
  (`inherited_f4.rs:942`, bug open on main). So the recorded claim does not
  replay across code paths — a genuine, reproduced reproducibility defect.

  **Correction (v2, 2026-10-23).** The measurement above stands (96 problems, twice);
  the attribution does not. On crypto `main` (a53a7e2) the BMI2 and portable rewrite paths
  charge identical work — #2043 (2026-10-17) made every host run the eager rewrite — as
  shown by a forced-path regression test and by re-running the claim-n13 spec on each
  path (224/224 records identical). The failure comes from a counted-work definition
  change: the session was recorded at `2245b75`, when a non-BMI2 host charged specialisation
  lazily; a temporary lazy build reproduces the recorded `ic-f4-m2` counts on 32/32 rows.
  See aburan28/crypto#2459 and REPORT §4.8.

## D. Exotic levers — toy falsification, each config many subspaces

- Rank-1 (factor-base subspace V): 295 subspaces across 5 configs (unsplit/split/
  m=2), solving degree min=median=max at every config; Frobenius-stable V identical.
  Lever **falsified**: d_solv is intrinsic to the ideal class (V is not a knob).
- Rank-8 (monomial/variable order): grevlex=grlex and all 4 random variable
  permutations give identical degree (3 unsplit / 5 split). Only pure lex differs
  (lex "degree" is not the solving degree — different algorithm regime).

## E. Engine-dependence of the measured "degree" (corroborates survey A)

Same ideal, two external engines (msolve F4 with fall reuse vs Singular std):

| system | nvars | Singular std | msolve F4 |
|---|---|---|---|
| split chain n8 ℓ2 m3 | 14 | 5 | 4 |
| split chain n10 ℓ2 m3 | 16 | 5 | 3 |
| split chain n8 ℓ3 m3 | 17 | 6 | 4 |
| unsplit S4 n8 ℓ2 m3 | 6 | 3 | 6 |
| unsplit S4 n10 ℓ2 m3 | 6 | 3 | 6 |
| unsplit S4 n9 ℓ3 m3 | 9 | 4 | 7 |

On the SPLIT chain, msolve-F4 reaches 3–4 where Singular reaches 5–6 — reproducing
survey A's "fall reuse keeps refutation degree 3–4 where XL needs 5–7". On the
UNSPLIT S4 the order reverses. Conclusion: the "degree" is engine- AND
encoding-dependent; only the coordinate-free last-fall / CM-regularity invariant
(Caminata–Gorla) is well defined — and that is exactly what resists cheap measurement.
