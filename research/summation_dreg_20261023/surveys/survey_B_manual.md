# Survey B (symmetry / Koblitz) — reconstructed manually

Surveyor B hit a safety classifier and produced no file. These facts were
extracted directly from the repo by the lead session (file + line cited), kept
to measured numbers and their honest status.

## Symmetrisation lowers per-variable DEGREE (the genuine lever)
- Torsion (2-/4-) symmetrisation rewrites S_m in FGHR/FHJRV invariants; for S_4 it
  lowers per-variable degrees `[4,4,4,4] -> [2,2,2,2,1]`.
  Source: `research/koblitz_symmetrised_e2e_20260927/README.md:21`. Tier T2/T3
  (reproduced; FGHR/FHJRV EUROCRYPT 2014).
- Dryło 2026/1299 certifies "FGHR symmetrisation = Semaev polynomial of E/G"
  (quotient curve), i.e. the gain is mostly a change of curve, not free. (survey_E §4.)

## D4-weighted S5 ("C8" / lead 5b) — the strongest symmetry study, measured honestly
Source: `research/d4_symmetrised_s5_20261015/README.md`, PR #1785.
- Group `D_4 = (Z/2)^3 ⋊ S_4`, order 192, on four summands' FGHR coord `Y_i`.
- Bézout `4096 -> 512`: **derived** two ways (weighted Bézout `8^4/(2·2·2·1)`;
  orbit count `|D4|/|S4|=8`). Labelled "derived" in-repo.
- **Measured**: quotient dimension = 512 on every solved target at `p>=29`
  (3/5 at p=17). Weighted Macaulay bound 26; matrix **4,200 columns** vs plain
  form's **40,920** (and C8's 11,628). Plain form Macaulay degree 29.
- **Measured speedup: 48× on the C_4 relation cost** (`C_4 = 2.54e10` at p=269),
  NOT the ≈630× projected. Explicitly class "engineering (a constant of one
  phase). No end-to-end S measured, and none claimed."
- CORRECTION to survey_E: survey_E called this "the repo's strongest measured
  dreg lever (Macaulay 29→14)"; the repo files show the *solved* object is the
  weighted system (bound 26, 4200 cols) with a **measured 48×**, and 29 is the
  plain-form Macaulay degree. The 630× is projected, not run. Tier: column/Bézout
  reduction T3 (measured quotient dim), 48× speedup T3, 630× T4 (projected).

## G2 — C_n-equivariant F4 on Koblitz (abelian Frobenius action)
Source: `RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md:186,191`.
- Character-block (Faugère–Svartz) decomposition; verified EXACT at n=7,11.
- Rank-only speedup up to **128×** at (n=23, ℓ=11): matrix 120379×273262 →
  block 5297×11880. Sequence 39×,105×,105×,128× matches the cost model.
- This is PER-DEGREE LINEAR ALGEBRA, not a d_reg reduction. Tier T3 (rank-only;
  full-solve speedup left open — "possibly mislabeled closed", survey_E §4).

## Koblitz end-to-end: symmetry helps refutation, can lose overall
Source: `research/koblitz_symmetrised_e2e_20260927/README.md:23-24`.
- Symmetrised system refutes non-decomposable targets **×300–6000 faster at n=15**
  but **LOSES ×1.6 at n=17, m=3** end to end. Honest mixed result. Tier T3.

## Structural immunity for prime n (the Koblitz target curves)
Source: `RESEARCH_KOBLITZ_INDEX_CALCULUS.md` (via survey_E §4, fact 3).
- Frobenius-stable F_2-subspaces are cyclic codes; their dimensions are subset-sums
  of 2-cyclotomic coset sizes. At n=131 the only sizes are {0,1,130,131} → empty at
  useful dims. Same for n=233,283 (prime, 2 of large order). Couveignes–Lercier
  elliptic-period bases also empty for n≥7. So the Frobenius-subspace and
  intermediate-field levers are STRUCTURALLY EMPTY for ECC2K-130/sect233k1/sect283k1.
  Tier T1-structural (a proof/cyclotomic-coset argument).

## Net reading for the report
Symmetry/automorphism quotients are the only proven large lever, and they act by
(a) lowering per-variable degree (4→2) and (b) cutting Bézout / column count by the
group order — both real and measured. But: they are per-phase constants that do not
bend the `d_reg(ℓ)` growth; the measured end-to-end effect is sometimes a loss
(×1.6 at n=17); and for prime n they are structurally empty. None crosses rho for
the target curves.
