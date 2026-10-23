# Survey D: the evidence and confidence map, and the prime-field picture

- **Date:** 2026-10-22/23 (UTC).
- **Repository:** `/home/user/aburan28/crypto` at `af663ad7c` (merge of #2382). It was read only: `git status` stayed clean through every run.
- **Machine-readable twin:** `survey_D_ledger.jsonl`, in the same folder. It has 77 rows, one per technique, row or replay.
- **Host for the replays in this survey:** Intel Xeon @ 2.10 GHz with 4 vCPU, 15 GB of RAM and BMI2, running Linux 6.18.44.
  - msolve 0.6.5-1build2.
  - `ecbench` at `/home/user/ecbench/target/release/ecbench`, sha256 `6ad57b87…`. It was built from the ecbench repo pinned to crypto `0dc9ef4` (#1845, 2026-10-14). Its outputs went to the scratchpad and never into the repo.

## 0. Tier rules as applied

| tier | rule used here |
|:--|:--|
| T1 | A gate-complete `VALIDATION.json`, one that `validation_gate.py evidence` reports as valid, or an equivalent machine-checkable certificate that covers the headline number. |
| T2 | The headline counts or scalars were re-executed and matched, either by a different host, binary build, engine or independent checker, or by the validation sweep's seeded replays. |
| T3 | A single measured study with raw data in the repo. This includes same-host, same-binary self-replays, failed replays and stopped runs. |
| T4 | Prose, a derivation or an extrapolation. |
| T5 | Never run. |

**Result.** The index-calculus side has **no T1 row**. The tiered rows of the JSONL break down as follows:

| tier | rows | rows tied to the summation-polynomial or PDP step |
|:--|--:|--:|
| T2 | 12 | 11 |
| T3 | 33 | 29 |
| T4 | 5 | 3 |
| T5 | 3 | 3 |

The other 24 rows are machinery, problem records and replay records, which carry no tier.

## 1. Evidence machinery

### 1.1 The validation gate

- **Sources:** `docs/VALIDATION_GATE.md` and `scripts/validation_gate.py`.
- **What admits a row:** a `VALIDATION.json` in the row's evidence directory that pins three things by sha256.
  1. A replay certificate, with `status` verified or pass, `independent: true`, `scalars_verified: true` and a non-empty list of runs.
  2. A validator report, with `status: pass`, `independent: true`, `raw_data` and a `headline` list of `{claimed, recomputed}` pairs within `rel_tol`.
  3. A rho reference matched to the protocol: the strongest available, in the same unit and binary, with the 0 / 0.1 / 1 GAE band.
- **Precedence:** a failing sweep verdict wins, then a complete `VALIDATION.json`, then a passing sweep verdict, then the pre-gate freeze, which renders as "unvalidated (pre-gate)".
- **Commands** (stdlib only, read-only; prefix them with `PYTHONDONTWRITEBYTECODE=1`):
  - `python3 scripts/validation_gate.py check [--base origin/main]`. This is the CI check. `--base` also refuses any freeze file that has grown.
  - `python3 scripts/validation_gate.py evidence DIR` explains one directory.
  - The `freeze` subcommand writes files, so it must not be run.
- **State on `af663ad7c`:**
  - `check` prints `validation gate: every row is validated or frozen pre-gate` and exits 0. Re-run in this survey.
  - The leaderboard's 28/28 board rows are `unvalidated (pre-gate)` (`docs/ic/leaderboard.json` `validation_gate.status_counts`), each with the reason "VALIDATION.json is missing".
  - `benchmarks/PRE_GATE.json` freezes 116 rows, all pre-gate.
  - **No real `VALIDATION.json` exists.** The only ones are `scripts/fixtures/validation_gate/positive`, for which `evidence` prints `valid`, and `negative-no-rho-reference`, for which it prints `rho_reference: missing`.
  - The only replay certificate outside the fixtures is `research/f6_ic_geometric_closure_20261003/small_cold/replay-certificate.json`, for the n9 control. It lacks the `independent` and `scalars_verified` fields, and `evidence` reports "VALIDATION.json is missing".
- **What the gate covers:** only `leaderboard.json["board"]`, that is Tables A–C, and `benchmarks/README.md`.
  - Every diagnostic panel sits outside it: n9 F6-IC, the F6-IC ladder, F4/F5/F6, n37, n41/n53, LA engines, pair-table, planner, conductor gap and N83.
  - Those panels are rendered as they stand.

### 1.2 The validation sweeps (`research/validation_sweep_20261019`)

- **Slice A**, `A-leaderboard/verdicts.jsonl`, has 55 records:

  | verdict | records |
  |:--|--:|
  | validated | 1 |
  | validated_with_caveats | 45 |
  | unvalidated | 8 |
  | contradicted | 1 |

  - The single validated record is the rho-optimized Callgrind row.
  - The 8 unvalidated records are three F6-IC ladder rows, four F4/F5/F6 rows and the n61 online bar.
  - The contradicted record is the frontier root-finder replay block.
  - Every Table A/B/C number recomputed exactly, and about 95 sampled scalars verified independently.
  - **The gate consumes none of these verdicts.**
    - `load_verdicts()` reads lines keyed `row`, but slice A writes `row_id`.
    - The verdict words `validated_with_caveats`, `unvalidated` and `contradicted` are in neither `VERDICT_PASS` nor `VERDICT_FAIL`.
    - The lines lack the three artifact keys.
    - A call made in this survey confirms it: `load_verdicts('leaderboard')` returns 0 entries.
- **Slice B**, `B-benchmarks/verdicts.jsonl`, is **not on main**.
  - The strike-through corrections in `benchmarks/README.md` link it on branch `claude/validation-sweep-b-benchmarks-20261019`, for rows 44.2, 104.1, 108.1 and 110.1.
  - The gate loads 0 benchmark verdicts.
- **Slice D**, `D-known-issues/verdicts.jsonl`, has 7 items: 3 medium and 4 low.
  - Item 1 is the root-finder replay mismatch.
  - Item 2: F4 `S` values recorded before #1774 are lower bounds. Re-run on main, the ladder's n13 F4 `S` moves from 16.27 to 488.89 and cross-family K13 from 16.49 to 634.12.
  - Item 3: `quartic_ic.rs` merges Pic⁰ classes.
  - Items 4–7: the KS-tie bug, the synthetic volcano study, the two N83 baselines, and the expired-deadline `verified=true` result.
- **Root cause of item 1, isolated on 2026-10-20** and replicated in `item1_cross_isa_20261020` and its rerun:
  - `inherited_f4.rs:942` `wordwise_rewrite()`, the BMI2 `pext`/`pdep` path, charges different word-XOR counts from the portable scatter path.
  - The method identity does not record which path ran. With a one-line BMI2-off patch, x86 reproduces the Apple counts on 168/168 replays.
  - The bug is **still present on main**. Counted F4 work is therefore host-path specific.

### 1.3 Replay and certificate mechanics (`docs/ecbench/README.md` §9)

- **`ecbench verify --dir SESSION --replay N | --replay-all [--out receipt.json] [--exit-code]`**
  - It re-executes measured runs in-process and requires bit-identical answer, totals, phase counts, counters, unpriced work and factor base.
  - It writes only `--out`, and the receipt's own sha256 is the "replay certificate" a claim cites. Running it from a scratch directory with `IC_ARTIFACT_CACHE=off` left the repo clean.
- **`ecbench claim build … --independent-receipt FILE --pointer WHERE` and `ecbench claim check --report claim.json`**
  - These build and check the `vs_rho` claim schema.
  - A PASS certifies pairing, timing fields and provenance, not a statistically established runtime gain.
  - The independent receipt must come from a different measurement environment and a different hardware class.
- **`ecbench bound check --root REPO --record docs/bounds/records/X.json`** re-derives a frontier bound from its sessions.
- **Replays are binary-bound.** A method's counts are frozen by its id, so a session replays only on a binary that registers its methods. Section 5 has an example of a planner session failing for this reason.
- **CI** (`.github/workflows/ecbench.yml`):
  - Pull requests run an integrity audit of every committed session with `--replay 0`.
  - A weekly lane replays every run of `research/ecbench_*/sessions` and others.
  - The workflow itself states that "PR correctness receipts are CI diagnostics, never full independent-replay certificates for claims".

## 2. Tier map of index-calculus rows

| row(s) | standing | sweep A verdict | tier | PDP / summation step |
|:--|:--|:--|:--|:--|
| Table A, Koblitz one-target, 6 rows | pre-gate | validated_with_caveats ×6 | **T2**: 3 seeded native replays per row, an independent GF(2ⁿ) scalar verifier, and 408/408 per-run `vs_rho` claim-checks PASS with `independent_validation: true` (`rule/v3/claims/summary.json`; scalar and pairing only, not cost) | relation collection by the PDP3 pair table (generic); no Gröbner |
| Table B prime, 8 rows | pre-gate | validated_with_caveats ×8 | **T2**: the whole ladder replays natively with identical counts and matched rho S | m = 2 MITM negation-folded pair table (generic); Relations (PDP) is 68–93 % of S |
| Table B random binary, 5 rows | pre-gate | validated_with_caveats ×5 | **T2** | the same generic oracle |
| Table C, Koblitz 32-target, 9 rows | pre-gate | validated_with_caveats ×9 | **T2**; no archived targets, so no independent [d]G = Q | generic |
| n9 F6-IC control (inherited F4, F6-IC, matrix F5) | outside gate | validated_with_caveats | **T2**: a replay certificate exists, but no validator report and no rho | **algebraic Boolean summation solve** |
| F6-IC E₀ ladder, m = 13/19/23 | outside gate | unvalidated ×3 | **T3**: 24/72 replays per size on main | **algebraic** |
| F4/F5/F6 cross-family, K13/R13/K19/R19 | outside gate | unvalidated ×4 | **T3**: 16/48, 24/40, 16/32 and 8/24 replays | **algebraic** |
| Root-finder vs inherited-F4 session claim-n13 | outside gate | contradicted | **T3** for the IC arms: 72/168 replays, reproduced here | **algebraic vs root-find** |
| n37 K8/K16, 4 panels | outside gate | validated_with_caveats ×4 | **T2**: 320 independent replays; the CI all-run replay of the shared-rank session also reproduced here, 15/15 | PDP3 MITM (generic) |
| n41 shared rank / n53 shared rank | outside gate | validated_with_caveats | **T2** (n41) / **T3** (n53: no IC run verified) | generic |
| N83 Tier 1 cold holdout | outside gate | validated_with_caveats | **T2** | — |
| n61 online bar / n83 dense K1798 online bar | outside gate | unvalidated / validated_with_caveats; sweep B contradicted both | **T4** (n61) / **T3** (n83), both superseded | — |
| LA engines; pair-table set-up | outside gate | validated_with_caveats / not swept | **T3** | — |
| Planner validation, 6 prime curves | outside gate | not swept | **T3**: same-binary self-audits of 12 replays each | the rootfind alternative is 150–637× the plan |
| Prime frontier IC bounds, 2 | frontier | — | **T2**: re-derived here; provenance audit of 12/12 replays | generic |
| 116 benchmarks rows | pre-gate | sweep B not on main | mostly **T3** | see the JSONL `bench-*`/`prime-*` rows |

- **Summation-polynomial-step results at T2.** Only three rows exercise a genuine algebraic summation-polynomial solve at T2:
  1. the n9 F6-IC control;
  2. the PKM tower solving degree at small N, where msolve and `f4_fp_tower` agree on 16/16 systems and five msolve rows were replayed exactly here;
  3. the inherited-F4 determinism finding itself.
- **Algebraic arms at T3.** Every larger Boolean F4/F5/F6 arm is T3, because none replays on current binaries: the F6 ladder, F4/F5/F6 and claim-n13.
- **Generic PDP rows at T2.** Tables A–C, n37, n41, the prime frontier bounds and the PKM-tower ecbench sessions are T2, but their PDP step is generic (a table or root-finder), so dreg never enters a ranked row.

## 3. The prime-field picture

### 3.1 Factor bases defined in the repo

| base | definition | where it is used | effect |
|:--|:--|:--|:--|
| Small-x | Points with x = 1, 2, 3, … where x³ + ax + b is a quadratic residue, canonical y (`ec_index_calculus.rs::build_factor_base`). In ecbench it is `prime-abscissa:size=K`. | Table B; the planner (K = 113–469); the a = −3 ladder (fb = 120) | the default |
| zeta-orbit, Eisenstein-smooth | `ec_index_calculus_j0.rs`; `glv-orbit:size` folds the j = 0 automorphisms | j = 0 curves | divides by \|Aut\| = 6; no exponent change |
| PKM Kummer tower | `prime-kummer:t`: V = μ_{2^t}, the zeros of x^{2^t} − 1; needs 2^t \| p − 1 | PKM rounds | — |
| PKM 2-isogeny tower | `prime-isogeny:t`: V = x(R′ + ⟨Q⟩), Q of order 2^t on an auxiliary E′/F_p; the curve search is charged | PKM ecbench | — |
| Dickson split / Dickson torus | V = {w + 1/w : w ∈ u·μ_{2^t}} in the norm-1 torus; needs 2^t \| p + 1 | PKM round 8 | — |
| Multiplicative subgroup | — | census only | — |
| Dickson torus at P-256 size | `p256_factor_base --factor-base dickson-torus:depth=18` | — | a factor-base inventory only, not an `ic.pipeline` measurement |

- **The PKM decomposition system** keeps one auxiliary variable per tower level, so N = m·t unknowns over a 2^N-point quotient ring. The m ≥ 3 systems use the chain presentation.
- **Subfield bases:** F_p has none. Embedding into F_{p^n} costs Õ(p^{2−2/n}), which is closed on paper.
- **Solution-count census** (`ic_next10_20261014/02`, T3): every algebraic base gives 0.93–1.03× the decompositions of a random base, over 172 cells.

### 3.2 Degree of regularity: observed against predicted (PKM tower)

**Predicted.**
- The pre-registered H0 was that D grows linearly in N.
- The Bardet–Faugère–Salvy degree of regularity is **N + 1**, because the top-part ideal has Hilbert series (1 + t)^N.
- XL/Macaulay gives **D_XL = N + 4 exactly for m = 2**: δ = N, with every one of the 2^N square-free coefficients of 1/S̄ non-zero. That runs from 16 at N = 12 to 30 at N = 26.
- A direct Macaulay rank gives 13 at m = 3, N = 9; at least 14 at m = 3, N = 12; and 12 at m = 4, N = 8.

**Observed F4 solving degree D:**

| m | N | D | tier |
|--:|:--|:--|:--|
| 2 | 2–20 (split and torus) | 3 4 4 4 4 5 5 5 5 6 | T2 at N = 12–18, T3 elsewhere |
| 2 | 22 | 6 | T3 |
| 3 | 6 / 9 / 12 / 15 / 18 | 5 / 6 / 6 / 7 / 7 | T2 at N = 9 and 12; T3 at 15 and 18 |
| 3 | 21, 24 | ≥ 7, inconclusive | T3 (stopped runs) |
| 4 | 8 / 12 / 16 | 6 / 7 / 7 | T2 at N = 8 and 12; T3 at 16 |
| 4 | 20 (msolve) | **≥ 8** | T3 |

- **Reading.** D sits far below both predictors; F4 escapes the XL bound through degree falls. It grows about one degree per eight tower unknowns over N = 12–20 at m = 4, so the plateau at 7 is refuted.
- **The m = 3, N = 21 wall:**
  - msolve with `-t 16` under a 108.8 GB cap exhausted degree 7 without the unit ideal, selected 398,409 degree-8 pairs, and died at 99.6 GiB.
  - The project records this as **inconclusive**.
  - *My interpretation, labelled as such:* if the system is inconsistent, its refutation needs a step of degree ≥ 8. That is the same argument the m = 4, N = 20 cell used.
- **Projected resources** for `f4_fp_tower`: about 130–150 GB for m = 3, N = 21, needing one process with at least 192 GB; 240–380 GB for m = 4, N = 20.
- **Engine checks:**
  - Signature F4 (F5/GVW) needs D_sig > D on 18 of 18 systems and 1.7–28× the memory.
  - Magma has never been run on these systems (T5).
  - v₂ ≥ 40 is blocked because the engines store 32-bit coefficients and require p < 2³².

### 3.3 Best measured cost against rho

| result | ratio to rho | tier |
|:--|:--|:--|
| Table B, m = 2 MITM (best row) | 5.47× at fp12 (log₂ r 11.9) and 13.5× at fp24 (log₂ r 23.4); the relation phase alone is 3.72–9.63× (calculated from the phase S values) | T2 |
| Planner | 12.6–25.8× at log₂ r 23.5–29.6 | T3 |
| PKM tower, ecbench, 20–28 bits | 29–141×; the tower base costs the same as a plain base (0.999–1.012) | T2 |
| PKM isogeny tower | 17–137× | T3 |
| Frontier, prime toy tier | ops 22.7× the floor; α = 0.84 [0.50, 1.02], against 0.40–0.45 for rho with the negation map; dominated by every rho and BSGS entry | T2 |
| a = −3 ladder | 19× → 403× wall over 16–24 bits | T3 |
| IC-3 at 32 bits | 20.5 s against 0.5 s for rho | T3 |
| secp160r1, S₃ solve per target | 79 bits worse | T3 |
| ECCp-131 two- and three-point probes | 0 relations; the method's coverage bound is 2^−111 per candidate at B = 512 | T3 |

**The gap is structural.**
- Every generic decomposition obeys S·T² ≳ n.
- m = 2 never beats rho and m = 3 at best ties. The first arity that can win is **m = 4**, which needs an S₅ small-root oracle at δ ≈ 1/5.
- Measured lattice reach:
  - S₃: Howgrave–Graham holds to δ = 0.10; the heuristic fades between 0.12 and 0.16. T3.
  - S₅: δ* = 0.0125, rising to 0.01398 with shared targets, under a proven ceiling of 1/32. T3/T4.
- The D₄-symmetrised S₅ lattice was **never run** (T5).
- Engineering moved constants by 10³ (PR #1560), and root-finders beat F4 per decomposition by 10⁹–10¹¹. Neither changes the exponent.

## 4. Binary and ECC2K-130 (summation-polynomial and dreg parts only)

- **Online cost is probe-rate bound.** The S₃ quadratic solve per probe dominates, at about 1.1 M probes/s on one core at n = 83. T3.
- **Higher-arity extraction is unbuilt.** The m = 5 extraction and the 6-sum S₄-per-probe design have never been built (T5). Exact S₅ machinery exists: 22,887 clauses at n = 53.
- **The 4-sum total work does not cross rho** (T4 extrapolation): r/(n²K) probes, which is about 2·10⁷ × rho at K = 10⁹.
- **RESEARCH_ECC2K130_IC_FEASIBILITY.md carries no dreg data.** The pricing lives in `research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` (T4 on T3 cells):
  - m = 4 costs 2^{+15.50} over rho if D stays at 6, and 2^{+112.32} if D follows ℓ + 4.
  - msolve F4 degree 4, where XL reads ≥ 7, gives +7.6 bits in the model before the K conversion.
  - The (4,6) cell needs about 70 GB and is unmeasured.

## 5. Problem records

**Records related to dreg, summation polynomials, PDP, Gröbner solving and reproducibility:**

| record | status | what would unblock it |
|:--|:--|:--|
| `eccp131-full-parameter-index-calculus` | open (10-22) | an oracle or base whose independent-row rate escapes the 8B³ support bound, then an independently replayed relation → rank → descent chain verifying [x]G = Q with full cold cost |
| `binary-boolean-system-64-variable-cap` | open | a monomial type for at least 128 unknowns in `build_decomposition_system`, inherited F4 and the SAT encoder, agreeing with the floor on 12 targets at m = 4, n = 23, d = 6 |
| `generalized-rootfinding-all-roots-8x` | open | complete canonical root output; 20 preregistered hard cases (baseline ≥ 1 s, 3 cold repetitions, A/A ≤ 10 %, worst slowdown ≤ 1.10); 8× against the fastest external F4/F5/F6; independent replay |
| `k5-joux-vitse-trial-cost-and-fpk-ecbench` | open | C_trial ≤ 10⁵ F_p multiplications (it is now degree-19 F4 at 5.9–9.2·10⁹) or a different exponent; an ecbench `F_{p^k}` instance; registry harvesting for extension curves |
| `rootfinding-complete-cost-calibration` | open | a frozen protocol with 8 targets × 24 executions per arm and `unpriced: []`, L2 isolation with A/A, and independent replay on a second host class |
| `release-suite-jv-nagao-wall-budget` | investigating | reproduce and instrument the parallel timeout without weakening the budgets |
| `n83-online-aa-noise` | investigating | a mirrored rescreen with adjacent A/B and C/D pairs within 10 %, plus a hashed IC1 manifest |
| `n83-recorded-query-sampler-contract-mismatch` | open | trace the streams, version any change, pass the regression |
| `ic-improve-04-genus3-2lp-study-missing` | investigating | the original tool and records at a pinned revision, and a matched comparison |
| `curve-traits-debug-subfield-trace` | open | reproduce, then fix the assertion without weakening it, then a byte-identical replay |
| `plane-quartic-jacobian-class-arithmetic-missing` | "resolved (toy smokes)"; sweep D says it should be partial | re-solve on canonical Pic⁰ arithmetic (#2227) |

**Two major blockers have no record**, even though AGENTS.md's "Name the obstruction" rule calls for one:
1. **The inherited-F4 BMI2 counter path and the IC replay regression.** It is tracked only in sweep D, PR #2290 and the leaderboard's "Open regression" row.
2. **The PKM/descent solving-degree memory wall.** It needs one process with at least 192 GB, or about 70 GB for the binary (4,6) cell. This is a resource limit, not a mathematical one.

## 6. Replays run in this survey, with outputs

Here `$S` is a scratch directory and `R=/home/user/aburan28/crypto`. All of these are cheap and deterministic: under a minute on 4 vCPU.

1. **Validation gate.**
   - `cd $R && PYTHONDONTWRITEBYTECODE=1 python3 scripts/validation_gate.py check` prints `validation gate: every row is validated or frozen pre-gate` and exits 0.
   - `... evidence scripts/fixtures/validation_gate/positive` prints `valid`.
2. **msolve PKM D replay.**
   - Command: `msolve -t 1 -v 2 -f $R/research/ic_next10_20261014/03-pkm-degree/systems/NAME.ms -o $S/NAME.out > $S/NAME.log`.
   - Parse with `awk '/^ +[0-9]+ +[0-9]+ +[0-9]+ +[0-9]+ x [0-9]+/ {n++; d=$1; if(d>D)D=d; if($4>R)R=$4; if($6>C)C=$6} END{print n,D,R,C}'`.

   | system | rounds | D | max rows × cols | output | wall |
   |:--|--:|--:|:--|:--|--:|
   | `pkm_kummer_p2013265921_m2_t6_x464833983` | 14 | 5 | 2803 × 2064 | `[-1]` | 0.08 s |
   | `…_m3_t3_x388393539` | 17 | 6 | 3827 × 2997 | — | 0.19 s |
   | `…_m4_t2_x92234565` | 22 | 6 | 2009 × 1906 | — | 0.08 s |
   | `…_m3_t4_x334226296` | 19 | 6 | 12900 × 12287 | — | 8 s |
   | `…_m4_t3_x259303215` | 30 | 7 | 26025 × 21665 | — | 59 s |

   All five are identical to `results.csv`.
3. **ecbench replay of the PKM tower prime sessions.**
   - Command: `cd $S && IC_ARTIFACT_CACHE=off RAYON_NUM_THREADS=1 ecbench verify --dir $R/research/pkm_tower_ecbench_20261013/sessions/p20 --replay-all --out $S/p20.json`.
   - p20 prints `audit OK: 192 records, 192 verified, 144 replays reproduced, 0 problems` in 5.9 s.
   - p28 prints `160 records … 120 replays reproduced` in 15.5 s.
4. **Root-finder claim-n13 replay, a negative control.**
   - Command: `ecbench verify --dir $R/research/ecbench_rootfind_engine_20261012/sessions/claim-n13 --replay-all --out $S/rf.json`.
   - It prints `audit FAILED: 224 records, 224 verified, 72 replays reproduced, 96 problems` in 8.4 s.
   - All 24 rows of each of `ic-f4-m2` and `ic-f4-m3` differ in total_gae, phase counts, gae and unpriced.
   - All 24 rows of each of `ic-rootfind-m2` and `ic-rootfind-m3` differ in phase counts only.
   - All 72 rho rows are identical.
   - This matches sweep A's result on current main, which means the drift was already present at crypto `0dc9ef4`.

Extras:
- `ecbench verify --dir $R/research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02 --replay-all` prints `audit OK: 18 records, 18 verified, 15 replays reproduced, 0 problems` in 0.96 s. This is a cross-ISA replay that does reproduce: the session was recorded on Apple arm64 and has no inherited-F4 arm.
- `ecbench bound check --root $R --record $R/docs/bounds/records/prime-ic-{mitm,subtract}.json` prints `re-derived ECBND1hbab6d6dc37fc OK` and `re-derived ECBND1h9d3beb6852c6 OK` in 0.3 s.
- `ecbench verify --dir $R/research/ic_planner_validation_20261020/sessions/r24 --replay 12` reproduces 4/12. Every failure is "oracle `mitm-prime-blocked` does not run on a prime curve": a binary-version limit, so these sessions need a current-main build.
- Table B's native replay needs the `ic` binary, which was not built here: `ic boundary --regime prime --prime-bits 10,…,24 --repeats 3 --seed 123212651130` (sweep A `scripts/replay_B.sh`).

## 7. Caveats

- **Tier boundary.** The T2/T3 line depends on what counts as independent. Same-host, same-binary self-replays were kept at T3, and the sweep's seeded replays were counted as T2.
- **Host class.** The replays in this survey ran on the same Xeon host class as several original runs, but with different binaries and containers. They are not the different-hardware-class receipts that ecbench claims require.
- **Sweep B.** Its verdicts were not readable (branch not on main, and the repo is not attached to this session's GitHub access). Benchmark-row verdicts are therefore taken from the README annotations only.
