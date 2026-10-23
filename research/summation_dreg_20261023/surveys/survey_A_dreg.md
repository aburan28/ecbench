# Survey A: levers on the degree of regularity in the PDP (summation-polynomial) step

Slice A covers the dreg / closure / XL studies, the semi-regular law, volcano regularity, XL bilinear/trilinear, F4 early stop, WDSat, Hamming ideals, and the `crossbred` / `degree_reduction*` modules. The repo `/home/user/aburan28/crypto` was only read. I followed AGENTS.md's reading discipline: located files with `rg`, checked sizes, read raw JSONL only through `jq`, and read sessions through their RESULTS/score files. The one exception was a build of examples into a target directory in my scratchpad, outside the repo; `git status` was clean afterwards.

Machine-readable version: `survey_A_dreg.jsonl` (44 objects; the fields are listed at the end). The rerun script is `repro_A.sh`. Both files are in this directory.

## How to read "degree" here

Every number below is a **refutation degree**: the first Macaulay degree D at which the constant 1 lies in the row space of an unsatisfiable (non-decomposable) instance. The engine sets which object that is:

* **XL** (plain bounded Macaulay): only the original equations times monomials.
* **MutantXL closure, F4 (msolve) or the m=2 mutant closure**: degree falls are multiplied again (fall reuse).
* **"FFD" in the dreg_* tables** is a *rank proxy*, not a certified first fall degree (degree_reporting_20260922).

## Headline findings

1. **The lever that worked was fall reuse.** The MutantXL closure (A01) refutes all 92 committed chained-S3 draws at degree **3–4**, where XL needs **5–7** (m=3,4; ℓ=1–5; N=11–32). The closure matrix is 2^1–2^7 narrower, and the gap widens with ℓ. msolve F4 (D07) reads 3–4 on 48 of the same draws, draw-for-draw identical to the closure on the 28 m=4 ℓ≤2 draws (10 of them at 3). The m=2 mutant closure (D13) gives d_F=3–4 against XL 5–8. The rising XL ladder (B01–B05: 5, 6, 6–7, 7 over ℓ=2..5, and ℓ+4 at m=4) therefore measures XL, not the ideal.
2. **Fall reuse has a limit.** Degree 4 holds to N=54 at m=4, ℓ=6 (A04, 2/2 draws) and at m=3 up to ℓ=6 (A02, A03). At **m=3, ℓ=7 (N=40, 42) the degree-4 closure saturates without refuting** on 8/8 draws (A02). Of the two surviving explanations, "m·ℓ≥21" was refuted by A04 and "ℓ=7" survives. Degree 5 at ℓ=7 is unmeasured: two hosts failed to allocate memory, at 169 GB and 405 GiB.
3. **Fixing variables brings ℓ=7 back to degree 4** (A05). Fixing one coordinate per summand refutes at 4 on 8/8 draws. Any single fixed variable does it on 6/8, and both branches of a single guess refute on 5/8. The registered question went unanswered because the null control also refuted.
4. **Rewriting the presentation helps.** The symmetric norm form `rr` refutes 3–4 degrees below the direct S4 descent `x4` on Koblitz K0/K1 (n=23–61; D01). For example, 6 vs 10 at K1/2^61 ℓ=6, with a 39,966×43,796 matrix in 24 s against 246,928×199,140 in about 45 min. Symmetrisation inside the degree-reduction thread lowers D* (5 vs 8–9 at ℓ=3) but widens the matrix (C04).
5. **Levers that did not lower the degree, or lowered it at no net gain:**
   * hybrid slicing (C01) and mutants plus guessing (C03): D* falls but the cost optimum degenerates to exhaustive search;
   * curve and isogeny moves (C05, and D10's isogeny arm): the degree is invariant across the class;
   * Hamming-weight ideals (D09): bounded generator degree, but call exponents of 1.9–2.5;
   * trilinear XL (D04);
   * the trace-zero pre-filter (D12): it never fires on odd-order IC targets;
   * the engineering paths (B07–B10, D08): same degree, lower cost.

   The F4 early stop (D05, prime-field Joux–Vitse cover) removes the last degree step (8 → 7) and cuts counted F4 multiplications by about 4–5× in the chunked Hilbert mode. Its performance gates are still open.
6. **Semi-regular prediction.** It predicts **XL** within ±1 in every cell (equal in 10/14 m=3 cells and 4/7 m=4 cells; B06). For x4 the affine d_aff is an upper bound, not exact (D02, round 2 unscored). It **overstates** the fall-reusing solvers by 3–5 degrees.
7. **No n=131 descent row beats rho** in any of these studies (E01, T4). The free-oracle floor binds independently of the degree.

## Tier counts (44 entries)

| T1 | T2 | T3 | T4 | T5 |
|--:|--:|--:|--:|--:|
| 1 | 11 | 24 | 1 | 7 |

* T1 is volcano (validation.json plus rank certificates; I reproduced it bit-identically here).
* T2 entries are either replayed with matching outputs (XL ladder identity checks over 299–376 committed rows; engineering identity checks; the n=7 cross-version build) or corroborated on the same draws by an independent engine (A01, via msolve F4 and the second sparse implementation). Rerun status in this session is in the "Six cheapest reruns" table below. Tiers record the repo's own evidence. My matching replays (JSONL field `replayed_in_survey`) cover A01, B02, B05, D01 and D10 in full, and one draw each of A02, A03 and A05. Counting those, A02, A03, A05 and D01 would also qualify as T2.

## Technique table

Verdicts concern the degree of regularity / refutation degree unless the entry says otherwise.

| id | technique | system | key parameters | measured | verdict | tier |
|---|---|---|---|---|---|---|
| A01 | MutantXL fall closure vs XL | binary Koblitz E0, chained S3, m=3,4 | 92 draws, N 11–32, ℓ 1–5 | closure 3 (16) / 4 (76) vs XL 5–7; gap 1 (26), ≥2 (58); 0/345 sat refuted | **lowers** 1–3 | T2 |
| A02 | closure at scale | E0 chain | N 26–42, ℓ 5–7, m=3,4, deg 4 | 32/40 refute at 4; m=3 ℓ=7: 8/8 saturate without 1 (rank 94,230/102,091) | lowers, bounded at ℓ=7 | T3 |
| A03 | closure driver (N vs ℓ vs m·ℓ) | E0 chain | m=3 ℓ=6 N 40,42; m=4 ℓ=5 N 42,44 | 16/16 refute at 4; XL 0/16 | lowers; N rejected | T3 |
| A04 | closure at (4,6) N=54; degree 5 at (3,7) | E0 chain | 341k cols; 184/471 GB caps | A1 2/2 refute at 4 (25.8 GB, ~24 min); C1 allocation failure twice | lowers/holds; deg-5 inconclusive | T3 |
| A05 | fix 1 or 3 vars, then degree-4 closure | E0 chain | 8 ℓ=7 draws, 4 arms | 3-fix 8/8; single fix 6/8; both branches 5/8 | lowers per branch (2^k cost) | T3 |
| A06 | sparse second closure implementation | E0 chain | 36 small comparisons | agrees; too slow at scale | neutral (instrument) | T2 |
| A07 | XL vs rank-proxy FFD, m=2/3 small | E0 (invariant subspaces) | n 5, 7 | XL 6 vs FFD 3 | baseline | T2 |
| B01 | fixed-surplus ladder (hold S) | E0 m=3 | 4 pairs S −2..+1 | 5,6 → 6,≥7: grows | raises (via ℓ) | T2 |
| B02 | (n,ℓ) grid at matched N | E0 m=3 | N 11–25 | ℓ sets degree; (10,5) ≥7 | raises with ℓ | T2 |
| B03 | lower surplus at fixed ℓ | E0 m=3 | S −5, −4 | (7,4) 7767; (8,4) 6767 | raises (keep S ≥ −2) | T2 |
| B04 | exact degree 7 at ℓ=5 | E0 m=3 | (10,5), 727k cols | 7777; 2.5 h, 11 GB/draw | raises (XL) | T3 |
| B05 | m=4 chain XL | E0 m=4 | ℓ 1–3, N 12–32 | 5 / 6 / ≥7; flat in n | raises (+1 vs m=3) | T2 |
| B06 | semi-regular Hilbert-series predictor | E0 | 14+7 cells | XL within ±1 everywhere | predictor | T3 |
| B07 | structured sparse elimination | E0 | n=5 | 6.03× end to end | neutral (cost) | T2 |
| B08 | dense finish | E0 | 376 rows | ~10×, 0 mismatches | neutral (cost) | T2 |
| B09 | matrix-F5 criterion rows | E0 | deg 7 | −20–26% rows | neutral (memory) | T2 |
| B10 | budgeted dense switch + `--d-min` | E0 | (10,5) deg 7 | fits 11 GB | neutral | T2 |
| C01 | L3 hybrid slicing (guess k) | generic binary m=2 | N 10–14 | D* → 2 but optimum = exhaustive; one-sided c=1/2 | lowers D*, no gain | T3 |
| C02 | L4 degree falls (mutants) | generic binary m=2 | N 10–20 | Random D* 4→2 (working 4→3), +1.44 bits; cascade dies at N≥18; 0/54 beat 2^N | lowers 1 degree, degrades | T3 |
| C03 | L3∘L4 | generic m=2 | N 10–14 | −2.87 bits vs guessing | raises cost | T3 |
| C04 | L2 symmetrisation (S4, m=3) | binary S4 | ℓ 2, 3 | sym D* 4,5 vs elim 8–9; never narrower | lowers degree, raises width | T3 |
| C05 | L5 isogeny / curve move | ordinary binary curves | 2,994 comparisons | 0 improvements | neutral | T3 |
| C06 | Δ_low early-defect screen | generic m=2 | 2,016 cells | ranks families only | diagnostic | T3 |
| D01 | rr symmetric norm vs x4 direct S4 | Koblitz K0/K1, n 23–61, m=3 | ℓ 4–6 | rr 5–8 vs x4 8–10 (3–4 lower) | **lowers** | T3 |
| D02 | semi-regular / d_aff / d_red rules | Koblitz | rounds 1, 2 | x4 in window 32/32; d_aff ≥ reading 42/42 (round 2, my tally) | predictor | T3 |
| D03 | bidegree XL (bilinear core) | E0 family, n 19–31 | 72 cells | reach past n/3; n=131 gap unchanged (71.82) | reach up, cost neutral | T3 |
| D04 | trilinear XL, t unknown | E0 family, n 11–17 | 54 cells | t-degree 2 needed; 87.16 > 71.82 | raises cost | T3 |
| D05 | F4 early stop (Bézout 64 / Hilbert) | prime-field JV genus-3 cover | p 53, 101, 251 | max deg 8→7; chunk1 ~4–5× fewer muls | lowers 1 (termination) | T3 |
| D06 | XOR-Gauss in CDCL (WDSat) | ECC2K family RR encodings, n 5, 9 | 8 targets | 1,277× fewer conflicts | n/a (SAT) | T3 |
| D07 | msolve F4 on committed draws | E0 chain | 56 draws to N 38 | D_F4 3–4 vs XL 5–≥7 | **lowers** | T3 |
| D08 | G2 Frobenius-equivariant pass | Koblitz barrel | n 7–31 | RREF identical; speed corrected to 2–3× | neutral | T3 |
| D09 | Hamming-weight ideals + MultiSolve | Koblitz K0, m=2 | n 7–19 | correct; exponent 1.9–2.5 | neutral/raises | T3 |
| D10 | S3 mod support polys; isogeny transport | F_103 j=0 + neighbour | 480 systems | d_reg_top 3→2 at B=2 only; neighbour 240/240 equal | lowers (tiny) / neutral | **T1** |
| D11 | crossbred (module) | K and R binary | engine suite | crossing at 26 unknowns vs scalar FES only | inconclusive | T3 |
| D12 | trace relation / trace-zero filter | binary m=2 | 112 + 29,624 | degree 2 on wrong-coset targets; 0 fires on IC targets | lowers only off-target | T3 |
| D13 | m=2 last fall degree (mutant closure) | binary m=2 | n 8–22 | d_F 3–4 vs XL 5–8 | **lowers** | T3 |
| E01 | n=131 descent cost rows | ECC2K-130 | — | 2^68.43 … 2^196.95; none below rho | extrapolation | T4 |
| F01–F07 | registered or named, never run | — | — | — | null | T5 |

The F rows are:

* F01: A2/B/C2, A3, and the C1 retry;
* F02: guess count at ℓ=8, 9;
* F03: (15,5) and (9,3) at degree 7;
* F04: a sparse closure for degree 5;
* F05: block-aware multihomogeneous Gröbner;
* F06: fall-degree experiments E3–E5;
* F07: the L2 solver at ℓ≥6 (withdrawn).

## Contradictions, retractions and corrections

1. **The "infeasible control" is satisfiable.** `random_control_system`/`ladder_control` has the zero root, so it can never be UNSAT (degree_reporting_20260922, confirmed in `koblitz_bench.rs::ladder_control`). Two readings were read as refutation comparisons, but they compare a pinning degree of a satisfiable system with refutation degrees:
   * RESEARCH_DREG_MEASUREMENT "structure not shape" (`ctrl(unsat)`);
   * fixed_surplus/ell_grid "the control resolves ≥1 degree above the Semaev draws".

   "FFD" is a rank proxy throughout.
2. **The XL growth was XL's own** (A01 vs B01–B05). The fits ⌈ℓ/2⌉+4 and ℓ+4 and the "grows at fixed surplus" verdict are XL properties. The fall-reusing solvers stay at 3–4.
3. **A missed prior result.** The mutantxl PREREGISTRATION (10-16) says "the coverage survey found no XL-vs-MutantXL experiment". ic_next10_20261014/06 had shown msolve F4 at 3–4 vs XL 5–≥7 on 48 of the same committed draws two days earlier.
4. **The "4, constant" cost row was revised four times.** The 10-16 mutantxl study added it (m=4 2^68.43). The scale study withdrew it for m=3. The driver study called it "contradicted" for m=4. Then m4l6 (10-18) refuted the m·ℓ≥21 explanation and kept degree 4 to N=54. Neither RESEARCH_DREG_MEASUREMENT.md nor RESEARCH_ECC2K130_DESCENT_DEGREE.md was updated after 10-17.
5. **m4l6 "single-threaded" claim.** README and PROTOCOL §8 call the closure single-threaded. `/usr/bin/time` recorded 1353–1429% CPU on A1 and 637–765% on C1.
6. **Memory projections missed in both directions.** The m4l6 model 3·cols²/8 projected 44 GB for A1 (measured 25.8 GB) and 217 GB for C1 (failed at 169 GB and again at 405 GiB). The driver prereg cites a 7.1 GB peak and the m4l6 protocol a 5.8 GB peak, both at 124,314 columns. Both numbers are real: they are different m in the scale study.
7. **Semi-regular law.** Round 1's verdict is "neither". The fixed-ℓ fit `x4 = ℓ+5` was retired, and rr matched the face-value window by accident (duplicated equations). Round 2 was measured but has no RESULTS.md. My tally of round 2 is: rule A fails (12/42 x4 draws read one below d_aff), A′ holds, and rule B fails (K0n53l5 reads 6 against a prediction of 5).
8. **ic_next10 corrections (ic_improve_20261014).** Item 6's 44–73× G2 speedup is 2–3× against `gf2_elim`. Item 8's 50% trace-zero refutation never applies to odd-order IC targets (0/29,624).
9. **The degree-reduction thread corrected itself.**
   * A syzygy was mistaken for a fall in iteration 1.
   * The R4b size-proxy flag on the FFD law was raised and then withdrawn.
   * EXP-R2's "censored at D≥5" was 25× under budget.
   * L4 was classified "relabelling": its gate passed while the 2^N margin worsened by 3.18 bits.
10. **Smaller corrections:**
    * fall_degree_bounds: as run, E1 falsified H1 (58/112) because of a constant mismatch. Amendment 1, made after the run, rescored it 112/112.
    * ic_gb_ladder errata: two "degree-4 at ℓ=5" readings in ic_rr_norm_ladder were constant equations.
    * The WDSat display printed a timeout as UNSAT; the counts were never wrong.
    * DREG Result 2/3: "n=7 did not complete" was retracted, the 494 s timing is unscoped (27× environment), and a circular argument was withdrawn.
    * degree7_ell5: amendment 1's 20.7 GB estimate for (13,5) assumed full rank; the revised figure is 47–50 GB.
11. **Tension, not contradiction.** The m=2 mutant cascade (C02) dies at N≥18 (D*=4 after one round), while the full MutantXL closure (A01–A04) holds at 4 to N=54. They are different objects: one degree-3 saturation followed by XL, versus saturation at every degree, on m=2 versus m=3/4 chains.

## Determinism guide (for reruns)

* **Bit-identical on any host:**
  * refutation or closure degree;
  * per-degree `xl_rank`, `closure_rank`, `columns`, `closure_products`, `saturated`;
  * seeded draw lines (draw, v_basis, x_r, solutions, outcome, ffd), excluding `secs`;
  * the `.sing` system bytes (sha256 in `dump.sha256`);
  * macaulay_dense `n_rows`, `n_cols`, `rank`, `refuted`;
  * dreg_score and score.py verdict tables;
  * volcano `results.json`/stdout;
  * Hamming counted units (calls, tame, wild, xor_words, tame_depths);
  * semireg.py output.
* **Identical on the same host and binary, possibly drifting across ISA or toolchain:** `closure_word_ops` and all inherited-F4 word-XOR counters. Thread count and hash salt were proven irrelevant. An x86 vs aarch64 drift is unexplained (inherited_f4_replay_determinism_20261018). In this survey, HEAD `closure_word_ops` on this x86 container equalled the 10-16/10-17 records on all 92 + 3 draws.
* **Nondeterministic:** `secs`, wall, `ms`, VmHWM/peak RSS, CPU%, queue timestamps. The F4 early-stop counts are deterministic in principle, but the HEAD sources differ from the receipt hashes, so a HEAD rerun may not match the historical files.

## Six cheapest reruns (each ≤10 min on 4 cores / 15 GB)

`repro_A.sh` in this directory implements all six. It writes only under `$OUT`. A crate build (outside the tree, via `CARGO_TARGET_DIR`) comes first and is not included in the times. **Status here:** all six ran once in this session (2026-10-23, this 4-core/15 GB container, HEAD `af663ad7` built out of tree with rustc 1.97.0; outputs in `reproA_run1/`, binary hashes in `binaries_head.sha256`). All six matched the committed records value for value. The only schema difference is that the HEAD `dreg_closure` writes an extra `"fixed":[]` field on unspecialised runs. `closure_word_ops` also matched the committed records, which came from older binaries on other hosts. The runtime column below gives the times measured here.

| # | lever | command (from repo root) | runtime / memory | expected deterministic output | status here |
|---|---|---|---|---|---|
| R1 | D10 encoding rewrite (T1) | `cp research/volcano_regularity_control_20260924/measure.py $OUT/ && python3 -B $OUT/measure.py > $OUT/run.log; sha256sum $OUT/results.json $OUT/run.log` | 0.13 s | `198bcef6ca7520f134088adc98b6a13c530eeefc5b2c121df465619d249a5a6b` results.json; `edd90185df3c7404b08ed4aaffaf29c55aa5619a6550d6fd8897269090f84fb4` run.log | **matched** |
| R2 | B02/B05 XL degree | `F4_F2_MAX_ROWS=50000000 F4_F2_MAX_COLS=50000000 $BIN/dreg_ladder --cells 12:2:7 --ffd-max 5 --seed 20260925 --unsat 4 --controls 0`; same with `--cells 5:3:7`; `KIC_SPARSE_DENSE_FINISH=1 KIC_SPARSE_F5=1 KIC_SPARSE_DENSE_BUDGET_MB=11000 $BIN/dreg_ladder --m 4 --cells 6:2:8 --ffd-max 5 --seed 20261006 --unsat 4 --controls 0` | 20 s here (14 + 3 + 3 s), <0.2 GB | (12,2): 4 rows, refuted at 5, FFD 3333. (5,3): 32 rows, unsat draws 9/14/22/31 at 6. m=4 (6,2): draws 0–3 at 6. Each line equals the committed `runs/` line except `secs` | **matched** |
| R3 | D01 rr vs x4 | `$BIN/rr_degree_ladder --a 0 --n 61 --ell 4 --unsat 4 --seed 20261005 --dump-dir $OUT/r3 --out $OUT/r3/K0n61l4.catalogue.jsonl`, then `macaulay_dense --file … --degree D --threads 4 --mem-gb 12` for D=1..9 | ~2 s here, <0.1 GB | 12 `.sing` sha256 equal to `runs/registered/dump.sha256`. rr first refutes at 5 (2801×2380, ranks 525/516/525/517); x4 at 8 (4819×3797, ranks 2509/2504/2508/2514) | **matched** |
| R4 | A01 MutantXL vs XL | loop over the 41 inputs listed in `research/dreg_closure_mutantxl_20261016/run.sh`: `$BIN/dreg_closure --max-cols 70000 <input>`, then `$BIN/dreg_score closure $OUT/*.jsonl` | 126 s here (4 m 15 s recorded), ≤1.0 GB | score table equal to `score-output.txt`; 92 draws; closure degree {3:16, 4:76}; Σclosure_products 4,902,150; Σclosure_rank 60,470; Σxl_rank 165,625; Σword_ops 101,006,104,256 (same-host) | **matched** |
| R5 | A02 + A05 | `prlimit --as=15000000000 $BIN/dreg_closure --to 4 --max-cols 130000 --draw 0 research/dreg_closure_scale_20261016/draws/m3-19-7.jsonl`; again with `--fix 6=0,13=0,20=0` | 130 + 60 s here, 3.6 / 3.1 GB | rise: closure_degree null; d4 xl 16,136, closure_rank 94,230, cols 102,091, products 533,998. fixB: closure_degree 4; d4 rank 1, cols 74,519, products 472,858 | **matched** |
| R6 | A03 | `prlimit --as=15000000000 $BIN/dreg_closure --to 4 --max-cols 155000 --draw 8 research/dreg_closure_driver_20261016/draws/m4-11-5.jsonl` | 84 s here, 3.5 GB | closure_degree 4; d4 xl_rank 10,776, closure_rank 1, cols 123,814, products 421,742 | **matched** |

**Build at HEAD.** HEAD is `af663ad7` (1-commit shallow clone, clean tree). An `--offline` build fails because `rusqlite` is not in the local registry. An online `CARGO_TARGET_DIR=<outside> cargo build --release --locked --example dreg_closure --example dreg_ladder --example dreg_score --example macaulay_dense --example rr_degree_ladder` succeeded in 13 min 09 s. `git status --porcelain` stayed empty; only `.git`'s index cache was refreshed.

**Other cheap reruns (not run here):**

* `python3 -B research/dreg_degree7_ell5_20260930/semireg.py` (bit-identical, <1 s);
* the Hamming PDP CI replay (`python3 replay_check.py --n 7 --encodings SUB,FC --targets 4`, after `cargo build --release --manifest-path hamming_pdp/Cargo.toml`; expect `replay differences: 0`). Copy the study directory out first, because the crate builds into its own `target/`;
* msolve F4 on small committed draws (D07, 0.1–4 s per draw);
* draw-generation replay `dreg_ladder --cells 15:6:3 --m 4 --ffd-max 3 --seed 20261018 --shard 2/4 --first-draw 40 --max-draws 44` (~1 min; reproduces A1 draw 42).

## Extras found by `rg` (in slice)

* `research/notes/index-calculus/RESEARCH_DREG_MEASUREMENT.md` (Results 1–14, the index of this line);
* `research/notes/index-calculus/RESEARCH_DEGREE_REDUCTION.md` (levers L1–L5 for degree_reduction*.rs) with `experiments/degree_reduction_*.json`;
* `research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md` (cost rows);
* `research/ic_next10_20261014/06-equivariant-f4-m4` (msolve F4) and `08` (trace-zero), with corrections in `research/ic_improve_20261014`;
* `research/fall_degree_bounds_20261008` (m=2 last fall degree, trace relation);
* `research/degree_reporting_20260922` (control and FFD correction);
* `research/ic_leads_rented_20261018` (A3/C1 continuation of m4l6);
* `research/inherited_f4_replay_determinism_20261018` (counter determinism);
* `problems/binary-boolean-system-64-variable-cap.md` (`build_decomposition_system` caps N at 64);
* `research/ic_framework_engines_20260922` (crossbred measurements).

Pointers only, owned by other slices: `ic_rr_norm_ladder_20260930`, `ic_gb_ladder_20261003`, `ic_dense_ladder_20261004` (rr/x4 predecessors); the `pdp_split_*` directories; the `pkm_*` directories (prime-field PKM msolve); and the `toy_*` isogeny audits (they agree with C05).

## JSONL fields

`id`, `name`, `lever_class`, `mechanism`, `curve_family`, `parameters`, `measured_outcomes`, `vs_baseline`, `vs_semiregular`, `verdict`, `confidence_tier`, `tier_rationale`, `retractions_contradictions`, `artifacts`, `reproduce_recorded`, `est_runtime_4core_15gb`, `deterministic`, `nondeterministic`, `replayed_in_survey`.
