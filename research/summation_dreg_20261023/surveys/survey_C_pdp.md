# Survey C: PDP splitting, oracles, root-finding vs Gröbner, SAT/SMT/WDSat/MITM, binary Weil-descent models

Machine-readable inventory: `survey_C_pdp.jsonl` in this directory (1 meta record + 41 technique records, C01–C46). **The inventory is partial.** The file stops at C46. The planned records for the n=83 F6 work (`f6_n83_*`), the `full_binary_norm_*` directories, the pricing-drift note referred to as C60 and the T5 proposals were not written. Records C05 and C17 still mention C60, but that record does not exist.

Repos were read-only. crypto checkout `af663ad7`; ecbench `5ca7976` with pinned `deps/crypto` `0dc9ef4a` (2026-10-14). The ecbench binary was built by the caller (sha256 `6ad57b87…`). This survey built `pdp_split_dreg` from `deps/crypto` into the scratchpad (`target-crypto/`, 10m23s). All replays ran on Linux x86-64 (Xeon @ 2.10 GHz, 4 vCPU).

## Tier counts (complete records)

| T1 | T2 | T3 | T4 | T5 |
|--:|--:|--:|--:|--:|
| 7 | 12 | 18 | 3 | 1 |

## Techniques

| id | technique | tier | verdict (short) |
|---|---|---|---|
| C01 | S₄ → S₃ chain (plan 2+1, one full-field aux u, 3ℓ+n unknowns) | T2 | solving degree 7→4, but 3.9–69× more XORs at ℓ=4; 6.6–20.8× fewer at ℓ=5, g=0 |
| C02 | two-way split m=4 (S₃ + S₄ link) | T2 | solving degree stays 7 (6 with guessing) |
| C03 | three-way split m=4 (three S₃, two aux) | T2 | solving degree 4; 113× cheaper (n=11) or 1.3× dearer (n=13) |
| C04 | BFP hybrid guessing | T2 | unsplit with g=4 is the cheapest variant in every m=3 cell |
| C05 | Frobenius-closure τ-slot base | T2 | 1.5–4.5× yield per system at the same degree; still 36×–25,257× strong rho (priced) |
| C06 | pdp-split m=2 on the invariant base (ecbench) | T3 | table oracle cheaper |
| C07 | round 3 (g=8, closure MITM bracket, ℓ-ladder) | T5 | registered, no results committed |
| C08 | F4 cost per residual unknown, n=23 | T4 | 2.25–5.99× per unknown (exploration, raw data not committed) |
| C10 | PairSumTable vs matrix-F4 (binary) | T3 | 383×/1535× per decomposition, but F4 wins single solves at n ≤ 13 |
| C11 | prime direct-subtraction root finding | T3 | 3×10³–5×10⁴× vs F4; Lean-checked identity |
| C12 | rf-mitm 2-sum table | T3 | slope at the pair-table floor; ~5.9×10⁴× F4 wall at (23,8) |
| C13 | rf-pairs (group law, memoryless) | T3 | slope at the enumeration floor |
| C14 | pairs-and-solve reference | T2 | digests reproduced; Θ(2ⁿ) whole-run cost |
| C15 | group-law root test kernel | T2 | 32–63× per target; exact |
| C16 | binary non-Koblitz MITM, ℓ=8–12 | T3 | 0 disagreements; 2^ℓ/2 fewer ops |
| C17 | rootfind-mitm/pairs in ic.pipeline | T1 | 350×–10⁴× F4 wall; IC still above rho |
| C18 | pdp3-koblitz root-finding engines | T3 | negfold 1.69× strong rho (n=13 smoke) |
| C19 | GPU root finding (sm_120) | T4 | correct in emulation; device never run |
| C20 | F4 root-shape replay | T2 | recorded null (0.95–0.98×) |
| C21 | e₃ enumeration | T3 | 553–34,901× the table |
| C22 | Macaulay redundancy ceiling (trace replay) | T3 | at most 4.3% saving |
| C23 | isogenous curves / twists | T3 | identical rank profiles |
| C24 | affine prescreen, m=2 | T2 | 5–15% fewer field ops |
| C25 | affine screen inside the m=3 chain | T1 | exact but 1.54–9.05× more work |
| C26 | E/⟨T⟩ system at m=2 ("eg") | T2 | 1.34–1.64× fewer XORs |
| C27 | parity cover s at m=3 | T3 | 1.6× faster wall, 1.9× more ops |
| C28 | Wu / Edwards coordinates | T3 | no coordinate is uniformly cheaper |
| C29 | Frobenius-power auxiliary variables | T3 | negative |
| C30 | convolution splitting of eᵢ | T3 | degree 3→2, but 1.21–2.07× more ops |
| C31 | Weil-friendly component charts | T3 | fixes ambient F4 timeouts; does not beat enumeration |
| C32 | toy Boolean elimination audit | T1 | verifier PASS |
| C33 | S₃ constant-coefficient audit | T1 | byte-identical replay |
| C34 | Gaudry S₄ specialisation engineering | T2 | local gain only |
| C35 | yield-law census vs Table 7.7 | T2 | byte-identical replay |
| C40 | native CDCL+XOR | T4 | ~1 conflict per triple; 5–62× slower than brute force |
| C41 | WDSat plain/symmetry/-x | T3 | brute force; -x incomplete |
| C42 | SMT (cvc5/Z3) | T3 | parity hypothesis refuted at n=13 |
| C43 | Hamming-weight NB base under SAT | T3 | nothing decided at n=23, m=3 |
| C44 | descent + native-XOR SAT in ecbench | T1 | S is a lower bound (conflicts not priced) |
| C45 | mitm-frobenius (incumbent table oracle) | T1 | 4.27× strong rho |
| C46 | DRAT-certified UNSAT branches | T1 | 5/5 verified |

## Contradictions and retractions (as recorded in C01–C46)

- **pdp_split_dreg, Amendment 1.** A solution read-out defect missed decompositions in cell 3:0:13:5. It is fixed in `pdp_split`; the library function is worked around, not fixed.
- **RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.** §1 runs F4 up to D=4 and gets the 383×/1535× headline. §7 runs F4 at its best max_degree (2), and F4 then wins single solves at n ≤ 13.
- **ecbench_sat_pdp pilot.** The audit failed on float GAE low bits, which led to `ic.pipeline_counted`.
- **WDSat `-x` mode.** It is incomplete: it gave a wrong UNSAT on a decomposable target. Its conflict figures are first-solution costs, not exhaustion costs.
- **SMT study.** cvc5 changed version between n=13 and n=19.
- **Thesis Table 7.7.** Its 2× is not reproduced.
- **lit2026 W8.** The repo's ×75 / 30–416× "symmetry" figure is mostly the change of curve to E/⟨T⟩.
- **pdp_split_slots.** The projection understates the measured cost (9× projected vs 36× measured).

### Pricing drift under the pinned ecbench (observed in this survey's replays)

The pinned crypto now charges F4 units whose name starts `word XORs (` at the pinned `ns_per_word_xor`. Committed sessions recorded those XORs as uncharged.

- Replays of F4 arms therefore fail the audit, but only on GAE. Every integer count is identical.
- The S direction reported for rootfind vs F4 flips. At n=13 the README gives ic-rootfind/ic-f4 = 2.56; repriced it is about 0.10.
- In the n=13 cycle-3 smoke session, the inherited-F4 word-XOR counts themselves differ: 8,952,839 recorded vs 9,603,153 replayed for seq 7. The cause is not isolated.

## ecbench commands (docs/ecbench/README.md)

- `plan --spec S`
- `run --spec S --out D [--cpus auto|inherit|none]`
- `verify --dir D [--replay N | --replay-all] [--allow-interrupted] [--out R] [--exit-code]`
- `compare --dir D --a X --b Y`
- `table --dir D…`
- `resources --dir D`
- `profile-input --spec S --seq N | exec`: re-run one measured child
- `claim build/attach/check`
- `db sql`
- `isolab-job`
- `bound fit/check`, `frontier build`, `challenge seal/spec/verdict`

## Replays run by this survey

`E=/home/user/ecbench/target/release/ecbench`. `R=/home/user/aburan28/crypto/research`.

| command | expected output (observed) | time |
|---|---|---|
| `$E verify --dir $R/ecbench_sat_pdp_20261009/sessions/count-qualified-main --replay-all` | audit OK: 240 records, 240 verified, 200 replays reproduced, 0 problems (cross-arch replay) | 28 s |
| `$E verify --dir $R/pdp_root_finder_ecbench_20261013/claim_round/sessions/n31 --replay-all` | audit OK: 128 records, 128 verified, 96 replays reproduced, 0 problems (n23 the same) | 16 s |
| `$E verify --dir $R/pdp_root_finder_ecbench_20261013/sessions/n13 --replay-all` | audit FAILED: 240/280 reproduced; the 40 failures are ic-f4 runs differing in GAE only | 12 s |
| `$E verify --dir $R/pdp_split_slots_20261012/sessions/A-n11 --replay-all` | audit FAILED: 16/48 reproduced; the 32 failures are pdp-split runs differing in GAE only | 104 s |
| `pdp_split_dreg --cells 3:0:13:4 --guesses 0,2,4 --targets 4 --seed 20261011 --budget-secs 300 --max-free 24 --lambda-targets 256 --out F --replay-summary` | runs 24, verified_total 6, word_xors_total 5680776626; profile identical to the cairn artifact | 32 s |
| same command with `--cells 3:1:11:4` | word_xors_total 1573482119, identical to runs-v2 | 12 s |
| `PYTHONDONTWRITEBYTECODE=1 python3 $R/toy_binary_formula_audit_20260924/audit.py --output F` | sha256 8b68bb62…a393, equal to results.json | 0.4 s |
| `PYTHONDONTWRITEBYTECODE=1 python3 $R/toy_binary_elimination_20260924/verify.py` | PASS | 0.4 s |
| `PYTHONDONTWRITEBYTECODE=1 python3 $R/ecc2k130_pdp_chain_holdout_20260925/verify.py …/receipt_power2048.json` | PASS | 39 s |
| `yield_census 24 4,5` (crate copied to the scratchpad) | sha256 572111a5…3338, byte-identical | 3 min |

**Determinism.** Integer counts, digests, verdicts and recovered logs are bit-identical. Wall time is never bit-identical. GAE/S depends on the pricing-table version.

Receipts and outputs are in `repro/` in this directory.
