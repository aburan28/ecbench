# Shared brief for the exotic-lever experiments (ranks 2–7)

Context: `../REPORT.md` (read §1, §4.3–4.7, §6, §7) and `../surveys/survey_E_literature.md` §5
(the ranked exotic ideas, each with its cheapest falsification experiment). Ranks 1 and 8
were already run and falsified (§4.6, §4.7). You run ONE of ranks 2–7.

This is academic measurement on toy instances: we measure the degree a Gröbner engine
reaches on point-decomposition systems. No claims about deployed curves.

## Tools you may reuse (certified — do not modify them)

`../harness/` (importable: `sys.path.insert(0, <abs path of ../harness>)`):
- `gf2.py`: GF(2^n) arithmetic (`gmul`, `gsq`, `ginv`, `IRRED` for n≤13), binary curve
  `y^2+xy=x^3+a2x^2+a6` group law `bin_add`, `bin_points`, numeric `S3_bin`.
- `boolean_weil.py`: symbolic Weil restriction with exact F2 monomial algebra (field
  element = list of n sets of frozenset monomials; `f_add f_mul f_sq f_const`,
  `var_elt(prefix,i,basis,n)`, `S3_bin_sym`, `S4_unsplit`, `weil_equations`).
- `pdp.py`: `random_subspace_basis`, `subspace_elements`, `find_decomposition`,
  `macaulay_rank_seq` (F2 multilinear first-fall detector), `run_singular` +
  `parse_singular_degree` (F2 only), `run_singular_order`, `to_singular_poly`.
- `exp_chain.py`: `build_chain`, `const_elt`, `free_elt`.
- `exp_engine_compare.py`: `to_msolve`, `msolve_maxdeg` (F2 systems).
- `semaev.py`: prime-field `ec_add`, `points`, `S3_closed` (char ≠ 2,3).
For prime fields / extension fields you will write your own Singular/msolve export
(msolve: first line variables, second line characteristic, then polynomials; prime
characteristic only, `-t 1 -v 2`, degree = max of the first column of the round table
`deg sel pairs rows x cols`. Singular: `ring r=p,(vars),dp; option(prot); std(I);`
degree = max `[k:d]` in the protocol; for GF(2^k) use `ring r=(2^k,a),(vars),dp;`).

## Measurement rules (non-negotiable)

1. **Write `PROTOCOL.md` first**: hypothesis, the exact systems compared, the falsifier
   (from survey_E §5 for your rank, sharpened), and what each outcome would mean. Do
   not edit the falsifier after seeing results; add a dated amendment instead.
2. **Correctness before degree**: any new summation polynomial / encoding must be
   checked against brute-force decompositions on concrete small curves (0 false
   positives / 0 false negatives, like `../harness/cert_*.py`), and each solved system
   must recover the known planted decomposition.
3. **Report the degree from BOTH engines** (msolve F4 and Singular std) wherever the
   field allows (msolve: prime fields only). The report showed the engines disagree and
   can flip comparisons; a lever only "lowers the degree" if it does so on both, or you
   say exactly which. Also report Bézout / `vdim` and variable counts.
4. **Run every deterministic measurement twice** and diff (system hash, degrees, vdim);
   record the diff result. Timings are not compared.
5. **Toy sizes; bounded resources** — 4 cores and 15 GB are shared by ~8 concurrent
   agents. Every solver process: `nice -n 10 timeout 300 prlimit --as=2500000000 …`,
   msolve `-t 1`. Total wall budget ≈ 90 minutes. A timeout is a recorded result
   ("≥ d, not reached"), never silently dropped.
6. **Honest verdict**: lowers / neutral / raises / inconclusive, with a tier
   (T1 machine-checkable certificate, T2 independently replicated, T3 single measured
   study with raw data, T4 derivation/extrapolation, T5 not run). Null and negative
   results are results. Keep measurement and interpretation in separate sections.

## Output (write ONLY inside your own directory `exotic/<your-slug>/`)

- `PROTOCOL.md` (written before running), your code, `raw/` outputs (JSON), run1/run2.
- `RESULT.md`: verdict, a table of measured degrees (both engines) / Bézout / vars per
  system, the twice-run diff, the falsifier outcome, tier, limitations, next step.
- One figure (`fig_<slug>.png` and `.svg`, matplotlib, palette
  `#2a78d6 #eb6834 #1baf7a #eda100`, text `#0b0b0b`, legible, no dual axes).
- Do NOT modify anything outside your directory, do NOT `git commit` or push.
