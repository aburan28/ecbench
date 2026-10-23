# Summation-polynomial degree-of-regularity study (2026-10-23)

An **independent, twice-reproduced** measurement study of the degree-of-regularity
problem in the point-decomposition (PDP) / Semaev summation-polynomial step of
elliptic-curve index calculus, for binary-Koblitz (ECC2K-130-type) and prime-Koblitz
(`j=0` / GLV) curves. It cross-checks the `aburan28/crypto` research program with a
from-scratch harness and two external solvers the program did not author.

## Read this first

- **`REPORT.pdf`** — the full report (14 pages, figures, tables). **`REPORT.md`** is its source.

## Headline results (all independently measured here)

- **Semaev polynomials certified correct** from scratch: prime `S_3`, binary `S_3`,
  binary `S_4` all vanish exactly when the points sum to `O`, with **0 false
  positives / 0 false negatives** across 15 concrete curves (`evidence/cert_*.json`).
- **Determinism:** a 20-system canonical grid run twice is **byte-identical**
  (`evidence/canonical_run{1,2}.json`, `digest_run{1,2}.json`).
- **The Huang–Kosters–Yeo gap reproduced:** first-fall degree stays flat (2–3) while the
  solving degree grows (4→6) — the first-fall-degree assumption fails.
- **The "degree" is engine-dependent:** on the identical ideal, `msolve`-F4 and
  `Singular`-`std` disagree, and the split-vs-unsplit degree order flips with the engine
  (`evidence/engine_compare.json`). Only the coordinate-free last-fall / CM-regularity
  invariant is well defined.
- **Two "exotic" levers falsified** on toy instances: the solving degree is identical
  across **295 factor-base subspaces** (`evidence/exotic_fb_consolidated.json`) and
  invariant under monomial/variable order (`evidence/exotic_order.json`). The barrier is
  intrinsic to the ideal class.
- **The program's PKM-tower msolve degrees reproduced twice, exactly**
  (`evidence/pkm/pkm_twice.txt`): degree, F4 rounds, matrix dims all bit-identical to the
  program's recorded `results.csv`.
- **A live reproducibility defect confirmed and re-diagnosed:** the `ecbench`
  `claim-n13` replay fails (96 problems, all F4 arms; twice-run evidence in
  `evidence/REPRODUCTION_LOG.md`). v1 repeated the program's "BMI2 host path" attribution;
  v2 found the host paths already agree on `main` and traced the failure to a counted-work
  definition (lazy → eager) that changed after the session was recorded — see REPORT §4.8
  and aburan28/crypto#2459.
- **v2: all eight exotic levers run** (`exotic/rank{2..7}_*/`, REPORT §6.1): four null, one
  inconclusive, one structurally empty for prime `n`, two (both symmetry quotients) lower
  the degree — and one of those costs more in linear algebra than it saves.
- **v2: validation sweep wired into the gate** (aburan28/crypto#2451): 0 → 55 verdicts
  consumed, no row promoted while it carries caveats.
- **Bottom line:** no tried lever bends the degree-growth curve enough to threaten the
  target curves, and no end-to-end index-calculus result in the program carries a
  completed validation certificate (tier T1). This is a map, not a break.

## Layout

```
REPORT.pdf / REPORT.md     the report
harness/                   from-scratch Python harness (Semaev, Weil restriction,
                           Macaulay rank, engine drivers, figures, md->pdf)
evidence/                  raw JSON outputs, figures (svg+png), REPRODUCTION_LOG.md,
                           and the PKM reproduce logs (run twice)
surveys/                   the five technique surveys (A dreg, B symmetry, C PDP/SAT,
                           D ledger/prime, E literature + .bib); machine-readable .jsonl
```

## Reproduce

External tools used: **msolve 0.6.5**, **Singular 4.3.3.0** (both `apt install`-able),
plus Python 3 with `sympy numpy matplotlib reportlab`. Host used: 4 vCPU / 15 GB.

1. **Correctness + determinism + experiments.** The scripts in `harness/` set a module
   variable `H` at the top to the harness directory's absolute path (they are the
   exact as-run artifacts, so `H` points at the original session scratchpad). To re-run,
   edit `H` in each script to this `harness/` directory, then:
   ```
   python3 cert_prime.py ; python3 gf2.py ; python3 cert_s4.py     # certificates
   python3 canonical_run.py run1 ; python3 canonical_run.py run2   # determinism
   python3 exp_split_vs_unsplit.py ; python3 exp_engine_compare.py
   python3 exotic_fb_search.py --grid A ; python3 exotic_order.py  # exotic levers
   python3 figures.py ; python3 build_report.py                    # figures + PDF
   ```
   Deterministic outputs (degrees, hashes, ranks) must match `evidence/`; only timings vary.
2. **Reproduce the program's PKM degrees.** With the `crypto` source fetched
   (`../../tools/bootstrap_crypto.sh`) at `deps/crypto`, the systems live under
   `deps/crypto/research/ic_next10_20261014/03-pkm-degree/systems/`. Run each through
   `msolve -t 1 -v 2`; the degree / rounds / matrix dims match `evidence/pkm/*.log`
   (both runs identical). `evidence/pkm/run_pkm.sh` is the driver used (adjust its repo
   path).

## Provenance

- Program tree surveyed: `aburan28/crypto` @ `af663ad7c`.
- Engines: msolve 0.6.5-1build2, Singular 4.3.3.0, PARI/GP 2.15.4.
- Surveys A/C/D/E were independent read-only passes; survey B was reconstructed by hand
  after its agent was interrupted by a content filter. Every cross-referenced number was
  checked against the program's own files. Nothing in the `crypto` tree was modified.
