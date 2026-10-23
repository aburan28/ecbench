## Abstract

The point-decomposition problem (PDP) at the heart of elliptic-curve index calculus
asks a Gröbner/XL engine to solve a Weil-descended Semaev summation-polynomial system.
Its cost is governed by the degree the solver must reach, and the central question for
both binary (Koblitz / ECC2K-130) and prime (j=0 "prime-Koblitz" / GLV) curves is
whether that degree stays bounded or blows up. This report (a) catalogues every
degree-control technique tried in the `aburan28/crypto` research program with an
explicit confidence tier, (b) rebuilds a **fully independent** measurement harness —
Semaev polynomials constructed from scratch and certified against real point
additions, solved by two external engines (`msolve` 0.6.5, `Singular` 4.3.3.0) the
program did not author — and (c) establishes reproducibility by running every headline
number **twice** and comparing. The independent harness reproduces the program's
central phenomenon (first-fall degree stays flat while the solving degree grows),
reproduces its strongest independently-replicable prime-field result (PKM-tower msolve
degrees, bit-identical on re-run), independently confirms a live reproducibility
**defect** on `main` (a BMI2 counter path that makes a recorded claim fail replay),
and tests two never-before-run "exotic" levers on toy instances — both falsified with
clean negative results. The honest bottom line is unchanged and now triple-checked:
**no technique bends the degree-growth curve enough to threaten the target curves, and
no end-to-end index-calculus result in the program carries a completed validation
certificate (tier T1) — only isolated sub-facts do.** The value delivered is not a
break; it is a reproducible, cross-validated map of what has been tried, what actually
happened, and how much to trust each claim.

> **Scope and honesty note.** This is a defensive cryptanalysis / benchmarking study
> on toy and small instances. Every "break"-shaped quantity is reported as a distance
> to Pollard rho, never as a threat to a deployed curve. Interpretation is labelled as
> such and kept separate from measurement, per the repository's own reporting rules.

---

## 1. The problem and the three degrees

For an elliptic curve `E/F_q` the `m`-th **Semaev summation polynomial** `S_m(X_1,...,X_m)`
vanishes exactly when there exist points `P_i=(X_i,*)` on `E` with `P_1+...+P_m=O`
(up to sign). Index calculus (Gaudry–Diem) fixes a **factor base** — for `GF(2^n)`, the
points whose `x`-coordinate lies in an `F_2`-subspace `V` of dimension `ell` — and
decomposes a random point `R` as a sum of `m` factor-base points by solving
`S_{m+1}(x_1,...,x_m,x(R))=0` with each `x_i in V`. Writing each `x_i` in a basis of `V`
and expanding the single `GF(2^n)` equation over an `F_2`-basis (Weil restriction) turns
this into a Boolean polynomial system in `m·ell` variables, which an F4/F5/XL/SAT engine
solves. **The entire difficulty of elliptic index calculus lives in the degree that
solve reaches.**

Three invariants are routinely conflated and must be separated (Caminata–Gorla 2023;
Huang–Kosters–Yeo 2015):

- **first-fall degree** `d_ff`: the lowest degree at which a non-trivial degree fall
  (syzygy) appears. For these descent systems `d_ff = m^2-m+1` (Kousidis–Wiemers 2015),
  small and roughly constant in `n`.
- **degree of regularity** `d_reg`: for a semi-regular sequence, the index of the first
  non-positive coefficient of the truncated Hilbert series. The Bardet–Faugère–Salvy
  prediction grows with the variable count.
- **solving / last-fall degree** `d_solv`: the degree actually needed to finish. The
  last-fall degree (HKY) is monomial-order-independent; `d_solv <= reg(I^{hom})`, the
  Castelnuovo–Mumford regularity of the homogenized ideal (Caminata–Gorla).

The binary subexponential claims (FPPR 2012, Petit–Quisquater 2012, Semaev 2015) **all
assume `d_reg ≈ d_ff`**. Whether that holds is the pivotal dispute, and it is the first
thing the independent harness was built to test.

---

## 2. Independent measurement harness

Trusting the program's own engines would defeat the purpose, so the harness shares no
code with it. It is ~700 lines of Python (exact `F_2` monomial algebra, so the field
equation `c^2=c` is built in) plus two external solvers. Everything is in
`harness/` with raw outputs in `harness/out/`.

**Certified-correct Semaev polynomials.** `S_3` (prime and binary) and `S_4` (binary,
built as the resultant of two `S_3` quadratics) are each checked to vanish *exactly*
when a sign choice of the points sums to `O`, over concrete curves — a real correctness
certificate, not a claim:

| polynomial | curves checked | false positives | false negatives |
| --- | --- | --- | --- |
| prime `S_3` | p = 13,17,19,23,29 (5 curves) | 0 | 0 |
| binary `S_3` | n = 3,4,5, incl. Koblitz a2∈{0,1},a6=1 (6 curves) | 0 | 0 |
| binary `S_4` (resultant) | n = 4,5,6 (4 curves) | 0 | 0 |

**External engines, two jobs.** `Singular`'s `std` with protocol gives the max step
degree (solving degree); `msolve`'s verbose F4 gives the max round degree; both confirm
zero-dimensionality (`dim=0`) and the solution count (`vdim`). Using *two* engines is
not redundancy — §4.5 shows they disagree, and that disagreement is itself a result.

**Determinism protocol.** A fixed 20-system grid is run twice; the emitted file records
only deterministic invariants (system SHA-256, `d_ff`, solving degree, `vdim`, `dim`,
the full Macaulay rank profile) — timings are excluded by construction. The two runs
must be byte-identical.

Host for every number below: 4 vCPU, 15 GB, Linux 6.18.44. `msolve` 0.6.5, `Singular`
4.3.3.0, `PARI/GP` 2.15.4. The program's own binary used for replays is `ecbench`
pinned to crypto `0dc9ef4`; the program tree surveyed is at `af663ad7c`.

---

## 3. What actually controls the degree — catalogue with confidence tiers

Tiers (applied uniformly): **T1** machine-checkable certificate / validation record;
**T2** independently replicated with matching outputs; **T3** single measured run with
raw data in the repo; **T4** prose, derivation or extrapolation; **T5** designed but
never run.

### 3.1 Levers that demonstrably move the degree

| lever | mechanism | curve family | measured effect | verdict | tier |
| --- | --- | --- | --- | --- | --- |
| Field equations `x^q=x` | force zero-dim / multilinear quotient; cap degree | all | the reason any of these terminate at low degree | essential, not optional | T1 (structural) |
| Fall reuse (MutantXL / F4) | re-multiply discovered degree falls | binary chain | refutes at **3–4** where plain XL needs **5–7**; holds to N=54 (m=4,ℓ=6) | **lowers** (vs XL) | T2 |
| Torsion symmetrisation (FGHR/FHJRV) | rewrite `S_m` in `S_m`- and `[±1]`-invariants | binary/prime | `S_4` per-variable degree `[4,4,4,4]→[2,2,2,2,1]` | lowers degree & Bézout | T2 |
| `D_4=(Z/2)^3⋊S_4` weighted `S_5` | solve in 4 weighted invariants | `F_{p^4}`, m=4 | Bézout 4096→512 (derived; quotient-dim 512 **measured**); matrix 40920→4200 cols; **48×** on C4 (not the 630× projected) | lowers (engineering constant) | T3 |
| Semaev `S_3`-chain split | `S_{m+1}` → m−1 copies of `S_3` + aux vars | binary/prime | keeps per-link degree 2; see §4.3 for the cost | **mixed** (raises degree at small m) | T2 |
| Variable fixing at ℓ=7 | guess one coord per summand | binary chain | brings the saturating ℓ=7 closure back to degree 4, at 2^k cost | lowers per branch | T3 |
| `rr` symmetric-norm presentation | norm-form coordinates vs direct `S_4` | Koblitz n=23–61 | refutes 3–4 degrees below direct descent (6 vs 10 at K1/2^61,ℓ=6) | **lowers** | T3 |
| F4 early stop (Hilbert-driven) | terminate at the last productive degree | prime JV cover | removes last degree step (8→7); ~4–5× fewer F4 muls | lowers (termination) | T3 |

### 3.2 Levers that did NOT bend the curve (honest negatives)

| lever | why it fails / net effect | tier |
| --- | --- | --- |
| Choice of factor-base subspace `V` | solving degree identical across 295 subspaces (§4.6) — intrinsic, not a knob | T2 (this report) |
| Monomial / variable order | solving degree invariant under grevlex/grlex + all permutations (§4.7) — last-fall degree is order-independent | T2 (this report) |
| Frobenius-power auxiliary variables | `x↦x^{2^k}` is `F_2`-linear → changes nothing (repo W8 §3) | T3 |
| Convolution `e_i` splitting (La Scala) | 1.2–2× *more* work | T3 |
| Wu `Z/4` state model | neutral | T3 |
| Isogeny / curve move within a class | degree invariant across 2,994 comparisons | T3 |
| Trace-zero / Artin–Schreier pre-filter | never fires on odd-order IC targets (0 / 29,624) | T3 |
| Hamming-weight ideals | bounded generator degree but call exponent 1.9–2.5 | T3 |
| Trilinear XL (unknown `t`) | cost rises (87.16 > 71.82 bits at n=131) | T3 |
| Crossbred / monomial dropping | margin shrinks with size; no stable `(D,k)` | T3 |
| Frobenius-invariant `V` for prime `n` | **structurally empty**: cyclic-code dims are cyclotomic-coset subset sums; at n=131 only {0,1,130,131} exist | T1 (proof) |

### 3.3 Designed but never run (T5)

`D_4`-symmetrised `S_5` **lattice** oracle; **Magma** F4/F5 on the PKM systems (so the
"beat the fastest external tool" rule is formally unmet for the summation systems);
order-6 automorphism quotient on `j=0` prime curves (§7); theta/Kummer quadratic
encoding; Lattès / isogeny-kernel factor base; non-abelian wreath block-decomposition.

![Confidence-tier distribution of the program's end-to-end index-calculus evidence (survey D's machine-backed leaderboard tally). Zero T1 at this level: no end-to-end IC result carries a gate-complete VALIDATION.json. Component-level machine-checkable certificates do exist separately (survey C counts 7: DRAT-certified UNSAT branches, Boolean-elimination verifier passes, byte-identical constant-coefficient replays) — but none certifies an end-to-end relation→solve→DLP claim.](out/fig/fig_tiers.png)

[PAGEBREAK]

## 4. Independent measured results

### 4.1 The first-fall vs solving-degree gap (the pivotal dispute), reproduced

On the split `S_3`-chain for `m=3`, the harness measures `d_ff` pinned at **2–3** while
the Singular solving degree climbs **4 → 5 → 5 → 6 → 6** as the system grows. This is a
clean, independent reproduction of the Huang–Kosters–Yeo phenomenon: **`d_reg` grows
while `d_ff` stays flat**, so the first-fall-degree assumption underlying the binary
subexponential claims does not hold for these systems.

![First-fall degree stays flat while the solving degree grows (independent harness, Singular F4). This is the Huang–Kosters–Yeo divergence, reproduced from scratch.](out/fig/fig_ff_gap.png)

### 4.2 Determinism: run once, then twice

The 20-system canonical grid was executed twice. The two output files are **byte-identical**
(0 mismatches across system hashes, first-fall degrees, solving degrees, `vdim`, `dim`,
and full Macaulay profiles). Every deterministic invariant this report rests on is
reproducible by construction; only wall-clock timings vary.

### 4.3 Splitting a point: Semaev's chain vs the unsplit resultant — and why the engine decides

Semaev's 2015 splitting trick replaces one big `S_{m+1}` by a chain of `S_3`'s to keep
each equation quadratic. Whether that *helps* turns out to depend entirely on the
engine, and getting this wrong is an easy way to publish a backwards conclusion.

Under **`msolve`'s F4** (the fall-reusing engine index calculus actually uses),
splitting **lowers** the degree: the unsplit `S_4` solves at degree 6–7 and the split
chain at 3–4. This is Semaev's intended benefit, and survey C's independent F4
measurement reports exactly it — solving degree **7 → 4** on every cell. Under
**`Singular`'s `std`**, the direction *reverses*: the unsplit `S_4` (few variables,
high degree) solves at degree 3 while the split chain's extra `n` auxiliary Boolean
unknowns push it to 5–6.

![The same ideal gives a different "solving degree" under different engines, and the split's effect flips sign: msolve-F4 (fall reuse) goes 6→4 as you split (Semaev's intended win), while Singular std goes 3→5. Any bare "d_reg = k" is under-specified without naming the engine.](out/fig/fig_split_engines.png)

So the honest statement is **not** "splitting is better" or "worse" — it is that the
split is a genuine degree *trade* whose sign depends on the solver, and under the
engine that matters (F4) it does what Semaev intended. The cost in field operations is a
third, separate axis: survey C measures the split at 3.9–69× *more* XORs at `ℓ=4` but
6.6–20.8× *fewer* at `ℓ=5`, and finds that **plain unsplit with guessing (`g=4`) beats
every split variant in every `m=3` cell** — no variant came within 22.8× of strong rho.
A three-way split at `m=4` reaches degree 4 (vs 7 for two-way) but is 113× cheaper at
`n=11` and 1.3× dearer at `n=13`. The degree, the XOR cost, and the distance to rho move
independently; a lever can help one and hurt another.

### 4.4 The semi-regular prediction massively overstates reality

For the prime-field PKM tower the Bardet–Faugère–Salvy degree of regularity is `N+1`
(top-part Hilbert series `(1+t)^N`) and XL gives `N+4`. The program's measured F4
solving degrees — four of which are independently reproduced here, §4.8 — are far lower:
`5`–`7` across `N=8..18` (e.g. `5,6,6,7,7` for `m=3` at `N=6,9,12,15,18`), growing
roughly **one degree per eight tower unknowns**. F4 escapes the XL bound through degree
falls. (The program's `m=4,N=20` cell reaching `≥8` refutes the "plateau at 7"
hypothesis; that cell is beyond this host's memory.)

### 4.5 The measured "degree" is engine- and encoding-dependent

Running `msolve`-F4 and `Singular`-`std` on the **identical** ideal gives different
degrees, and the ordering flips with the encoding:

| system | variables | Singular `std` | msolve F4 |
| --- | --- | --- | --- |
| split chain, n=8, ℓ=2 | 14 | 5 | 4 |
| split chain, n=10, ℓ=2 | 16 | 5 | 3 |
| split chain, n=8, ℓ=3 | 17 | 6 | 4 |
| unsplit `S_4`, n=8, ℓ=2 | 6 | 3 | 6 |
| unsplit `S_4`, n=10, ℓ=2 | 6 | 3 | 6 |
| unsplit `S_4`, n=9, ℓ=3 | 9 | 4 | 7 |

On the chain, `msolve`'s fall-reusing F4 reaches **3–4** where Singular reaches 5–6 —
exactly reproducing the program's "fall reuse keeps the refutation degree low" finding.
On the unsplit system the order *reverses*. **Conclusion:** the headline "degree of
regularity" is not a single well-defined number; it depends on the engine's fall-reuse
strategy and on the encoding. Only the coordinate-free last-fall / CM-regularity
invariant is intrinsic — and that is precisely the quantity that resists cheap
measurement. Any claim stated as a bare "`d_reg = k`" is under-specified.

### 4.6 Exotic lever #1 — choosing the factor-base subspace (FALSIFIED)

The cheapest never-tried idea in the literature (survey E, rank 1): treat `V` as a free
design variable and pick it to minimise the solving degree (Hilbert-series / CM-regularity
targeting; the Caminata–Gorla frame is cited nowhere in the program). Across **295
subspaces** (random + Frobenius-stable) in five `(model, n, ℓ)` configurations, the
solving degree is `min = median = max` *every time*; Frobenius-stable `V` give the
identical value.

![Exotic lever #1 falsified: across 295 subspaces the solving degree has zero spread at every configuration. V is not a design knob — the degree is intrinsic to the ideal class, as Caminata–Gorla CM-regularity predicts.](out/fig/fig_fb_distribution.png)

**Verdict: falsified, and informatively so.** The degree-growth is a property of the
ideal class, not of the subspace, so `V`-shaping cannot soften the barrier. This is the
"publishable negative" the frame predicted.

### 4.7 Exotic lever #2 — monomial/variable ordering (FALSIFIED as a lever)

Under grevlex, grlex, and all four random variable permutations, the solving degree is
invariant (3 unsplit / 5 split). Only pure lex reaches a different number, and lex
"degree" is not the solving degree (FGLM-style elimination is a different regime). This
confirms the last-fall degree's order-independence (HKY) and closes ordering as a lever.

### 4.8 Reproducing the program's own results — twice

**PKM-tower msolve degrees (the program's one genuine T2 algebraic-summation result).**
Each of four systems was solved **twice** with `msolve`; the degree, F4 round count, and
max matrix dimensions match the recorded `results.csv` bit-for-bit, both runs:

| system | recorded (D / rounds / rows×cols) | run 1 | run 2 |
| --- | --- | --- | --- |
| m2,t6 (N=12) | 5 / 14 / 2803×2064 | identical | identical |
| m3,t3 (N=9) | 6 / 17 / 3827×2997 | identical | identical |
| m4,t2 (N=8) | 6 / 22 / 2009×1906 | identical | identical |
| m3,t4 (N=12) | 6 / 19 / 12900×12287 | identical | identical |

**The replay that FAILS — and why.** The `ecbench` replay of the root-finder
`claim-n13` session fails on this host: 224 records verified but **96 problems**, all in
the IC (F4) arms; all 72 rho rows match. Run twice, the failure is itself deterministic
(identical 96 problems). The cause, isolated in the program's own sweep and confirmed
here, is a counter bug: the BMI2 `pext`/`pdep` path in `inherited_f4.rs:942` charges
different word-XOR counts than the portable path, so a claim recorded on Apple arm64
does not replay on an x86 host. **This is a real, open reproducibility defect on `main`**
— exactly the "run it twice, get different numbers" failure mode, caught and explained.
By contrast the clean PKM-tower ecbench session (`p20`) replays OK both runs (192/192,
0 problems).

> Reproducibility verdict: the summation-step *science* is reproducible (every degree
> here is bit-identical on re-run), but the program's *cost accounting* has a live
> host-dependent defect that invalidates at least one recorded claim's replay. Both
> facts are now independently established.

### 4.9 Alternative solvers: root-finding and SAT do not change the verdict

Two families sidestep the Gröbner degree entirely, and neither crosses rho (survey C):

- **Root-finding oracles** (pair tables, rf-mitm, group-law root tests) beat F4 in wall
  time by `10^2.5`–`10^4×` per decomposition and are correctness-exact (one identity is
  even Lean-checked). But they sit at the pair-table enumeration floor `S·T^2 ≳ r`, so
  the resulting index calculus stays above rho (the best binary row is `1.69×` strong rho
  at a smoke size; the prime incumbent `4.27×`). A caveat survey C flags: the headline
  "`383×`/`1535×` root-finding beats F4" uses F4 capped at degree 4; run at F4's best
  setting, F4 is *faster* for whole solves at `n ≤ 13`. Apples-to-apples matters.
- **SAT / SMT / WDSat** behave as **brute-force enumerators** — about one conflict per
  candidate tuple, no algebraic cancellation captured — and are 5–62× slower than direct
  enumeration. A correctness landmine worth recording: **WDSat's `-x` mode is incomplete**
  and returned a *wrong UNSAT* on a decomposable target. A solver that can call a true
  instance unsatisfiable is a reminder that "no relation found" is not "no relation
  exists" without a certificate (which is exactly what the DRAT-certified UNSAT branches,
  §5, provide and what WDSat `-x` does not).

[PAGEBREAK]

## 5. Confidence and reproducibility — the "no trust-me-bro" audit

- **Zero end-to-end T1, but component certificates exist.** The validation gate exists
  and passes, but **no real `VALIDATION.json` exists** — every one of the 28 leaderboard
  rows and 116 benchmark rows renders as "unvalidated (pre-gate)"; the only passing case
  is a test fixture. So no *end-to-end* IC result (relation collection → linear algebra →
  DLP) is certificate-grade. At the *component* level, survey C counts 7 machine-checkable
  certificates (DRAT-certified UNSAT branches, a Boolean-elimination verifier PASS,
  byte-identical `S_3` constant-coefficient replays); these certify sub-facts, not a
  break. The distinction matters: "we proved this UNSAT branch" is T1; "this pipeline
  beats rho" is not certified at all.
- **The validation sweep is disconnected from the gate.** Sweep A produced 55 verdicts
  (1 validated, 45 with caveats, 8 unvalidated, 1 contradicted), but `load_verdicts()`
  reads the key `row` while the sweep writes `row_id`, and the verdict words are in
  neither pass nor fail set — so the gate consumes **0** of them. Real work, not wired in.
- **Only three rows** genuinely exercise an algebraic summation-polynomial solve at T2:
  the n9 F6-IC control, the PKM-tower degrees (reproduced here), and the determinism
  finding itself. Every larger Boolean F4/F5/F6 arm is T3 because none replays on
  current binaries (the BMI2 bug).
- **This report's own additions** sit at the top of the tier scale for *what they
  certify*: the Semaev correctness certificates are machine-checkable (T1-grade for
  correctness), the determinism certificate is byte-identical twice (T1-grade for
  reproducibility of the pipeline), and the PKM reproduction is an independent replay
  matching the repo (T2). None of these is a rho-beating claim, and none is presented
  as one.
- **The degree/science results DO reproduce.** Independently of this report's harness,
  survey A re-ran six of the program's dreg results at HEAD and all six matched the
  committed records value-for-value — including the MutantXL closure degree histogram
  ({3:16, 4:76} over 92 draws), the sum of closure products (4,902,150), and the
  `rr`-vs-`S_4` refutation degrees (5 vs 8). The algebraic content is solid and
  reproducible; it is only the *priced cost* layer that drifts.
- **Pricing drift, from a second angle.** Survey C independently hit the same class of
  failure: the pinned `ecbench` now charges F4 word-XORs at a fixed calibration rate while
  committed sessions left them uncharged, so F4-arm replays fail the audit *even though
  every integer count is identical* — and it reverses a headline (at `n=13`, recorded
  root-finding `2.56×` F4 in the `S` unit becomes `≈0.10×` when repriced). Survey C also
  found one **unexplained** F4 XOR-count mismatch it could not localise (8,952,839 recorded
  vs 9,603,153 replayed). So two distinct accounting defects exist, one understood (BMI2
  path) and one not — a strong argument that no cost-unit claim should be trusted until
  the pricing and counting are deterministic across hosts.
- **Determinism map** (what is safe to compare across runs): refutation/solving degree,
  Macaulay ranks, system hashes, `vdim`/`dim`, msolve round counts and matrix dims — all
  bit-identical. Word-XOR counters are host-path-specific (the bug). Wall-clock, RSS,
  CPU% are nondeterministic and never compared.

---

## 6. Literature placement and the exotic frontier

Three facts from the literature constrain what *can* work (survey E):

1. **Solving degree is invariant under linear / Möbius / `F_2`-Frobenius coordinate
   changes** (Caminata–Gorla; FHJRV Lemma 6: every degree-2 coordinate is a homography).
   So no reparametrisation lever — tropical ordering, Frobenius/linearised auxiliaries,
   Edwards/Hessian frames — can lower the intrinsic degree. **Confirmed here** by §4.6–4.7.
2. **`d_reg ≈ d_ff` is false/doubtful**; `d_reg` grows while `d_ff` stays constant
   (HKY; Kosters–Yeo; and §4.1 here). A lever must bend the *growth in ℓ*, and measuring
   `d_ff` alone proves nothing.
3. **Symmetry/isogeny quotients are the only proven large lever, and even maxed they do
   not beat rho** for the target curves; for prime `n` the Frobenius-subspace and
   intermediate-field levers are structurally empty.

Against this frame, the eight ranked exotic ideas and their status:

| # | exotic lever | cheapest test | status here |
| --- | --- | --- | --- |
| 1 | `d_reg`-minimising subspace `V` | 295-subspace sweep | **run → falsified** (§4.6) |
| 8 | tropical / order-based | multi-order sweep | **run → falsified** (§4.7) |
| 2 | order-6 automorphism quotient on `j=0` | S_3 quotient on ζ-stable base | designed (§7); T5 |
| 4 | Artin–Schreier trace elimination | residual-degree compare | repo: 0 fires on IC targets; T3 negative |
| 3 | wreath / full-reflection equivariant | invariant-ring solve over `F_{p^2}` | obstruction found: `Z/2` part is modular in char 2; T5 |
| 5 | theta / Kummer quadratic encoding | S_4 via theta vs chain | T5 (strong null prior) |
| 6 | Lattès / isogeny-kernel factor base | degree compare over `F_p` | T5 |
| 7 | intermediate-field τ-orbit descent | composite-n only | structurally empty for prime n (131/233/283) |

The two cheapest and most decisive (1 and 8) were run and **both falsified** — a
genuine, reproducible contribution: the two most plausible "shape the problem to tame
`d_reg`" ideas do not work, for a reason the theory (CM-regularity invariance) predicts.

---

## 7. The named targets: binary Koblitz and prime-Koblitz (j=0)

**Binary Koblitz / ECC2K-130 (prime `n`).** Every degree lever above either is
structurally empty (Frobenius subspaces: cyclotomic-coset argument), or lowers a
constant without bending the `d_reg(ℓ)` growth, or (fall reuse) holds degree 3–4 only
until it saturates at `ℓ=7` where the next degree is unmeasured for want of ~200 GB.
The program's extrapolated ECC2K-130 cost rows span `2^68…2^197` and none drops below
rho; the free-oracle floor `S·T^2 ≳ r` binds independently of the degree (T4).

**Prime-Koblitz (`y^2=x^3+b`, `p≡1 mod 3`, j=0, order-6 automorphism).** The order-6
automorphism `ζ:(x,y)↦(ωx,-y)` is currently used only as a GLV `÷6` fold on the factor
base (`glv-orbit`), *not* as a symmetry quotient of the PDP ideal. The PKM-tower route
(reproduced in §4.8) shows the prime-field solving degree is low (5–7 to N=18) but the
cost is `29–141×` rho and the first arity that could win is `m=4`, needing an `S_5`
small-root oracle at reach `δ≈1/5` — against a measured reach of `0.014` and a proven
ceiling of `1/32`. **The concrete open experiment** (survey E rank 2) is to build `S_3`
in `⟨ζ⟩`-invariants (power sums of `x_i^6`) on a `ζ`-stable factor base and measure the
quotient's solving degree and Bézout against the plain system — the `D_4` recipe applied
to the extra prime-field automorphism. By the `D_4` precedent a `Z/6` grading should cut
Bézout by up to 6; whether it lowers the *degree* is the unanswered question. It is
cheap (`p≤2^16`, `m=2`) and specified; it was not run here.

---

## 8. Open problems, what would unblock them, and the bottom line

**Open problems carried forward** (from the program's own records, with unblock criteria):

- `binary-boolean-system-64-variable-cap`: the system builder caps at 64 unknowns;
  128 unknowns agreeing with the floor on 12 targets at `m=4,n=23,d=6` would unblock
  the `ℓ=7` degree-5 measurement (currently out of memory at 169–405 GB).
- The **PKM / descent memory wall**: `m=3,N=21` needs a single ~192 GB process; **no
  problem record exists** for it (an AGENTS.md gap).
- The **BMI2 replay regression** (`inherited_f4.rs:942`): no problem record exists;
  tracked only in a sweep. Until fixed, recorded F4 word-cost claims do not replay
  across ISAs. **A one-line BMI2-off reproduces the reference counts** (168/168), which
  is the obvious interim fix.
- **Magma never run** on the summation systems (T5): the "beat the fastest external
  tool" gate is formally unmet; an msolve-vs-Magma parity run on the PKM systems would
  close it.

**The bottom line, now triple-checked.** The degree of regularity of the Semaev PDP
systems is not tamed by any tried lever: factor-base shaping and ordering are falsified
(this report), symmetry/automorphism quotients lower constants and Bézout but not the
`d_reg(ℓ)` growth, fall-reuse holds a low *refutation* degree only up to a saturation
point, and the semi-regular prediction overstates everything. For the named curves —
binary Koblitz with prime `n`, and prime-Koblitz `j=0` — the structural levers are empty
or uncompetitive, and **no end-to-end result is certificate-grade** (only component
sub-facts are). What this study adds is not a
break but the thing the program most needed: an **independent, twice-run, cross-engine
measurement** confirming where the barrier is, why it is intrinsic (CM-regularity
invariance), and which of the field's open claims actually reproduce.

The single most valuable next action is not another lever but **wiring the validation
sweep into the gate and fixing the BMI2 counter** — because an index calculus you cannot
replay cannot be improved, and right now the program's own best results do not replay.

---

## Appendix A — reproduce everything

All commands assume the harness directory. External tools: `msolve` 0.6.5, `Singular`
4.3.3.0. Deterministic outputs are noted; only timings vary.

```
# correctness certificates (prime + binary Semaev, ~2 min)
python3 cert_prime.py        # -> out/cert_prime.json   (0 FP / 0 FN)
python3 gf2.py               # -> out/cert_binary.json  (0 FP / 0 FN)
python3 cert_s4.py           # unsplit S4 certificate   (0 FP / 0 FN)

# determinism: run the canonical grid twice, diff (must be identical)
python3 canonical_run.py run1 ; python3 canonical_run.py run2
diff <(jq -S . out/canonical_run1.json) <(jq -S . out/canonical_run2.json)

# split vs unsplit, first-fall gap, engine comparison, V/order independence
python3 exp_split_vs_unsplit.py
python3 exp_engine_compare.py
python3 exp_fb_search.py --grid A        # subspace independence
python3 exp_order.py                     # order independence

# reproduce the program's PKM-tower msolve degrees (twice)
bash out/pkm/run_pkm.sh                   # all four match results.csv, both runs

# figures
python3 figures.py
```

Raw evidence: `out/cert_*.json`, `out/canonical_run{1,2}.json`,
`out/split_vs_unsplit.json`, `out/engine_compare.json`,
`out/exotic_fb_consolidated.json`, `out/exotic_order.json`,
`out/REPRODUCTION_LOG.md`, and the survey files `survey_{A,B,D,E}*.{md,jsonl}`.

## Appendix B — provenance

- Program tree: `aburan28/crypto` @ `af663ad7c`; CLI: `aburan28/ecbench` (crypto pin `0dc9ef4`).
- Host: 4 vCPU, 15 GB RAM, Linux 6.18.44. `msolve` 0.6.5-1build2, `Singular` 4.3.3.0, `PARI/GP` 2.15.4, Python 3 with `sympy`, `numpy`, `matplotlib`, `reportlab`.
- Figures use the validated accessibility palette (light mode); every number traces to a file in `out/`.
- Surveys A/C/D/E were compiled by independent read-only passes over the program tree (survey C's inventory is partial — it stops at record C46 — but covers the split / root-finding / SAT slice used here); survey B was reconstructed directly after its agent was interrupted by a content filter. Every cross-referenced survey number was checked against the program's own files.
