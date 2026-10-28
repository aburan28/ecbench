# Survey E — Index calculus for ECDLP via Semaev summation polynomials: controlling the degree of regularity / solving degree of the point-decomposition systems

Compiled 2026-10-22 for the research report. Companion BibTeX: `survey_E_literature.bib` (same directory).

**Curves of interest.** (a) binary `GF(2^n)`, especially Koblitz / anomalous binary curves with Frobenius `tau` (ECC2K-130 `n=131`, sect233k1 `n=233`, sect283k1 `n=283`); (b) prime-field curves with extra automorphisms: `j=0`/GLV ("prime Koblitz", `y^2=x^3+b`, `p≡1 mod 3`, order-6 automorphism) and `j=1728` (`y^2=x^3+ax`, order-4).

**The quantity in focus.** For the PDP one Weil-descends the `(m+1)`-th summation polynomial `S_{m+1}` over a factor-base subspace/variety and solves a zero-dimensional polynomial system by F4/F5 (or SAT/XL). The cost is governed by the highest Macaulay degree the solver must reach. Three distinct invariants are conflated in the literature and must be kept apart:
- **first fall degree** `d_ff` — the lowest degree at which a non-trivial degree fall (syzygy) appears;
- **degree of regularity** `d_reg` — (semi-regular) the index of the first non-positive Hilbert-series coefficient;
- **solving degree** `d_solv` / **last fall degree** `d_lf` — the degree actually needed to finish; `d_lf` is monomial-order-independent (Huang–Kosters–Yeo) and `d_solv ≤` CM-regularity of the homogenized ideal (Caminata–Gorla).
The subexponential binary claims all assume `d_reg ≈ d_ff`; the dispute (§3) is whether that holds.

---

## 1. Chronology table

| year | work | venue / id | setting | the dreg-relevant claim |
|---|---|---|---|---|
| 2004 | Semaev | ePrint 2004/031 | any `E` | Defines `S_m` (symmetric, degree `2^{m-2}` per variable, `m≥3`). PDP = bounded roots of `S_{m+1}`. Conjectures low solving cost; gives **no** zero-finding algorithm. |
| 2006/09 | Gaudry; Diem | ANTS / JSC 44(12) | `E(F_{q^n})`, `n` fixed | Index calculus via `S_{m+1}` + Weil descent; Gaudry heuristic `Õ(q^{2-2/n})`. |
| 2011 | Diem | Math. Comp. 80(273) | `E(F_{q^n})`, `n≥2` fixed | Rigorous `Õ(q^{2-2/n})`. No dreg heuristic needed (the descent degree is bounded by the fixed `n`). |
| 2011 | Diem | Compositio 147(1) | `n→∞`, `n/log q→0` | Subexponential `(q^n)^{o(1)}`. The strongest **proven** result; does not cover prime-degree binary. |
| 2012 | Faugère–Perret–Petit–Renault (FPPR) | EUROCRYPT 2012 / 2012/146 | `F_{2^n}` | Reformulates PDP as a Boolean system in a vector subspace `V⊂F_{2^n}`; conjectures `d_reg` barely above the equation degree → heuristic `2^{c n^{2/3} log n}`. **First "first fall degree" style analysis.** |
| 2012 | Petit–Quisquater | ASIACRYPT 2012 | `F_{2^n}` | Under the first-fall-degree assumption, Diem-style IC is subexponential `O(2^{c n^{2/3} log n})`, `c<2`; crossover with generic methods only near `n≈2000`. |
| 2013 | Huang–Petit–Shinohara–Takagi (HPST) | IWSEC 2013 | `F_{2^n}` | Symmetrized variables cut time+memory; report `d_ff` assumption **holds in their experiments**, deduce improved complexity. |
| 2013 | Shantz–Teske | ePrint 2013/596 | `F_{2^n}` | Experimental test of Petit–Quisquater; confirms the whole chain "rests on unproven heuristic assumptions on Gröbner performance." |
| 2013 | Ding–Schmidt | LNCS 8260 | general | Separates solving degree, `d_reg`, first fall degree over `F_q`; warns `d_ff` can be forced artificially low by a planted syzygy. |
| 2014 | Faugère–Gaudry–Huot–Renault (FGHR) | J. Cryptology 27(4) | `F_{q^n}`, Edwards | Symmetric-group **and** torsion invariants → "exponential speedup" in PDP. The symmetry quotient lowers Bézout/Macaulay degree. |
| 2014 | Faugère–Huot–Joux–Renault–Vitse (FHJRV) | EUROCRYPT 2014 / hal-00935050 | `F_{q^n}` | 2-torsion symmetrization; computes up to `S_8` symmetrized. **Lemma 6: every degree-2 coordinate is a homography** `φ(P+T)=f_T(φ(P))` (Prop. 8 normal forms `t↦-t` odd char, `t↦t+1` char 2). Full `E[2]` reflection group `(Z/2)^{n-1}⋊S_n`, factor base ÷8. |
| 2014 | Galbraith–Gebregiyorgis | INDOCRYPT 2014 / 2014/806 | `F_{2^n}`, `n` prime | Binary-Edwards variables invariant under a large group **lower the summation-polynomial degree**; symmetry-breaking factor base; SAT vs Gröbner. Verdict: **rho still much faster** at reasonable sizes. |
| 2015 | Semaev | ePrint 2015/310 / arXiv 1504.01175 | `F_{2^n}` | Splits `S_m` via resultants into a **chain of `S_3`'s** (fewer vars, degree stays low); heuristic `2^{c√(n ln n)}`, `c≈1.69`, under a first-fall assumption (`d_reg ≤ 4` claimed for the cubic Boolean system). |
| 2015 | Karabina | ePrint 2015/319 / arXiv 1504.02347 | `F_{2^n}` | Independent "split the point" `S_3`-chain PDP encoding. |
| 2015 | **Huang–Kosters–Yeo (HKY)** | **CRYPTO 2015 / 2015/573** | `F_{2^n}` | Introduces the **last fall degree** (order-independent). **Computational evidence casting doubt on the first-fall-degree assumption**; constructs a summation-polynomial Weil-descent system where it is "unlikely to hold." |
| 2015 | Kosters–Yeo | arXiv 1503.08001 | `F_{2^n}` | `d_ff` of `S_3` systems; the `m=2` case "somewhat pathological"; evidence `d_reg` **grows with `n`**. |
| 2015 | Kousidis–Wiemers | ePrint 2015/1121 | `F_{2^n}` | `d_ff` of the `S_{m+1}` descent bounded by `m^2−m+1` (improves earlier `m^2+1`); expected sharp except `m=2`. |
| 2015 | Huang–Petit–Shinohara–Takagi | ePrint 2015/358 | `F_{2^n}` | "Generalized first fall degree assumptions"; refines when the assumption can be relied on. |
| 2016 | Courtois | ePrint 2016/003 | `F_{2^n}` | "On splitting a point with summation polynomials" — critical look at the `S_3`-chain. |
| 2016 | Petit–Kosters–Messeng (PKM) | PKC 2016 | **`F_p`** | Factor base `{(x,y):L(x)=0}`, `L` a composition of low-degree rational maps → low-degree decomposition system. "Practical only for small params"; **asymptotics "limited by our understanding of Gröbner basis algorithms."** |
| 2016 | Galbraith–Gaudry | DCC 78(1) / 2015/1022 | survey | Prime-degree binary curves **"not threatened … still safe."** Reviews the complexity claims critically. |
| 2016 | Momonari–Kudo | (verify) | `F_p` | Improves PKM16 (cited secondhand). |
| 2017 | McGuire–Mueller | ePrint 2017/1262 | `F_p` | Fast `S_m` evaluation, **no Gröbner basis**; sidesteps `d_reg` by brute evaluation (exponential). |
| 2018 | Amadori–Pintore–Sala | FFA 51 / 2017/609 | `F_p` | Needs a **single Gröbner basis computation**; claims to beat Semaev's and PKM16. |
| 2018 | Huang–Kosters–Yang–Yeo | JSC 87 / arXiv 1505.02532 | Weil descent | Non-heuristic last-fall-degree method, between full Gröbner and `d_ff`. |
| 2020 | Trimoska–Ionica–Dequen | AFRICACRYPT 2020 / 2019/313 | `F_{2^n}` | WDSat + **symmetry breaking → `m!` speedup** on the `(m+1)`-th summation polynomial; up to 300× Magma F4. |
| 2020 | Trimoska–Ionica–Dequen | arXiv 2001.11229 | `F_{2^n}` | XOR/parity reasoning (XG-ext) captures a cancellation algebraic methods use but SAT-GE misses. |
| 2020 | Galbraith–Granger–Merz–Petit (GGMP) | SAC 2020 / 2020/1315 | subfield/Koblitz | Frobenius-invariant factor base: ~`1/n` fewer systems, `n^2` linear-algebra speedup. Polynomial-solving savings **experimental only**. |
| 2020 | Huang et al. | JMC 14 | `F_{2^n}` | Quasi-subfield-polynomial factor bases; conditional speedup if good polynomials exist. |
| 2021 | Euler–Petit | FFA 75 | `F_{2^n}` | New quasi-subfield families; **negative** for speedups. |
| 2021 | Tsakou–Ionica | ePrint 2021/721 | genus 1,2 | Effective **endomorphisms shrink the factor base**; net gain despite per-element cost. |
| 2021 | Caminata–Gorla | WAIFI 2020 (LNCS 12542) / 2017/593 | general | **`d_solv ≤` CM-regularity** of the homogenized ideal; field equations ⇒ generic coordinates. |
| 2023 | Caminata–Gorla | JSC / 2021/1611 | general | Links `d_solv`↔last fall degree and `d_reg`↔CM-regularity rigorously. |
| 2024–25 | Hu | ePrint 2024/1923, 2025/015 | `F_p` | Low-impact; concludes **no subexponential method known** for prime-field ECDLP. |
| 2025 | Malhou–Perret–Lauter (HATSolver) | ePrint 2025/2219 | generic Gröbner | ML (hierarchical-attention transformers) solve larger zero-dim systems; **not** applied to summation systems. |
| 2026 | Dryło | ePrint 2026/1299 | any `E` | A degree-2 `G`-compression `= f∘Φ` with `Φ` a separable isogeny with kernel `G`. Certifies "FGHR symmetrisation = Semaev polynomial of `E/G`." |
| 2026 | Wu; Bariant; La Scala–Marchesin–Tiwari | arXiv 2608.01675; ePrint 2026/1927; arXiv 2609.18866 (ids per repo W8, **verify**) | models / splittings | Binary `Z/4` state model; Frobenius-power auxiliary variables; convolution splitting of `e_i`. Repo W8 measured **all three as neutral or negative** on binary descent systems (§4). |

---

## 2. Per-paper dreg claim (the exact statement relevant to solving degree)

- **Semaev 2004 (2004/031).** `S_{m+1}=0` has a root iff `R` decomposes over the factor base. `S_m` symmetric, `deg_{x_i} S_m = 2^{m-2}`. No solving-degree claim is proven; the paper only *hopes* the systems are easy.
- **Gaudry 2009 (JSC 44).** For fixed small `n`, the descent system has bounded degree, so `d_reg` is bounded and `Õ(q^{2-2/n})` follows heuristically. The key is that `n` (not the field) is the degree driver.
- **Diem 2011 (Math. Comp.).** Makes Gaudry rigorous for fixed `n≥2`: the Macaulay degree is a function of `n` only, so no dreg heuristic is required. **This is the one regime where dreg is provably controlled.**
- **Diem 2011 (Compositio).** Subexponential for `n→∞` slowly; again the control comes from the geometry of the restriction-of-scalars variety, not a dreg assumption. Does **not** reach prime-degree binary or prime fields.
- **FPPR 2012 (2012/146).** The Boolean descent of `S_{m+1}` over `V` (dim `n'`) has `d_ff` close to the input degree; **assumes** `d_reg ≈ d_ff`. This is the load-bearing assumption of the whole binary subexponential program.
- **Petit–Quisquater 2012.** Formalizes the assumption: `d_reg` grows like `~m^2` (constant in `n`), giving `2^{c n^{2/3} log n}`. If `d_reg` instead grows with `n`, the result collapses.
- **HPST 2013 / 2015 (2015/358).** `d_ff` assumption "seems to hold" in their symmetrized experiments; complexity improvements are **conditional** on it.
- **FGHR 2014 / FHJRV 2014.** The dreg lever that actually works: quotient by `S_m` (symmetric functions) and by torsion/`[-1]`. Groups of order `2^{m-1} m!` (2-torsion) and `4^{m-1} m!` (4-torsion). Rewriting in invariants lowers the **Bézout number and the Macaulay degree** of the system (not merely the variable count). FHJRV Lemma 6 is the structural bound: *every* degree-2 coordinate is a homography of `x`, so coordinate choice alone cannot help beyond passing to a quotient curve.
- **Galbraith–Gebregiyorgis 2014 (2014/806).** Binary-Edwards invariant variables **lower the degree of the summation polynomials**; the symmetry-breaking factor base raises relation probability. Measured conclusion: these gains do not overturn rho for prime-degree binary.
- **Semaev 2015 (2015/310).** The `S_3`-chain keeps each link degree-2-per-variable; claims the cubic Boolean system has `d_reg ≤ 4` under a first-fall assumption, giving `2^{1.69√(n ln n)}`.
- **Karabina 2015 (2015/319) / Courtois 2016 (2016/003).** Same chain idea; Courtois is skeptical of the implied savings.
- **HKY 2015 (2015/573).** The central disputed result (§3): `d_lf` (order-independent) is the right invariant, and the `d_reg ≈ d_ff` assumption is **doubted with explicit counter-construction**.
- **Kosters–Yeo 2015; Kousidis–Wiemers 2015 (2015/1121).** `d_ff` of the `S_{m+1}` descent `= m^2−m+1` (proven bound). Kosters–Yeo give evidence `d_reg` **grows with `n`** — i.e. `d_reg > d_ff` asymptotically, exactly what breaks FPPR/PQ.
- **PKM 2016.** Over `F_p` the factor base is `L(x)=0` for `L` = composition of small rational maps; the decomposition system's degree is controlled by the **depth/degree pattern of `L`**, but the asymptotic `d_reg` is openly unknown.
- **Amadori–Pintore–Sala 2018 (2017/609).** One Gröbner basis; the `d_reg` of that single system is the whole cost and is not bounded theoretically.
- **Caminata–Gorla 2021/2023.** The rigorous frame: `d_solv ≤ reg(I^{hom})`, and (2023) `d_solv`↔`d_lf`, `d_reg`↔CM-regularity. **Any dreg claim should be stated as a CM-regularity or last-fall-degree statement**, not a first-fall one.
- **GGMP 2020; Trimoska et al. 2020.** Do not change `d_reg`; they cut the *number* of systems (Frobenius ÷`n`) and the per-degree linear-algebra cost (`m!` symmetry breaking, XOR reasoning). Orthogonal to dreg.

---

## 3. The disputed claim: `d_reg ≈ d_ff` (the first-fall-degree assumption)

**What is assumed.** FPPR 2012, Petit–Quisquater 2012, HPST 2013, Semaev 2015 all need: the degree the Gröbner solver actually reaches (`d_solv`, governed by `d_reg`) is within `O(1)` of the first fall degree `d_ff`, which for these systems is small and roughly constant in `n` (`d_ff = m^2−m+1`, Kousidis–Wiemers). If true, the binary PDP is subexponential.

**Why it is doubted.**
- **HKY (CRYPTO 2015, 2015/573):** builds a summation-polynomial Weil-descent system in which the assumption is *unlikely to hold*, and gives computational evidence in char 2 that `d_reg` and `d_ff` diverge. Proposes the **last fall degree** as the honest, monomial-order-independent invariant; rigorously bounding it for Weil-descent systems is left open.
- **Kosters–Yeo (1503.08001):** direct evidence that `d_reg` of the `S_3` systems **increases with `n`**, so the gap `d_reg − d_ff` is not `O(1)`.
- **Ding–Schmidt 2013:** general warning that `d_ff` can be driven artificially low by a planted syzygy, so a small measured `d_ff` is *not* evidence of a small `d_solv`.
- **Caminata–Gorla:** the clean replacement — `d_solv` is bounded by CM-regularity, and equals the last fall degree under their hypotheses. This makes "measure `d_ff`, conclude subexponential" logically invalid.

**Repo's own measurements agree with the skeptics.** The repo measured, on chained-`S_3` binary descent systems (`research/notes/ecc2k130/RESEARCH_ECC2K130_DESCENT_DEGREE.md`, `RESEARCH_KOBLITZ_SCALING_TARGET.md`): `d_ff = 3` everywhere (the Kosters–Yeo trace equation), while the **refutation/solving degree grows with the factor-base dimension `ℓ`**: `m=3` fit `⌈ℓ/2⌉+4`; `m=4` fit `ℓ+4`; degree reached 7 at `(n,ℓ)=(15,5)`. This is a concrete, reproduced instance of `d_reg ≫ d_ff` and growing — the single most important empirical constraint for this report.

**Status.** No proof *or* definitive refutation of the assumption for summation systems has appeared since 2015. It is "doubted, not dead." The consensus (Galbraith–Gaudry survey) is that prime-degree binary curves are safe because the assumption probably fails.

---

## 4. What the repo already has, and the gaps it missed

**Already covered (do not duplicate):**
- Frobenius-invariant factor bases (GGMP) with the structural kill for prime `n`: Frobenius-stable `F_2`-subspaces are cyclic codes, whose dimensions are subset-sums of 2-cyclotomic coset sizes; at `n=131` the only sizes are `0,1,130,131` → empty at the sizes that matter. Same for `n=233,283` (prime, 2 of large order). (`RESEARCH_KOBLITZ_INDEX_CALCULUS.md`.)
- Couveignes–Lercier elliptic-period bases: closed twice — need an auxiliary `E'/F_2` with `n | #E'(F_2)≤5`, empty for `n≥7`.
- FGHR/FHJRV symmetrization fully reproduced and *correctly re-attributed* (`RESEARCH_EXOTIC_COORDINATES.md` §19): every degree-2 coordinate is a Möbius/homography frame (their Lemma 6); the measured 30–416× "symmetry" gain is mostly the change of curve to `E/⟨T⟩` (certified against Dryło 2026/1299, W8 §1), plus a genuine `s`-cover at `m=3`.
- The `S_3`-chain ("splitting trick") attributed to HKY/HPST/Karabina/Semaev 2015.
- Equivariant F4 over the Koblitz `C_n` action (lead "G2"): Faugère–Svartz character-block decomposition, verified exact at `n=7,11`, rank-only speedups ≤128× at `n≤23` — **left open, possibly mislabeled closed.**
- `D_4=(Z/2)^3⋊S_4` symmetry on `S_5` over `F_{p^4}` (lead "C8"): **Macaulay degree 29→14, Bézout 4096→512, `C_4` ÷47–48×** — the repo's strongest measured dreg lever; crossover model moves `2^150→2^135`.
- Crossbred/mutant-XL (monomial dropping): measured, margin shrinks with size, no `(D,k)` filter found.
- Coppersmith/lattice small roots for `S_m` over `F_p`: exact LP ceiling `δ* = 1/6, 1/24, 1/80` (`S_3,S_4,S_5`); "the Coppersmith route to the `m=4` oracle is closed" for single-polynomial triangular lattices.
- Frobenius-power auxiliary variables (Bariant-style): **neutral** over `F_2` because `x↦x^{2^k}` is `F_2`-linear (W8 §3). Convolution `e_i` splitting (La Scala): **negative**, 1.2–2× more work (W8 §4). Wu `Z/4` model: neutral (W8 §2).
- PKM tower oracle over `F_p`: `d_solv` plateaus at 6–7 for `N≤22`, memory-bound at `N≥21`; inconclusive whether `d_solv` stays bounded.

**Gaps the repo appears to have missed (checked against `research/lit2026_*`, `papers/*/refs.bib`, `research_notes/`):**
1. **No Hilbert-series / `d_reg`-minimizing search over the factor-base subspace `V`.** The repo chose `V` for Frobenius-invariance or pipeline cost, never to minimize the *measured solving degree*. Caminata–Gorla's CM-regularity frame is cited nowhere in the repo bib (absent from all four `.bib` files; `Caminata` appears in only 7 files, `Gorla` 18, none as a design principle). **This is the cheapest unexplored lever** (§5 idea 1).
2. **No full invariant-ring (secondary-invariant) rewrite for the non-abelian wreath/reflection groups.** The repo does `C_n`-equivariant F4 (abelian, Faugère–Svartz) and reproduces FHJRV invariants, but does not push to `S_m ≀ (Z/2)` or the full `E[2]`-reflection group *as a Macaulay-matrix block decomposition* (Hubert–Labahn / Steidel style for non-abelian groups). 
3. **The `j=0`/`j=1728` extra automorphism is used only as a GLV `ψ`-fold on extension fields (`÷3`) and as a 3-torsion frame (negative at `m=3`).** It is **not** used as an order-6/order-4 *symmetry quotient of the prime-field PDP ideal* à la C8. The C8 result shows this class of lever is the repo's best; it was never tried on the named prime-Koblitz / `j=1728` curves at `m=2,3`.
4. **No Artin–Schreier trace-elimination pre-reduction.** The repo knows `d_ff=3` comes from the Kosters–Yeo trace equation `z^2+z=c` but never *uses* AS theory to eliminate the trace part and solve on the trace-zero residual.
5. **No theta/biextension (Kummer-variety) PDP encoding.** The repo's splittings are resultant-chains and convolution; the degree-2 Riemann/theta addition law (keeps per-variable degree 2 for all `m`) is untried. Gorla–Massierer trace-zero compression is in the bib-adjacent notes but not as a PDP encoding.
6. **Lattès maps / isogeny-kernel polynomials as the rational map `L` in PKM.** The repo tried PKM Kummer/Dickson towers and isogeny fibers (bounded ×4) but not `L` = kernel polynomial of a small isogeny or a Lattès self-map of `P^1`.
7. **AG/algebraic-geometry-code decoding view of relation search** — absent entirely (only the cyclic-code = Frobenius-subspace connection via quasi-subfield polynomials is present).
8. **Tropical / Newton-polytope-guided monomial and variable ordering** — absent (the repo uses MinTermWeight/MOM for SAT branching only).

---

## 5. Exotic levers, ranked by (plausibility of `d_reg` reduction × cheapness to test)

Each lever states: the idea, closest prior art (novelty), expected effect on `d_reg`/`d_solv`, risk, and **the single cheapest falsification experiment** on toy sizes (`n≤40` binary or `p≤2^20` prime). Ranking is the product of (plausibility it reduces solving degree) × (cheapness), so cheap-and-decisive beats expensive-and-maybe.

### Rank 1 — `d_reg`-minimizing factor-base subspace search (Hilbert-series targeting)
- **Idea.** Treat the factor-base subspace `V` (binary) or the rational map `L` (prime) as a *free design variable* and pick it to minimize the measured solving degree / to push the descended ideal toward semi-regularity (flat Hilbert series). Quantify each candidate by its Hilbert series truncated at low degree, selecting the `V` whose first Hilbert-series sign change is latest.
- **Novelty.** Genuinely open. GGMP and the repo chose `V` for invariance/cost; FPPR fix an arbitrary `V`. No paper optimizes `V` against `d_reg`, and Caminata–Gorla's CM-regularity frame is unused in this area.
- **Expected effect.** If any `V` of fixed dimension `ℓ` yields `d_solv` even 1 below the median, the `⌈ℓ/2⌉+4` growth is not intrinsic and the barrier softens; if the min over hundreds of `V` equals the median at every `ℓ`, the growth is a property of the ideal class and the whole family is closed more firmly than before (a publishable negative).
- **Risk.** Medium–high that min = median (growth likely intrinsic). But either outcome is informative.
- **Cheapest falsification.** Binary `n∈{12,15,18}`, `m=3`, `ℓ∈{3,4,5}`: sample ~300 random `V` (and all Frobenius-stable ones), build the descended `S_3`-chain ideal, compute the solving/last-fall degree with the repo's `solving_degree`. **Falsifier:** `min_V d_solv = median_V d_solv` for every `(n,ℓ)`. ~1 day on existing Rust harness; no new math.

### Rank 2 — Order-6 / order-4 automorphism quotient of the PDP ideal on `j=0` / `j=1728` (prime field)
- **Idea.** On `y^2=x^3+b` (`p≡1 mod 3`) the automorphism `ζ:(x,y)↦(ω x, -y)` (order 6 with `[-1]`) acts on summation variables by `x_i ↦ ω x_i`. Build `S_{m+1}` in `⟨ζ⟩`-invariants (power sums of `x_i^6` etc.), on a factor base `L(x)=0` chosen `ζ`-stable, and solve the quotient ideal — the C8 recipe but with the extra prime-field automorphism instead of `S_4`.
- **Novelty.** The automorphism is used in the literature only for rho (GLV) and, in the repo, as a `ψ`-fold on *extension* fields (`÷3`) and a (negative) 3-torsion frame. Using it as a Bézout-reducing symmetry quotient of the *prime-field* PDP ideal is close to C8/FGHR in spirit but not instantiated for these curves.
- **Expected effect.** By the C8 precedent (`D_4`: Macaulay 29→14), a `Z/6` grading should cut Bézout by up to `6` and the Macaulay degree measurably. Will **not** cross rho (prime-field PDP loses by orders of magnitude) but is a real `d_reg` reduction, which is what the report asks for.
- **Risk.** Medium. The factor base must be both `ζ`-stable and large enough; the `Z/6` grading costs F4 multiplications (the repo saw `Z/3` cost 12–30×), so net wall-time may not drop even as degree does.
- **Cheapest falsification.** `p≡1 mod 3`, `p≤2^16`, `j=0` curve, `m=2`: compare `d_solv` and Bézout of the plain `S_3` vs the `⟨ζ⟩`-quotient `S_3` on a `ζ`-stable base. **Falsifier:** quotient `d_solv` ≥ plain `d_solv` and Bézout drop `< 2×`. ~1–2 days in Sage/repo F4.

### Rank 3 — Wreath-product / full-reflection-group equivariant solve (push symmetrization past FHJRV)
- **Idea.** Decompose the Macaulay matrix of the symmetrized `S_{m+1}` into isotypic blocks for the **full** acting group — `S_m ≀ (Z/2)` from `[-1]` on each summand, or the `E[2]`-reflection group `(Z/2)^{n-1}⋊S_n` — using representation-theoretic block diagonalization (Faugère–Svartz for the abelian part, Hubert–Labahn / Riener-style for the non-abelian part), and work in secondary invariants so the *effective* degree in the invariant generators is minimized.
- **Novelty.** Faugère–Svartz (abelian) and FHJRV (writes invariants, computes `S_8`) exist; the repo's G2 does the abelian `C_n` block split. The non-abelian wreath block-decomposition of the summation Macaulay matrix is only partly realized.
- **Expected effect.** Block diagonalization per se reduces per-degree linear-algebra cost (like G2's ≤128×), but solving in the invariant ring can lower the *degree* when the group order is large (FHJRV ÷8 factor base; C8 degree 29→14). Plausible 1–2 degree reduction at `m=4`.
- **Risk.** Medium. Non-modular condition: need `char ∤ |G|` for clean semisimple blocks — fails for the `Z/2` parts in char 2 (modular!). So on binary curves only the `S_m` (odd-order-safe) part block-splits cleanly; the `Z/2` torsion part is modular and needs Specht/James modular representation theory. This is a real obstruction the repo has not confronted.
- **Cheapest falsification.** `F_{p^2}` (odd `p≤100`) so the group order is invertible, `m=3`, compare `d_solv` of `S_4` solved (a) plainly, (b) in `S_3`-invariants only, (c) in full `S_3≀(Z/2)` invariants. **Falsifier:** (c) `d_solv` = (b) `d_solv`. ~2–3 days (needs an invariant-ring rewrite in Sage).

### Rank 4 — Artin–Schreier trace elimination → solve on the trace-zero residual (char 2)
- **Idea.** The descended binary `S_3` system's low `d_ff=3` is the Kosters–Yeo trace equation `z^2+z=c`. Use Artin–Schreier theory to *solve/eliminate* the additive-trace component in closed form (it is `F_2`-linear algebra), leaving a smaller residual system on the trace-zero subvariety, whose `d_reg` may be lower because the "easy" linear part has been removed before Gröbner.
- **Novelty.** Trace-zero IC exists (Gorla–Massierer, Diem); using AS structure to *pre-reduce the descended ideal* is implicit in Kosters–Yeo's analysis but not an algorithmic lever in the repo (which treats `d_ff=3` as "not a weakness").
- **Expected effect.** Uncertain sign: either removes the `d_reg`-driving part (good) or removes the only easy part and leaves the hard core unchanged (null). Worth knowing which.
- **Risk.** Medium–high (could be null by construction).
- **Cheapest falsification.** Binary `n∈{12,15}`, `m=3`: build the descended system, split variables into AS-trace vs trace-zero, eliminate the trace block linearly, measure `d_solv` of the residual vs the full system. **Falsifier:** residual `d_solv` = full `d_solv`. ~1 day on existing harness.

### Rank 5 — Theta / biextension (Kummer-variety) addition-law PDP encoding
- **Idea.** Encode "`P_1+…+P_m = R`" not via the resultant-chained `S_m` (degree `2^{m-2}` per variable) but via the **degree-2 Riemann relations** of the theta/Kummer model, introducing `O(m)` auxiliary Kummer coordinates but keeping *every* equation quadratic for all `m`. Attacks the `2^{m-2}` degree blow-up at its root.
- **Novelty.** The `S_3`-chain already keeps links quadratic via *point* auxiliaries; the theta/biextension encoding is a genuinely lower-degree alternative not in the IC literature.
- **Expected effect.** If the quadratic Kummer system has `d_reg` comparable to the chain's but far lower per-variable degree, the `m=4,5` oracles (where `2^{m-2}` dominates) could drop sharply. 
- **Risk.** High. The repo's convolution-split (W8 §4) and Frobenius-aux (W8 §3) results are a strong prior that *adding auxiliaries that are algebraically dependent does not lower `d_reg`* — it rose or stayed flat. Theta differs (truly lower degree, not a relabel), but the null prior is real.
- **Cheapest falsification.** `p≤2^12`, `m=3`: build `S_4` via (a) resultant chain, (b) Kummer theta relations; compare `d_solv` and Bézout. **Falsifier:** theta `d_solv` ≥ chain `d_solv` with no Bézout drop. ~3–4 days (theta addition law must be coded).

### Rank 6 — Lattès-map / isogeny-kernel non-linear factor base (prime field, PKM-style)
- **Idea.** In PKM's `{L(x)=0}`, choose `L` = the kernel polynomial of a small-degree isogeny `E→E'`, or a Lattès map (the degree-`ℓ^2` self-map of `P^1` induced by `[ℓ]` on `E`). The decomposition system inherits the isogeny's low-degree factorization, giving a *shallow* composition `L = L_k∘…∘L_1` with controlled degree growth.
- **Novelty.** PKM use generic composed rational maps; specializing `L` to isogeny/Lattès structure is new. Repo closed Couveignes–Lercier (empty by Hasse) and isogeny fibers (bounded ×4) but not this.
- **Expected effect.** PKM argue depth-`k` composition gives degree-`d^k` reach with `d_reg` tied to the *small* `d` per layer; isogeny kernels give the cleanest such layers. Could keep `d_reg` bounded as `|F|` grows — the open PKM question.
- **Risk.** High. PKM already explored the general mechanism and found it "practical only for small params"; isogeny structure may not beat generic `L`.
- **Cheapest falsification.** `p≤2^16`, `m=2`: pick `L` from a degree-`ℓ` isogeny kernel (`ℓ∈{2,3,5}`) and from a Lattès `[2]`-map; measure `d_solv` vs a generic-composition `L` of equal total degree. **Falsifier:** isogeny/Lattès `d_solv` = generic `d_solv` at equal degree. ~2–3 days in Sage.

### Rank 7 — Intermediate-field Galois descent with `tau`-orbit "orbit-variables" (composite `n` only)
- **Idea.** For composite `n=n_1 n_2`, descend `E(F_{2^n})` as `E(F_{(2^{n_1})^{n_2}})` and solve the PDP over the *small* extension degree `n_2`, with the factor base taken as `tau`-orbit representatives (one orbit-variable per Frobenius orbit, i.e. work with traces/norms to `F_{2^{n_1}}`). This is the regime where Diem's bounded-degree control holds.
- **Novelty.** For composite `n` this is essentially Gaudry–Diem + GGMP's orbit idea combined; the explicit "single orbit-variable = relative trace" formulation is a small step beyond.
- **Expected effect.** Large and *proven-in-spirit*: the solving degree is controlled by the small `n_2`, not `n` (Diem 2011). 
- **Risk.** Low that it works on composite `n`; **but off-target** — ECC2K-130 (`n=131`), sect233k1 (`n=233`), sect283k1 (`n=283`) all have **prime** `n`, so the lever is empty exactly for the named curves. Include it to (a) validate the harness against Diem's theorem and (b) document precisely why prime `n` is immune.
- **Cheapest falsification.** Binary `n=16=4·4` and `n=12=3·4`, `m=3`: solve the PDP descended to `F_{2^4}`/`F_{2^3}` with orbit-variables vs descended all the way to `F_2`. **Falsifier (of the *applicability*, not the mechanism):** confirm it needs composite `n` by checking `n=13,17` (prime) give no intermediate field. ~0.5 day.

### Rank 8 — Tropical / Newton-polytope-guided variable and monomial ordering
- **Idea.** Use valuations / Newton-polytope (tropical) data of the summation system to choose a monomial order and variable order that minimize the degree the solver reaches.
- **Novelty.** High — no IC paper uses tropical ordering; the closest is learned orderings (Kera NeurIPS 2024; HATSolver 2025/2219), untried on summation systems.
- **Expected effect.** Can lower the *practical* `d_solv` of a specific F4 run toward the floor — but that floor is the **last fall degree, which is provably order-independent** (HKY). So the lever cannot lower the intrinsic barrier; it can only close the gap between a bad run and `d_lf`.
- **Risk.** High that it touches the intrinsic barrier (theory says it cannot). Low cost, so it is a fast, decisive boundary test of "does ordering matter here at all."
- **Cheapest falsification.** Binary `n=15`, `m=3`: run the repo F4 under ~20 monomial/variable orders (tropical-chosen vs lex/grevlex/random), record the max degree reached. **Falsifier:** all orders reach the same max degree (confirming order-independence of the floor). ~0.5 day — the cheapest of all; best run first as a sanity gate.

**Deliberately demoted (no `d_reg` effect, or already measured):** structured/Toeplitz Macaulay linear algebra (lowers per-degree cost, not degree — this is G2's lane); AG-code list-decoding view (novel but almost certainly reframes the generic `S·T^2≳r` bound rather than lowering `d_reg`); crossbred/monomial-dropping (measured, fades with size); x-only/`E/⟨T⟩` 2-torsion quotient (maxed out — proven = Semaev of `E/⟨T⟩`, Dryło 2026/1299 + repo W8).

---

## 6. The three literature facts that most constrain what can work

1. **Solving degree is invariant under linear (incl. `F_2`-linear Frobenius) and Möbius coordinate changes.** Caminata–Gorla: `d_solv ≤` CM-regularity of the homogenized ideal, a coordinate-independent invariant in generic coordinates (which the field equations guarantee). FHJRV Lemma 6: *every* degree-2 coordinate is a homography of the `x`-line. Repo W8 §3: Frobenius-power auxiliaries are `F_2`-linear and change nothing. **Consequence:** no reparametrization lever — tropical ordering (rank 8), Frobenius/linearized auxiliaries, Möbius/Edwards/Hessian frames — can lower the intrinsic `d_reg`. Only genuine *quotients* (isogeny/torsion lowering per-variable degree; automorphism/`S_m` symmetry lowering Bézout) or genuinely *lower-degree encodings* (theta, rank 5) can.

2. **The `d_reg ≈ d_ff` first-fall-degree assumption is false or doubtful; `d_reg` grows while `d_ff` stays ~constant.** HKY (CRYPTO 2015), Kosters–Yeo, and the repo's own measurements (`d_ff=3` fixed vs `d_solv ≈ ⌈ℓ/2⌉+4` growing in `ℓ`) all point the same way. **Consequence:** the barrier is the *growth of `d_reg` in the factor-base dimension `ℓ`*, not a constant; a lever must bend that growth curve. And **measuring `d_ff` is not evidence of success** — every falsification experiment above must measure the full solving/last-fall degree, never just the first fall.

3. **Symmetry/isogeny quotients are the only proven large `d_reg` lever, and even maxed they do not beat rho for the target curves.** FGHR/FHJRV (`÷` group order), Dryło (symmetrisation = Semaev of `E/G`), and the repo's C8 (`D_4`: Macaulay 29→14, Bézout 4096→512) establish the mechanism; Galbraith–Gebregiyorgis + the repo establish that for prime-degree binary and for prime fields, rho remains much faster regardless. For prime `n` (ECC2K-130, sect233k1, sect283k1) the Frobenius-subspace and intermediate-field levers are *structurally empty* (cyclic-code dimension argument; `2` of large multiplicative order mod `n`). **Consequence:** a realistic `d_reg` *result* (not a break) should target either (i) capping the `d_reg(ℓ)` growth via factor-base shaping (rank 1), or (ii) a new, larger symmetry group from the `j=0`/`j=1728` extra automorphisms (rank 2) — and should be reported honestly as a solving-degree reduction inside a family already uncompetitive with rho, not as a threat to the curves.
