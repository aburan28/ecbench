# ecbench: golden polynomial and elliptic-curve benchmarks

A reproducible battery for polynomial **index calculus** experiments: relation
yield, Gröbner bases, regularity diagnostics, and measured solve time. It includes
**69 explicit polynomial systems across all seven fields**, exact signed point decompositions, planted
and independently sampled targets, and recorded Buchberger/F5b baseline runs.
The current systems search for two-factor-base-point decompositions using the
third summation polynomial; they are relation-search experiments, not a complete
discrete-logarithm solver.

The [complete PDF report](docs/ecbench_complete_report.pdf) combines the
field/curve analysis with a six-page appendix on the benchmark battery. The
[appendix](docs/ecbench_benchmark_addendum.pdf) is also available separately.

The arithmetic foundation contains **7 fields, 56 explicit curves, 6,931 known-answer vectors,
418 primality certificates, 741 named points, and 52 represented isogeny
classes**. The canonical records include complete curve coefficients, point
orders, subgroup generators, group structure, twists, embedding degrees, and
endomorphism metadata.

Separate fixtures classify all 56 curves for GLV2/Frobenius use: 30 supported
cases with 480 scalar-decomposition vectors. An additional curve over an
82-bit quadratic field supplies genuine GLV/GLS4 tests on an 80-bit prime
subgroup, with proven 20-bit component bounds. This supplement does not change
the primary 56-curve corpus.

These field sizes are for arithmetic experiments and benchmarks. They do not
provide modern production cryptographic security. A large embedding degree or
CM discriminant does not change that limit.

## Quick start

Run commands from this directory with Python 3.10 or later. The reference
arithmetic, CLI, and default verification use only the Python standard library.
Working from the checkout needs no installation. To install it for use
elsewhere, run:

```sh
python -m pip install .
```

Source installation uses the declared `setuptools>=68` build backend. A supplied
wheel can instead be installed directly:

```sh
python -m pip install /path/to/ecbench_golden-1.0.0-py3-none-any.whl
```

After installation, the `ecbench` command is also available. The module runner
works both from the checkout and after installation:

```sh
python -m ecbench list
python -m ecbench show --curve p83-S-0
python -m ecbench verify
python -m ecbench verify --with-pari
python -m ecbench bench --curve p83-S-0 --samples 5 --iterations 3 --output results.json
```

Install the optional polynomial baseline and run a bounded experiment:

```sh
python -m pip install '.[algebra]'
python -m ecbench algebra-cases
python -m ecbench algebra-verify --output fixture_verification.json
python -m ecbench algebra-bench --method buchberger --timeout 5 --max-cases 6 --output polynomial_results.json
python -m ecbench algebra-bench --curve S83-a1 --method f5b --timeout 10 --output f5b_results.json
python -m ecbench endomorphism-verify
python -m ecbench endomorphism-bench --curve p83-E4-2 --curve K83-a0 --output glv2_results.json
python -m ecbench glv4-verify
python -m ecbench glv4-bench --output glv4_results.json
python -m ecbench battery-verify --output verification.json
```

For a supplied wheel, use `python -m pip install '/path/to/wheel.whl[algebra]'`.
Each polynomial solve runs in an isolated subprocess. A timeout is censored;
it is never recorded as a completed zero-relation solve. Successful roots are
lifted to rational points and checked against the exact selected subgroup and
the signed target equation.

The recorded `build/algebra_buchberger.json` and `build/algebra_f5b.json` runs
used a 2.5-second wall budget per case. Buchberger completed 53 of 55 cases;
F5b completed 51. Every completed result matched its golden truth. See
[measured results](docs/BENCHMARKS.md), [case definitions](docs/POLYNOMIAL_CASES.md),
and [metric definitions](docs/ALGEBRA_METRICS.md) before comparing timings.
The additional 14 coverage cases completed with both methods and bring the
default battery to 69 cases. Their separate result files preserve the original
55-case measurements. [Coverage details](docs/POLYNOMIAL_COVERAGE.md) identify
the exact curves and constructions.

Final basis degree is reported separately from the highest-degree-input-part
ideal's regularity diagnostic. SymPy does not expose its solving degree, so
that metric is explicitly unavailable. Solve timing starts from the exported
sparse fixture; dataset generation and golden truth enumeration have their own
recorded scope.

A concrete Sage/Singular adapter is included for characteristic-two systems:

```sh
python -m ecbench backend-export --case K83-a0-f3-low_monomials-d2-planted --destination external_case
sage -python external_case/solver.py > external_case/basis.json
python -m ecbench backend-replay --case K83-a0-f3-low_monomials-d2-planted --result external_case/basis.json --output external_case/replay.json
```

Sage/Singular was not installed for the supplied baseline measurements. This
adapter explicitly rejects characteristics above `2^31−1`; the 83-bit prime
cases use the arbitrary-characteristic SymPy baseline. External replay checks
the exact finite-domain solution set and point relations; it does not claim a
certificate of polynomial-ideal equality. See [backend protocol](docs/BACKEND_PROTOCOL.md).

The optional `--with-pari` check requires PARI/GP. Put `gp` on `PATH`, or set
`CURVE_REPORT_GP` to its executable path. It independently recomputes curve
cardinalities, rational group structures, named-point orders, subgroup-generator
orders, and successful curve-arithmetic answers. Field, encoding, and membership
vectors are checked by the standard-library reference verifier.

`p83-S-0` is the screened prime-order example
`y² = x³ − 3x + 154` over `p = 2^83 − 55`. A binary example is `S83-a1`,
with `y² + xy = x³ + x² + u^251` over the recorded degree-83 polynomial field.
The exact binary coefficient `u^251` has integer bitmask `9925396` in that basis.

## Fields and curve families

| Field ID | Definition |
|---|---|
| `p61` | Prime field, `p = 2^61 − 1` |
| `p64` | Prime field, `p = 2^64 − 2^32 + 1` |
| `p83s` | Prime field, `p = 2^83 − 2^59 − 1` |
| `p83` | Prime field, `p = 2^83 − 55` |
| `b53` | `F₂[u]/(u^53 + u^6 + u^2 + u + 1)` |
| `b79` | `F₂[u]/(u^79 + u^9 + 1)` |
| `b83` | `F₂[u]/(u^83 + u^7 + u^4 + u^2 + 1)` |

The 42 prime-field records cover special-`j` families (`E3` for `j=0`, `E4`
for `j=1728`), generic examples (`G`), their twist representatives, and the
screened 83-bit prime-field pair (`S`). The 14 binary records comprise six
Koblitz examples (`K`, coefficient `b=1`), six B-style examples (`B`, `b=u`),
and the screened degree-83 pair (`S`, `b=u^251`). These are explicit examples;
the dataset does not enumerate every curve over these fields, and the B-style
examples are not NIST parameter sets.

Every curve records all five Weierstrass coefficients `a1,a2,a3,a4,a6`. The
supported models are prime-field short Weierstrass curves and ordinary binary
curves `y²+xy=x³+a2*x²+a6`.

## Files and representation

The authoritative package data is in `ecbench/data/`:

| File | Contents |
|---|---|
| `manifest.json` | Dataset version, counts, payload hashes, curve-record hashes, proof boundaries |
| `fields.json` | Field definitions, encodings, and proof references |
| `curves.json` | Full curve parameters, invariants, selected subgroups, named points, twist metadata |
| `vectors.json` | Field, curve, encoding, membership, and map known answers and rejection cases |
| `relations.json` | Isogeny, geometric-isomorphism, endomorphism-order, and twist groupings |
| `certificates.json` | Flat graph of 418 Lucas primality certificates, with child references by integer value |
| `algebra_cases.json` | 55 sparse polynomial systems, domains, factor bases, targets, and exact golden decompositions |
| `algebra_coverage_cases.json` | 14 additional systems covering the other five fields; the CLI combines both catalogs |
| `endomorphism_cases.json` | All-curve GLV2/Frobenius classification, eigenvalues, kernel lattices, bounds, and vectors |
| `glv4.json` | Supplementary quadratic-field GLV/GLS4 curve, maps, lattices, certificates, and known answers |
| `research_manifest.json` | Research-catalog counts, SHA-256 payload hashes, and coverage limits |

`golden_ecbench.json` combines the canonical records into one document;
`golden_ecbench.jsonl` exposes them as line-oriented records. The published
shape and string formats are specified by `schemas/dataset.schema.json`.
Polynomial cases and measured results have their own schemas in `schemas/`.
The 6,931 vectors comprise 203 field, 2,447 curve, 3,520 encoding, 553
membership, and 208 map vectors. Original arithmetic evidence is retained in
`provenance/`; intermediate vector generator outputs belong in `build/`. See [docs/QUALITY.md](docs/QUALITY.md)
for the detailed conformance contract.

Authoritative mathematical integers use **canonical decimal strings**, including
orders, traces, discriminants, cofactors, prime factors, embedding degrees, and
scalar inputs. This preserves exact values in applications whose JSON numbers
cannot represent these integers. Bounded metadata such as a degree or factor
exponent can be a JSON number.

Field elements and finite-point coordinates use `0x` followed by fixed-width,
lowercase, big-endian hexadecimal. Coordinate width is recorded per field.
A finite point is `{"x": "0x...", "y": "0x..."}`; infinity is JSON `null`.
For a binary field, bit `i` represents the coefficient of `u^i`; field addition
is XOR, and multiplication reduces modulo the exact recorded polynomial.
This is not arithmetic modulo the integer `2^m`.

SEC 1 octet strings use lowercase hex **without** `0x`. The accepted forms are
compressed (`02`/`03`), uncompressed (`04`), and infinity (`00`). For binary
compression, the selector is the low bit of `y/x` when `x≠0`, and is zero when
`x=0`. Decoding strictly rejects hybrid forms, unsupported prefixes, wrong
lengths, noncanonical coordinates or high unused bits, off-curve points,
unsolvable compressed coordinates, and the binary `x=0` selector `03`.
Decoding does not silently reduce coordinates.

## Point orders and subgroup semantics

For a rational point group `Z/d1 × Z/d2`, records give `d1 | d2`, the group
order `N=d1*d2`, and the group exponent `d2`.

* `G` has exact selected prime order `r`. Its deterministic construction uses
  projection multiplier `d2/r`, with ascending encoded `x` and a specified
  choice of `y`. The first nonzero successful projection is retained.
* `R` has exact order `d2`. It generates the **whole rational group only when
  the group is cyclic**, meaning `d1=1`. It does not constitute a full basis
  for a noncyclic group.
* Other named points have independently verified exact orders. Infinity has
  order one. The binary point `T=(0,sqrt(b))` has order two.

The general whole-group scalar multiplication vectors include negative, zero,
small, order-boundary, field-bit-width, and larger-than-field-size scalars.
An arbitrary point's scalar must not be reduced modulo the selected `r`:
`[r]P` can be nonzero when `P` is outside the selected subgroup.

The package distinguishes on-curve membership, rational `r`-torsion, membership
in the particular cyclic subgroup `<G>`, and public-key validity. If rational
`r`-torsion has rank one, on-curve membership plus `[r]P=O` determines `<G>`.
For the recorded rank-two cases, the reference implementation uses exact
baby-step/giant-step membership in `<G>`; an independent `r`-torsion point need
not belong to that subgroup. Infinity belongs to a subgroup mathematically,
but the public-key predicate rejects it. These predicates are experimental
conformance rules, not a complete protocol specification.

## Arithmetic evidence and mathematical limits

Field primality is supported by Lucas certificates using complete
factorizations of `n−1`. The 418 certificates are stored as a flat proof graph: a
prime factor refers to the certificate with that integer value, and verification
follows those references down to the base case 2. Binary moduli are checked with Rabin irreducibility
identities. Curve cardinalities and rational group structures originate from
PARI/GP computations; Koblitz extension traces also have an exact Frobenius
recurrence. Independent Python vector generators use modular arithmetic and
Tonelli–Shanks over prime fields, and polynomial arithmetic with odd-degree
half trace over binary fields. The packaged reference verifier is a separate
implementation, and the optional PARI oracle adds another check.

The recorded discriminants have different meanings: the curve discriminant
`delta` is a field element, `D_pi=t²−4q` is the integer Frobenius-order
discriminant, `D_K` identifies the quadratic field, and the endomorphism-order
discriminant includes its actual conductor. A large Frobenius conductor does
not imply the same conductor for the full endomorphism ring. Ordinary
endomorphism rings are identified for these records. Class numbers are `null`
where they were not computed; `null` does not mean zero.

For supersingular records, the geometric endomorphism algebra is quaternionic,
while the endomorphisms defined over the recorded base field form the stated
quadratic order. The dataset identifies the geometric quaternion algebra and
maximality but does not supply an explicit quaternion-order basis.

The 52 represented isogeny classes use **field ID and Frobenius trace**. This
is a class-membership criterion over a fixed finite field. Separate map
vectors give operational known answers for the supplied explicit isogenies
and endomorphisms, including their target coefficients and exceptional
points; they do not compute paths between every pair of isogenous records.
Same `j`, same endomorphism order, abstract group isomorphism, and isogeny are
different relations. Twist partners may be isomorphic
representatives of the standard twist; the metadata records the representative
coordinate scale where required for special-`j` families. Self-isomorphic
quadratic twists are retained. Binary relative Frobenius maps `(x,y)` to
`(x²,y²)` on the target with `b` replaced by `b²`. This target is the same
curve for the Koblitz examples; otherwise it is an auxiliary model, not an
additional record among the 56 curves.

The 208 operational map vectors comprise:

| Map | Vectors |
|---|---:|
| Prime-field degree-two isogenies | 20 |
| Gaussian automorphisms | 40 |
| Eisenstein automorphisms | 78 |
| Binary relative Frobenius | 70 |

Their source points, target coefficients, parameters where needed, and exact
point outputs are recorded explicitly. Auxiliary target models do not change
the counts of primary curves, named points, or represented isogeny classes.

Embedding degrees are exact multiplicative orders `ord_r(q)` for prime
subgroups away from the field characteristic. Characteristic factors have no
such degree. Large embedding degrees do not compensate for small subgroup
orders, weak twist factors, or the low field sizes here.

## Rebuilding and benchmarking

The bundled provenance is sufficient to regenerate vectors and canonical data:

```sh
python tools/prime_vectors.py --pari-crosscheck
python tools/binary_vectors.py
python tools/prime_maps.py --pari-crosscheck
python tools/binary_maps.py
python tools/build_dataset.py
python tools/build_algebra_cases.py
python tools/build_algebra_coverage_cases.py
python tools/build_endomorphism_cases.py
python tools/glv4_cases.py --pari-crosscheck
python tools/build_research_manifest.py
python -m ecbench verify --with-pari
python -m ecbench battery-verify
python tools/check_reproducibility.py --gp /path/to/gp
```

The vector generators and assembler use the standard library; the requested
prime cross-check and final independent verification require PARI/GP.
Recomputing the original point counts, factorizations, and source certificates
uses the original computation scripts and may additionally require SymPy.
The supplementary GLV/GLS4 generator also requires SymPy; its optional PARI
cross-check records the actual PARI version. The isolated rebuild check compares
canonical bytes and therefore expects that source-proof version metadata to
match the release.
To regenerate the PDF, install `.[report]` and run
`python tools/build_benchmark_report.py`. The layout uses the DejaVu Sans fonts
from the system font directory. PDF generation is independent of ordinary
arithmetic verification and solver runs.
Neither ordinary package use nor regeneration of the bundled vectors requires
Matplotlib. Generator source paths default to this project's `provenance/`
directory; their output paths default to this project's `build/` directory.

Benchmarks use the supplied deterministic points, including the maximal-order
point `R`, and write machine-readable timing results. Keep Python version,
machine, sample count, iteration count, and selected curve fixed when comparing
runs. The baseline is variable-time Python affine arithmetic: its timings do
not measure optimized implementations, constant-time behavior, side-channel
resistance, or end-to-end protocol performance. They are useful for checking a
local reference baseline and for designing comparable implementation tests.
