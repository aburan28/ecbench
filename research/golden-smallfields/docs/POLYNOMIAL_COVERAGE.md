# Polynomial coverage across all seven fields

The default polynomial battery contains 69 cases across all seven primary
fields. It combines the original 55-case complexity catalog with a separate
14-case field-coverage extension. The original catalog and its existing
benchmark results are unchanged.

| Field | Definition | Original cases | Added cases | Added curve |
| --- | --- | ---: | ---: | --- |
| `b53` | \(\mathbb F_{2^{53}}\) | 0 | 4 | `K53-a1` |
| `b79` | \(\mathbb F_{2^{79}}\) | 0 | 4 | `K79-a1` |
| `b83` | \(\mathbb F_{2^{83}}\) | 49 | 0 | — |
| `p61` | \(p=2^{61}-1\) | 0 | 2 | `p61-G-0` |
| `p64` | \(p=2^{64}-2^{32}+1\) | 0 | 2 | `p64-G-0` |
| `p83s` | \(p=2^{83}-2^{59}-1\) | 0 | 2 | `p83s-G-0` |
| `p83` | \(p=2^{83}-55\) | 6 | 0 | — |
| Total | Seven fields | 55 | 14 | Five added curve selections |

The existing `b83` cohort contains `K83-a0`, `S83-a1` and `B83-a0`; the existing
`p83` cohort uses `p83-S-0`. These larger cohorts exercise factor-base sizes,
basis choices and empty-factor-base diagnostics. The new cases provide small,
deterministic field coverage with the same two-point raw polynomial contract.

```bash
python tools/build_algebra_coverage_cases.py
python -m ecbench algebra-cases
python -m ecbench algebra-verify --cases build/algebra_coverage_cases.json
python -m ecbench algebra-bench --cases build/algebra_coverage_cases.json --method buchberger --order grevlex --timeout 2.5 --output coverage-buchberger.json
python -m ecbench algebra-bench --cases build/algebra_coverage_cases.json --method f5b --order grevlex --timeout 2.5 --output coverage-f5b.json
```

Without `--cases`, the CLI merges the packaged 55-case catalog and packaged
14-case extension. An explicit `--cases` path selects that catalog alone.
For example, `--cases build/algebra_cases.json` selects only the original 55
cases, while `--cases build/algebra_coverage_cases.json` selects only the new 14.
`--curve K53-a1` selects its four cases from the default combined battery.

## Binary extension cases

Each of `K53-a1` and `K79-a1` has four cases: dimensions two and three, with
one planted target and one deterministic subgroup target at each dimension.
The F2-linear basis starts with the explicitly verified independent elements
\(x(G),x([2]G)\). At dimension three, the generator adds the first independent
low polynomial-basis element.

Both unknown x coordinates are parameterized in this basis. The binary Semaev
equation is expanded into all field-coordinate equations over F2, with Boolean
field equations and explicit selected-subgroup membership indicators retained.
The selected point pool contains both signs of every accepted x coordinate.
For these eight fixtures, each selected pool has two distinct x coordinates
and four signed points; no pool is empty.

The planted target is \([3]G=G+[2]G\). Its exact truth is two ordered signed
point decompositions, or one decomposition modulo swapping. The seeded targets
use `SHA256(seed NUL case-label)`, reduced modulo \(r-1\) and incremented by one.
They are deterministic test targets, with no uniform-randomness claim. All four
seeded binary cases have zero relations in their declared factor bases.

## Prime extension cases

Each of `p61-G-0`, `p64-G-0` and `p83s-G-0` receives two cases. Their factor
base consists of the four distinct coordinates \(x([k]G)\), \(1\le k\le4\),
with both signs allowed. The short-Weierstrass Semaev equation is paired with
the two degree-four membership polynomials
\(\prod_{x\in X}(x_i-x)\).

The planted target is \([3]G\). Its exact truth has four ordered signed
decompositions, or two modulo swapping: the pairs \(G+[2]G\) and
\(-G+[4]G\), together with their swaps. The seeded target uses the same
SHA-256 derivation as the existing prime helper, with the actual selected curve
ID in the hash input. Each of the three seeded prime cases has zero relations
in its factor base.

## Generation and validation

`tools/build_algebra_coverage_cases.py` reuses the existing binary polynomial
and point-enumeration helpers and the independent prime arithmetic/Semaev
helpers. It leaves their original defaults intact. The new catalog is
`build/algebra_coverage_cases.json`; it includes exact sparse polynomials,
factor-base points, target derivations, golden roots and signed decompositions.

Generation first cross-checks each Semaev formula with independent point sums.
It then exhaustively re-evaluates the complete sparse system on every declared
finite-domain assignment. The runtime reference arithmetic separately verifies
the signed point sums and exact membership in the selected cyclic subgroup.
All 14 cases pass these checks. Source hashes bind the extension to the primary
field/curve records, and the original 55-case catalog's byte hash is retained.
The builder checks that the original catalog remains byte-for-byte unchanged.
Rebuilding with the same seed produces the same catalog bytes.

Generation and golden-truth validation do not measure Gröbner solving or
degree of regularity. The solver extracts roots from its computed basis and
uses the golden answers only for post-solve conformance checking. The extension
has seven planted successes and seven seeded zero-relation fixtures; those
cohorts should be reported separately. They provide coverage, not an unbiased
relation-yield estimate, a complete index-calculus attack, or evidence about
higher-arity systems.

Keep the original 55-case measurements and the 14-case coverage measurements
as separate result artifacts when comparing solver complexity. The new cases
share the established timing scope: fixture generation is excluded, actual
basis/root/relation work is measured, and the isolated worker timeout also
covers the degree-of-regularity diagnostic when that diagnostic is enabled.
