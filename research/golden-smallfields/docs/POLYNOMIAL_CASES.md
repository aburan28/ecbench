# Curated two-point Semaev polynomial cases

`build/algebra_cases.json` contains **55 deterministic systems** for testing
polynomial solvers and relation enumeration. It is a separate experimental
fixture battery. The seven primary fields, 56 curve records, and their
arithmetic vectors remain unchanged.

Every case seeks a relation `P + Q = R` with two factor-base points and one
fixed finite target. The relation arity is **2**; the summation polynomial is
**f3**, whose third argument is the target's fixed x-coordinate. These cases do
not supply higher-arity systems, a complete index-calculus attack, or an
asymptotic hardness estimate.

## Exact selection

The canonical curve and field records are embedded in each case. The selected
curves are `K83-a0`, `S83-a1`, and `B83-a0` over `b83`, and `p83-S-0` over the
prime field `p83 = 2^83 − 55`.

| Curve | Factor-base style | Dimensions or requested x count | Targets | Cases |
|---|---|---|---|---:|
| `K83-a0` | `low_monomials` | Binary d=2,3,4,5,6 | Planted and seeded | 10 |
| `S83-a1` | `low_monomials` | Binary d=2,3,4,5,6 | Planted and seeded | 10 |
| `B83-a0` | `low_monomials` | Binary d=2,3,4,5,6 | Seeded; empty-base diagnostics | 5 |
| `K83-a0` | `subgroup_seeded_basis` | Binary d=2,3,4,5 | Planted and seeded | 8 |
| `S83-a1` | `subgroup_seeded_basis` | Binary d=2,3,4,5 | Planted and seeded | 8 |
| `B83-a0` | `subgroup_seeded_basis` | Binary d=2,3,4,5 | Planted and seeded | 8 |
| `p83-S-0` | `small_generator_multiples` | Prime-field x count=4,8,16 | Planted and seeded | 6 |

There are **49 binary and 6 prime cases**, comprising **25 planted and 30
seeded-target cases**. The five empty-base diagnostics are included in the
seeded cohort. Dimensions are factor-base dimensions, not field degrees or
rational group invariant factors.

## Factor-base construction and membership

For a binary dimension d, each x-coordinate lies in an explicit d-dimensional
F₂-linear span. `basis_elements` and the two `coordinate_groups` record the
actual field elements, in the exact polynomial basis from the field record.
They are serialized as fixed-width field hex values.

* `low_monomials` uses `1,u,...,u^(d−1)`.
* `subgroup_seeded_basis` starts with `x(G),x([2]G)` and completes an independent
  basis with the first eligible `1,u,u²,...`. It deliberately includes two
  known selected-subgroup points. It is a constructed performance fixture,
  not an unbiased factor-base sample.

For every candidate x, independent polynomial-field arithmetic solves the
curve equation. A point is admitted only if it lies in the selected rational
prime-order subgroup. These three binary records have rank-one rational
r-torsion, so exact `[r]P=O` testing identifies the selected subgroup. Both
signs are retained. The artifact records accepted x values, their binary
encoding masks, all actual points, and counters for candidates with no curve
points or with points outside that subgroup.

The low-basis accepted x counts in this release are:

| Curve | d=2 | d=3 | d=4 | d=5 | d=6 |
|---|---:|---:|---:|---:|---:|
| `K83-a0` | 1 | 3 | 5 | 8 | 14 |
| `S83-a1` | 2 | 3 | 6 | 10 | 21 |
| `B83-a0` | 0 | 0 | 0 | 0 | 0 |

No planted relation is invented for an empty factor base. The `B83-a0`
low-basis cases contain an explicit constant membership equation `1=0`, and
are unsatisfiable diagnostics. The subgroup-seeded basis supplies genuine
planted cases for the same curve; it retains two accepted x values at each
selected dimension.

For `p83-S-0`, the requested x values are `x(kG)` for `k=1,...,4`, `1,...,8`,
or `1,...,16`, with both point signs admitted. These are exact finite lists
of subgroup points. There is no claimed prime-field bit restriction or binary
linear-span interpretation.

## Polynomial systems

For the ordinary binary model `y²+xy=x³+a*x²+b`, use

```text
f3(x1,x2,z) = (x1*x2 + x1*z + x2*z)^2 + x1*x2*z + b.
```

Substitute the finite target `z=x(R)` and expand x1 and x2 in their recorded
binary bases. The result gives all m=83 coordinate equations over F₂. Squares
of Boolean variables are reduced using `v²=v`; explicit field equations
`v²+v=0` are also included. The resulting f3 coordinate equations have degree
at most two, although the original summation polynomial has total degree four.

Separate `factor_base_membership` equations restrict each x to the actual
admitted subgroup x values. They use Boolean indicator polynomials: a tuple
indicator is a product of `v_i` or `1+v_i`, and membership is
`1 + sum(allowed tuple indicators) = 0`. Their degrees can reach d. A measured
input or regularity degree for the full system must not be attributed solely
to the Semaev coordinate equations. Roles and degrees are recorded separately
for `semaev_f3_coordinate`, `boolean_field_equation`, and
`factor_base_membership`.

For the prime short-Weierstrass model `y²=x³+A*x+B`, use

```text
f3(x1,x2,z) = (x1-x2)^2*z^2
              - 2*((x1+x2)*(x1*x2+A) + 2*B)*z
              + (x1*x2-A)^2 - 4*B*(x1+x2).
```

Substitute the same finite-target convention. The remaining two variables
are constrained by `product(x1-x_i)=0` and `product(x2-x_i)=0`, using the
recorded finite factor-base list. Those membership polynomials have degrees
4, 8, or 16; their degree is separate from f3's degree four.

Polynomial coefficients are canonical decimal strings modulo
`coefficient_modulus` (`2` or the recorded prime). A term's nonnegative integer
exponents are in the exact `variables` order. An empty term list denotes zero.
For binary cases the order is the x1 bits followed by the x2 bits; prime cases
use `x1,x2`. Coordinate groups specify how assignments become field elements.

## Target cohorts and independent truth

A planted target is the sum of two actual retained points. Subgroup-seeded
binary bases and prime small-multiple bases use `G + [2]G = [3]G`. Low-basis
binary cases retain a concrete available pair, allowing equal points when
needed, and record that planted decomposition explicitly.

The other target is a deterministic seeded scalar multiple of G, generated
independently of relation search. The artifact records the seed, scalar, and
sampling method. SHA-256 output is reduced modulo `r−1` and incremented, so
this is reproducible pseudorandom sampling with reduction bias, not a claim
of exact uniform sampling. Planted and seeded cohorts must remain distinct.

The default battery's 30 seeded targets have zero decompositions. This is a
result for the specified tiny factor bases and targets, not an estimated
cryptographic relation probability. Prime planted cases have 4, 12, and 28
ordered decompositions; a diagonal binary relation can have only one. Counts
must distinguish ordered solutions from the quotient under swapping P and Q.

Golden truth is generated without a Gröbner solver:

1. Enumerate all admitted x pairs and evaluate f3 directly.
2. Enumerate all ordered signed factor-base point pairs and add them using
   independent affine curve arithmetic.
3. Match the sums to the full target point, verify subgroup membership, and
   cross-check the x solutions with the polynomial system.
4. Record exact variable assignments, actual point decompositions, ordered
   counts, and swap-orbit counts.

Boolean variable domains contain `0,1`; there are at most `2^(2d)` raw tuples,
with d≤6. Membership reduces these to the actual accepted x pairs. Prime
finite domains contain the recorded 4, 8, or 16 x values, hence at most 256
ordered x pairs. `expected` is a post-solve conformance oracle, never a solver
output or a shortcut that may be reported as a solved system.

## Reproduction and generation timing

From the project directory:

```sh
python tools/build_algebra_cases.py
python tools/measure_algebra_generation.py
```

`build_algebra_cases.py` reads the primary canonical field/curve records and
uses independent point arithmetic. It supports `--seed` and `--output`.
Changing the seed changes the seeded target cohort and its truth, while the
fixed factor-base rules remain explicit. Default regeneration is byte-identical;
the artifact also records hashes of its primary source files.

`measure_algebra_generation.py` performs an actual complete build in a
temporary directory, times it, and requires its SHA-256 and bytes to match
`build/algebra_cases.json`. It writes
`build/algebra_generation_benchmark.json` without modifying the canonical
case payload. The generator and timing tool use the Python standard library.

The measured generation time **includes golden truth enumeration**. Prime
truth enumeration is separately timed; binary polynomial construction and
truth enumeration are timed together per case. Factor-base construction is
also separately timed. Timing wrappers change no payload values, and their
small overhead is included. Interpreter startup, imports, post-build integrity
comparison, and timing-result assembly are excluded.

This measurement is fixture-preparation cost. It is separate from solver
construction, Gröbner computation, solver-root enumeration, and relation
verification on solver output. No generation or golden-truth time may be
presented as an index-calculus solve time. Backend benchmarks must state their
own algorithm, ordering, limits, timeout behavior, and regularity definition.
