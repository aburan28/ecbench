# GLV2 and binary Frobenius fixtures

The supplementary catalog classifies all 56 primary curves. It supplies 30
verified two-term decompositions: 24 prime-field GLV2 cases and six binary
Koblitz Frobenius cases. The other 26 records explain why the selected map or
selected cyclic subgroup is unavailable. The primary parameters and generators
are unchanged.

```bash
python tools/build_endomorphism_cases.py
python -m ecbench endomorphism-cases --supported
python -m ecbench endomorphism-cases --curve p64-E3-2
python -m ecbench endomorphism-verify
python -m ecbench endomorphism-bench --curve p83-E4-2 --curve K83-a0 --samples 3 --output results.json
```

Generation writes `build/endomorphism_cases.json`. The CLI uses the packaged
catalog. Rebuilding the supplementary catalog does not change the primary
dataset. `build/endomorphism_verification.json` records verification of all 480
scalar vectors. `build/endomorphism_benchmark.json` records three timing samples
per operation for seven representative prime and binary cases.

## Maps and scalar eigenvalues

For ordinary prime-field E4 curves, the map is

\[
\phi(x,y)=(-x,iy),\qquad i^2=-1.
\]

The selected unit is the smaller integer square root of \(-1\) in the base
field. Its characteristic polynomial is \(X^2+1\). For ordinary E3 curves,

\[
\phi(x,y)=(\zeta x,y),\qquad \zeta^2+\zeta+1=0.
\]

Here \(\zeta=g^{(p-1)/3}\), where \(g\) is the primitive root in the primary
field record, and the characteristic polynomial is \(X^2+X+1\).

The builder computes the roots of the relevant polynomial modulo the selected
prime subgroup order \(r\), then checks

\[
[\lambda]G=\phi(G).
\]

This equality establishes the scalar action on the entire cyclic group
\(\langle G\rangle\). Checking only \([r]\phi(G)=O\) would establish rational
\(r\)-torsion membership, which is insufficient when rational \(r\)-torsion
has rank two.

For the six binary K curves, \(a\in\{0,1\}\) and \(b=1\), so the curve is
defined over \(\mathbb F_2\). The coordinate-squaring map

\[
\tau(x,y)=(x^2,y^2)
\]

is an endomorphism of the same curve. With \(\mu=2a-1\), verification checks

\[
\tau^2-\mu\tau+[2]=0,\qquad
\lambda^2-\mu\lambda+2\equiv0\pmod r,\qquad
\tau(G)=[\lambda]G,\qquad \tau^m(G)=G.
\]

The eigenvalue also matches
`provenance/binary_data.json:F2_frobenius_scalar_on_largest_subgroup`.
This is a two-term lattice decomposition using \(1,\tau\). It is not a
specialized tau-adic NAF implementation.

## Exact lattice decomposition

The kernel lattice starts with the rows

\[
L=\begin{pmatrix}r&0\\-\lambda&1\end{pmatrix}.
\]

Exact Lagrange-Gauss reduction produces a short basis. Its determinant has
absolute value \(r\), and every row \((u,v)\) satisfies
\(u+\lambda v\equiv0\pmod r\). The verifier separately checks the Gauss
reduction inequalities.

Exact rational nearest-plane reduction decomposes any supplied integer \(k\)
as signed components \(k_1,k_2\), with

\[
k\equiv k_1+\lambda k_2\pmod r,
\qquad [k]G=[k_1]G+[k_2]\phi(G).
\]

The binary cases use \(\tau\) in place of \(\phi\). Rounding chooses the
nearest integer, with ties toward positive infinity. Only this operation on
the selected subgroup reduces \(k\) modulo \(r\); arbitrary-point reference
multiplication retains its existing scalar semantics.

Each supported record contains its generator, map image, eigenvalue, reduced
basis, determinant, exact component bounds and 16 scalar vectors. The vectors
cover zero, one, minus one, \(r-1,r,r+1\), fixed 80- and 128-bit integers, and
eight deterministic subgroup residues. Verification checks both congruence
and signed point reconstruction using ordinary scalar multiplication and an
interleaved two-point multiplication.

Universal nearest-plane bounds are computed from the Gram-Schmidt basis:
each coordinate is bounded by half the sum of the absolute corresponding
Gram-Schmidt coordinates. The catalog records these bounds as exact rational
numbers and integer ceilings. The largest component observed in this battery
is 40 bits; observed values do not replace the per-case proven bounds.

Mathematical integers, including signed components and point coordinates, use
canonical decimal strings. Case point coordinates use decimal strings rather
than the primary dataset's padded hexadecimal storage. Source hashes bind the
catalog to the canonical primary field and curve records.

## Unsupported selections

| Reason | Curves | Meaning |
| --- | ---: | --- |
| `gaussian_unit_not_in_base_field` | 4 | The required square root of \(-1\) is absent from \(\mathbb F_p\). |
| `eisenstein_unit_not_in_base_field` | 2 | A nontrivial cube root of unity is absent from \(\mathbb F_p\). |
| `selected_cyclic_subgroup_not_stable` | 2 | The rational map exists, but its image of \(G\) lies outside \(\langle G\rangle\). |
| `relative_frobenius_targets_distinct_curve` | 8 | Coordinate squaring changes a coefficient, so it maps to a distinct conjugate curve. |
| `no_explicit_nontrivial_endomorphism_in_this_battery` | 10 | This battery supplies no nontrivial efficient map for that model. |

In `p64-E3-2`, the selected order is 65537; in `p64-E3-3`, it is 2. Both have
rational rank-two torsion. Their Eisenstein images are rational \(r\)-torsion,
but are outside the selected cyclic subgroup. In both cases \(X^2+X+1\)
has no root modulo \(r\), so no scalar eigenvalue exists on the selected
line. The verifier uses exact subgroup membership, not torsion membership,
to reject these cases.

The eight binary B/S records are classified specifically for the single
coordinate-squaring map supplied here. The classification does not rule out
maps over other coefficient fields or extensions.

## Timing interpretation and GLV4 boundary

Five operations are measured for every selected supported case:

| Operation | Measured scope |
| --- | --- |
| `endomorphism` | Apply the map to \(G\). |
| `scalar_decompose` | Exact rational nearest-plane scalar decomposition. |
| `scalar_mul_subgroup_baseline` | Ordinary multiplication by the same residue modulo \(r\). |
| `glv2_joint_mul_precomputed` | Signed interleaved multiplication, with components and map image precomputed. |
| `glv2_scalar_mul_full` | Decomposition, map application and signed interleaved multiplication together. |

Dataset loading and full fixture verification happen before the measured
intervals. Each operation is warmed up. Reports include individual samples,
medians, environment information and the exact scalar/components. Three
samples are a smoke benchmark; small timing differences can reflect ordinary
runtime noise. There is no constant-time, production-performance or security
claim.

These cases do not construct GLV4 by adding powers of one Frobenius map. The
separate `glv4-show`, `glv4-verify` and `glv4-bench` commands use the genuine
supplementary quadratic-extension GLV+GLS construction, whose maps are
\(1,\phi,\psi,\phi\psi\). Consult that catalog for its independent parameters,
maps and four-dimensional lattice.

The standard-library API is available as
`verify_cases(data, cases)`,
`benchmark_cases(data, cases, ids, samples, iterations)` and
`scalar_mul(curve, map_description, lattice_basis, scalar)` in
`ecbench.endomorphisms`. Callers can supply a catalog dictionary or catalog
path. The source-tree default points to the generated build artifact; the CLI
explicitly supplies its packaged catalog path.
