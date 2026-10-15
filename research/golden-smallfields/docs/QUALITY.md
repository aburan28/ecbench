# Dataset quality and conformance requirements

This dataset contains 56 explicit elliptic-curve records over seven fixed finite
fields. Its purpose is reproducible arithmetic and implementation conformance at
small sizes. A matching test vector establishes the tested behavior, not a
production-security assessment or a proof that every implementation path is
correct. The records must distinguish certified integer primality, exact
point-count computations, mathematical deductions, and generated known answers.

## Representation contract

Every authoritative mathematical integer is a canonical decimal JSON string.
This includes field cardinalities, prime moduli, traces, group orders, subgroup
orders, cofactors, embedding degrees, conductors, discriminants, scalar inputs,
and prime factors. Never emit these as JSON numbers: many exceed JavaScript's
exact integer range. Bounded metadata such as degree 83, a bit length, a factor
exponent, or a schema version component may use JSON numbers.

A decimal integer uses `0` or a nonzero initial digit followed by digits, with a
single initial minus sign only where signed values are permitted. Reject `-0`,
`+1`, `01`, whitespace, decimal points, and exponent notation. Algebraic-field
coefficients represented as integer strings are the canonical representative,
not an unreduced or negative convenience expression.

Operational finite-point coordinates are `0x` followed by exactly twice the coordinate byte
length in lowercase hexadecimal. The bytes are big endian. Infinity is JSON
`null`; it is not `(0,0)`, a pair of null coordinates, or an arbitrary invalid
point. A finite point contains both coordinates and no infinity sentinel.
SEC 1 octet strings use lowercase hexadecimal **without** `0x`. These two string
formats have different roles and must not be mixed. Algorithm-selection traces
may retain decimal coordinate arrays as explanatory metadata; they do not
replace the operational point objects.

For a prime field, a canonical field element satisfies `0 <= value < p`; the
coordinate byte length is `ceil(bit_length(p)/8)`. Decoders reject values `>= p`
rather than reducing them modulo p. For a binary field of degree m, the byte
length is `ceil(m/8)` and the unsigned bit pattern must be `< 2^m`.

For binary fields, state the polynomial basis and its exact irreducible modulus.
The integer bit pattern `sum(c_i * 2^i)` represents `sum(c_i * u^i)`. Bit i is
the coefficient of u^i; integer addition is **not** field addition. Bytes still
use big-endian serialization of that integer bit pattern. All field values
needed for arithmetic, including b, j, and the curve-model discriminant, must
have an actual canonical value. Expressions such as `u^251` are explanatory
metadata, not operational field elements. The full modulus includes its
leading degree-m term. Do not silently swap a reciprocal polynomial, change
basis, or treat an extension field as integers modulo `2^m`.

## Minimal robust record structure

A schema should require the following logical fields; exact naming may follow
the published package schema. JSON Schema validates shapes and encodings, while
the reference validator enforces cross-record and arithmetic invariants.

| Record | Required content |
|---|---|
| Dataset | Schema/version, dataset identifier, provenance/version of computation tools, seven field records, 56 unique curve records, test-vector records, proof/certificate references. |
| Field | Stable identifier; kind `prime` or `binary`; characteristic; cardinality; coordinate byte length; prime modulus or binary degree, polynomial basis, and irreducible-polynomial coefficients/bit pattern. |
| Curve | Stable identifier; field identifier; exact equation/model; canonical coefficients; characteristic-dependent invariants; trace and order; complete factorization; invariant factors d1,d2; ordinary/supersingular classification; discriminant and endomorphism metadata; twist/isogeny relation references. |
| Selected subgroup | Certified prime r; h=N/r; explicit finite generator G; generator-order evidence; group exponent d2; projection multiplier d2/r; multiplicity of r in N; whether rational r-torsion is rank one; precise membership policy. |
| Vector | Stable identifier; curve or field identifier; operation; exact inputs; exact expected point/field element/bytes or expected rejection category; applicable domain and provenance. |
| Integer proof | Prime being certified, completely factored n-1, witnesses, child-certificate references, and certificate-method version. |
| Relation | Source and target identifiers, base field, relation kind, evidence; explicit map and exceptional-point behavior only when an actual map is supplied. |

A point schema is a `oneOf` between `null` and an object with exactly x,y in the
fixed coordinate format. Curve models must be explicit: `y^2=x^3+A*x+B` over
prime fields, and `y^2+x*y=x^3+a*x^2+b` over binary fields. Avoid inferring a
model from coefficient count or reusing prime-field addition formulas in
characteristic two. Reject duplicate identifiers and unresolved references.

A subgroup generator is a generator of the named cyclic subgroup, not a
claim to generate the whole rational point group. Any floating-point estimates
or bit-security summaries are derived annotations, never authoritative input.

## Mathematical invariants and proof boundaries

For every curve verify:

- `N=q+1-t`, `t^2 <= 4*q`, and that the model is nonsingular.
- The complete order factorization multiplies to N; every reported factor is
  prime with a checked certificate.
- `d1*d2=N`, `d1|d2`, and `d1|(q-1)`. The group exponent is d2, not necessarily N.
- The characteristic does not divide t exactly for the ordinary records.
- The chosen r is the declared certified prime factor, h=N/r, and the named
  finite G is on the curve, nonzero, and satisfies `[r]G=O`. These conditions
  establish the exact order r because r is prime.
- If r differs from the characteristic, `q^k=1 mod r` and
  `q^(k/ell)!=1 mod r` for every prime ell dividing k. Those tests establish the
  exact pairing embedding degree. This is not automatically the extension
  degree defining full r-torsion, especially at k=1.

Keep these discriminants distinct:

| Quantity | Domain and meaning |
|---|---|
| Curve-model Delta | A field element. For prime models, `-16*(4*A^3+27*B^2)`; for the stated binary model, b. Nonzero means nonsingular. It changes by a twelfth power under a coordinate change. |
| Frobenius D_pi | The signed integer `t^2-4*q`, the discriminant of `Z[pi]` when pi generates a quadratic order. |
| Fundamental D_K | The fundamental discriminant of the quadratic Frobenius field, with `D_pi=f_pi^2*D_K`. |
| Ordinary D_End | The actual endomorphism-order discriminant `f_E^2*D_K`, where `f_E|f_pi`. A trace alone generally does not determine f_E. |
| Supersingular rational-order discriminant | A quadratic-order discriminant for `End_Fp(E)` in the t=0 prime-field examples. It must be labeled rational, not geometric. The geometric ring is a maximal quaternion order; its basis/discriminant is not supplied by this quadratic number. |

The supersingular distinction between `Z[sqrt(-p)]` and the maximal quadratic
order for p=3 modulo 4 depends on whether **all of E[2]** is rational, not on
whether merely one point of order two is rational. For ordinary records,
maximality requires either fundamental D_pi, a known maximal-CM model, or a
separate torsion/isogeny computation. Preserve the proof method and evidence.
The core verifier checks discriminant identities and allowed conductors; actual
model-specific endomorphism-ring identifications are source-proven metadata.
It does not repeat the degree-480 31-torsion calculation or every special-CM
ring argument on each run. The original mathematical computation programs in
provenance make those arguments reproducible. Known small class numbers are
checked; the large Koblitz Frobenius-order class numbers remain exact
source-proven values whose formula is not recomputed by the core verifier.

Prime certificates prove integer primality. Rabin tests prove irreducibility
of the binary reduction polynomials. PARI/GP point counts and group computations
are exact computations whose reproducible commands and tool version are part
of the evidence. Special-CM counts additionally have exhaustive trace-candidate
arguments and rejection witnesses; Koblitz counts have a base-field Frobenius
recurrence. Testing `[N]P=O` on a handful of points is a consistency check and
must not be presented as a point-count or full-group-structure proof. Complete
benchmark domain parameters include a chosen prime-subgroup generator; this
does not imply that a full rational-group basis or an explicit supersingular
geometric quaternion-order basis has been supplied. State those exclusions.
Independent r-torsion H on the three rank-two examples covers subgroup
selection, but is not a full-group basis for those higher-exponent groups.

Known-answer vectors must be independently reproduced, preferably using a
separate implementation/tool where available. A vector generator and validator
that share arithmetic code can certify internal consistency while sharing the
same bug. State this boundary explicitly rather than promising independent
verification that was not performed.

## Subgroup projection and membership traps

Use the group exponent to project a trial point:

`G = [d2/r]P`.

Since `[d2]P=O`, any nonzero result has order r. Search for a nonzero result and
verify its exact order; a zero result is not a generator. Using `[N/r]P` is not
an interchangeable rule for noncyclic groups.

For `p61-E4-0`, `N=2^61` and the group is `Z/2 x Z/2^60`, so the exponent is
`2^60`. The usual `N/r=2^60` with r=2 kills every point. The exponent-based
multiplier `d2/r=2^59` can produce a valid order-two generator.

This same curve has rational two-torsion of rank two. All three nonzero
order-two points satisfy `[2]P=O`, but only one lies in the chosen cyclic
subgroup `<G>`. Thus `[r]P=O` is an r-torsion test, not a universal proof of
membership in the chosen `<G>`. There are three rank-two selected-r cases in this dataset:

| Curve | Selected r | Rational group |
|---|---:|---|
| p61-E4-0 | 2 | Z/2 x Z/2^60 |
| p64-E3-2 | 65537 | Z/4294967295 x Z/4294967295 |
| p64-E3-3 | 2 | Z/2^32 x Z/2^32 |

For r=2, membership in the chosen subgroup is exactly `P=O` or `P=G`. For
r=65537, a complete baby-step/giant-step membership test in `<G>` is practical
and can establish the distinction exactly. Include another independent
r-torsion point H as a negative chosen-subgroup membership vector on each of
these three records: H is on-curve, H!=O, `[r]H=O`, and H is not in `<G>`.
Where r does not divide d1, rational r-torsion is the unique cyclic subgroup of
order r, and on-curve validation plus `[r]P=O` suffices. A protocol can separately reject
infinity even when a mathematical membership predicate accepts it.

Curve-level scalar multiplication must not silently reduce a scalar modulo r
for arbitrary on-curve input points. Such reduction is valid only within the
specified r-subgroup. Include a point outside that subgroup whenever one
exists; an odd-cofactor point or the binary order-two point is useful. A
prime-order full rational group has no such point, so mark that coverage case
inapplicable rather than fabricating it.

## SEC 1 encoding and decoding

The supported strict formats are infinity `00`, compressed `02/03 || X`, and
uncompressed `04 || X || Y`. Infinity is exactly one byte, compressed length is
`1+L`, and uncompressed length is `1+2L` for coordinate length L. Reject other
prefixes, trailing bytes, short/long inputs, odd-length hex, malformed hex,
noncanonical field values, and off-curve finite points. If hybrid 06/07 forms
are unsupported, document that choice and reject them explicitly.

For prime fields, the compressed tag selects the low bit of y. Decode by
solving `y^2=x^3+A*x+B`, reject a nonsquare right side, and select the root with
the requested low bit. If the only root is y=0, tag 03 must be rejected: there
is no odd root. Do not turn that root into the noncanonical integer p.

For binary fields with x nonzero, the tag selects the low bit of `z=y/x`, not
the low bit of y. Recover z from

`z^2+z = x+a+b/x^2`.

The right side must have absolute trace zero. Its two solutions differ by 1,
so the tag selects exactly one; then y=x*z. All three degrees in this dataset
are odd, so a correctly implemented half-trace can solve the equation. Test
both tags and reject trace-one inputs. If x=0, there is one possible ordinate,
`y=sqrt(b)=b^(2^(m-1))`, and SEC 1 requires tag 02; tag 03 is invalid even when
the actual y bit pattern is odd. Include this exception for every binary model.

Binary negation is `(x,x+y)`. For nonzero x, negation toggles the compressed
tag because `(x+y)/x=y/x+1`. Include a point with low bit of x equal to zero:
its y and x+y have the same low bit, so a codec incorrectly tagging y can be
exposed. For x=0, the unique point is self-inverse and has order two.

Uncompressed inputs also require on-curve validation; syntactically valid
coordinates are not automatically a valid point. A decoder should produce
one canonical point object and an encoder should produce one canonical octet
string. Round trips alone are insufficient: two equally wrong implementations
can round-trip their own incompatible convention.

## Conformance coverage matrix

"Every curve" means all 56 records unless the stated feature is mathematically
inapplicable. Record inapplicability explicitly. Tests use exact expected
values, not just expected equality between two outputs from the same code.

| Area | Required cases | Coverage/expected behavior |
|---|---|---|
| Field arithmetic | zero, one, addition/subtraction, multiplication, squaring, inverse, distributivity | All seven fields; binary addition is XOR and reduction uses the exact polynomial. |
| Field rejection | inverse of zero, out-of-range/noncanonical input, invalid string syntax | Deterministic rejection; no silent reduction by decoder. |
| Curve model | coefficients, Delta, j, ordinary/supersingular | Every curve; independently match metadata. |
| Point membership | selected G, general valid P, invalid finite point | Every curve; coordinate bounds and equation both checked. |
| Infinity | O+P, P+O, O+O, -O, [0]P, [k]O | Every curve; expected infinity uses JSON null. |
| Negation | -P, P+(-P), -(-P) | Every curve; characteristic-two formula explicitly covered. |
| Addition | distinct points P,Q with nonexceptional denominator | Every curve; independently generated expected sum. |
| Doubling | P+P and [2]P; order-two doubling | Every curve; order-two case where present; no divide-by-zero crash. |
| Scalar multiplication | 0,1,2,3,r-1,r,r+1 and selected large scalar | Every selected subgroup; exact expected points. |
| Scalar domain | negative scalar support or explicit rejection; no implicit mod-r for general P | Apply documented policy; include outside-subgroup point where present. |
| Generator/order | G!=O, on curve, [r]G=O | Every selected subgroup, checked certified prime r. |
| Projection | [d2/r]P yields the published nonzero G | Every curve; explicit M61 noncyclic regression against N/r. |
| Chosen-subgroup membership | valid G, infinity policy, outside subgroup | Include independent r-torsion H on each of the three rank-two selected-r records. |
| Group metadata | factor product, invariant factors, exponent; supplied R of exact order d2 | Every curve; prove R order by [d2]R=O and [d2/ell]R!=O for every prime ell dividing d2; R generates the full group only when d1=1. |
| Prime SEC 1 | exact compressed/uncompressed/infinity bytes; both signs | Every prime curve; y=0 exception where present. |
| Binary SEC 1 | exact bytes; both y/x tags; x=0; y-low-bit trap | Every binary curve; unique x=0 ordinate encoded with tag 02. |
| SEC 1 rejection | invalid prefix, lengths, hex, bounds, nonsquare/trace-one x, off-curve uncompressed | Both field types; invalid cases belong to stated curve. |
| Twist metadata | same j, opposite trace, N+Ntwist=2(q+1), equal D_pi | Every declared quadratic pair, including documented self-isomorphic special-j twist. |
| Explicit twist map | map exists over quadratic extension | Do not invent a base-field coordinate map from metadata; extension-field vectors only if implemented. |
| Isogeny class | same field and same trace | All cross-record relations; different j can share class, same j can have opposite traces. |
| Explicit 2-isogeny | phi(x,y)=(x-1/x,y*(1+1/x^2)); target y^2=x^3+4x; kernel O,(0,0) | Supplied map vectors cover generic P, kernel, infinity, and on-target property. Target is an auxiliary model, not an extra primary record. |
| Gaussian automorphism | (x,y)->(-x,i*y), i^2=-1 | Ordinary prime-field j=1728 models; validate parameter and expected image; phi^2=-1 and phi^4=1 give useful composition checks. |
| Eisenstein automorphism | (x,y)->(zeta*x,y), zeta!=1, zeta^3=1 | Ordinary prime-field j=0 models; validate parameter and image; phi^3=1 and 1+phi+phi^2=0 give useful composition checks. |
| Binary relative Frobenius | (x,y)->(x^2,y^2), b->b^2 for a=0,1 | Supplied map vectors cover all 14 binary records; target curve changes unless b lies in F2. |
| Embedding degree | modular equality and prime-factor minimality | Every selected r!=characteristic; k distinct from field degree and subgroup size. |
| Integer/field proofs | recursive Lucas certificates; Rabin irreducibility | All declared certificate subjects and all binary field moduli. |
| Integrity | unique IDs, references, record/vector counts, deterministic hashes | Whole bundle; stable regeneration and independent checks. |

An isogeny-class label proves class membership from the fixed-field trace
criterion. It does not supply an explicit isogeny or a shortest path. An
explicit map has a named source and target equation, including its exceptional
point behavior; equality of counts alone must not be described as a tested map.
Twists with opposite nonzero traces are not isogenous over that base field.
A twist partner identifier can name an Fq-isomorphic representative of the
standard twist, rather than the literal coefficients `(d^2*A,d^3*B)`. This is
necessary when special-j classes are indexed modulo fourth/sixth powers. For
example, the d=37 twist of `p61-E4-0` has A=-37^2; its partner representative
has A=-1 because these models are Fp-isomorphic. Do not apply a base-field
coordinate map merely by equating partner coefficients, and document this
representative convention in the schema/package.
Special-j self-isomorphic quadratic twists must not be forced into a generic
"two distinct curves" rule.

## Release acceptance

The release should pass its schema validation and reference verifier, all
published positive and negative known-answer vectors, every certificate check,
and field/curve/subgroup cross-record identities. Independently reproduced
vectors and the actual checks run must be listed in release evidence. Repeated
execution should produce the same canonical JSON and hashes. A deliberate
mutation of a coordinate, expected sum, subgroup order, SEC 1 tag, field modulus,
and certificate witness should cause a deterministic failed check.

The 53-83-bit sizes remain experimental: large embedding degree or CM
fundamental discriminant does not make these curves comparable to modern
128-bit-security parameters. No timing, constant-time behavior, side-channel
resistance, or protocol safety claim follows from functional conformance.

## Structural schema and explicit map vectors

The Draft 2020-12 schema is `schemas/dataset.schema.json`. It constrains
canonical strings, object shapes, finite points/null infinity, exclusive
expected-result tags, proof/factorization structures, field/curve records, and
each vector domain/operation. Field-dependent byte widths, numerical ranges,
references, counts, hashes, proofs, and arithmetic belong to the semantic
verifier. It supports the `map` domain with `degree2_isogeny`,
`gaussian_automorphism`, `eisenstein_automorphism`, and `relative_frobenius`.
A map vector supplies P, exact target coefficients, and a field-element
parameter when the automorphism requires one. Targets are auxiliary models;
they do not change the primary corpus of 56 curves. Expected images were
supplied independently, and the core checks the source, defining parameter,
exact target coefficients, image membership, and the expected image.

Operational data and canonical certificates use exact string encodings.
The `provenance` directory can also contain legacy original JSON snapshots
with integer JSON numbers; these require an exact-integer JSON parser. They
are historical reproduction evidence, not the canonical interface for
JavaScript numeric parsing.

## Algebra fixture and measured-result contracts

`schemas/algebra_cases.schema.json` constrains the exact two-summand sparse
systems, finite variable domains, factor-base point pools, target derivations,
and golden root/relation records. It accepts both the original 55-case catalog
and the separate 14-case field-coverage catalog, without hard-coding either
count. `schemas/algebra_results.schema.json` constrains completed, timed-out,
and error results, backend identity, unavailable intermediate solving degree,
finite-domain root encodings, nonnegative stage times, and all-attempt
aggregates. Large field values and root assignments remain decimal strings.
These contracts allow bounded numerical counts and floating-point runtime
measurements; they do not allow large mathematical integers as JSON numbers.

The schemas reject additional properties on authoritative records. They do
not prove that exponent-vector lengths match the variables, field coefficients
are canonical residues, factor-base bases are independent, sparse equations
are equivalent to the stated Semaev system, or the published roots are
complete. They also do not certify a returned Gröbner basis or the truth of a
timing claim. Those require semantic checks and the recorded proof boundaries.
The original-system fixture verifier exhaustively evaluates all 36,528
declared assignments across the 69 cases, then checks rational lifts, signed
point sums, and membership in the chosen cyclic subgroup. This certifies those
finite-system golden answers; it does not certify ideal equality for an
external backend's returned basis.

Release schema checks pass both catalogs and four measured result collections
(138 result rows), alongside the primary 6,931-vector dataset. Twenty deliberate
algebra mutations are rejected, including noncanonical/numeric coefficients,
negative or fractional exponents, malformed points, additional properties,
invalid statuses, fractions exceeding one, and replacing an unavailable
algorithm solving degree with the final basis degree. Evidence is recorded in
`build/algebra_schema_verification.json`.
