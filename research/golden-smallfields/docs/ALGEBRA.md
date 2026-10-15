# Algebraic elliptic-curve benchmark methodology

The primary algebraic battery supplies exact finite-field polynomial systems
and measures Gröbner computation, finite-domain root extraction, point lifting,
and relation verification. The recorded SymPy target timings start from prebuilt
fixtures; system construction is excluded. A separate generation measurement
records preparation of the original 55-case catalog, including its golden
truth enumeration. These measurements do not provide one generation-inclusive
per-target solve time. Scalar multiplication, isogenies, Frobenius, and GLV/GLS
decompositions are separate benchmark families; their timings do not substitute
for relation yield or polynomial-solver work.

A small factor base makes reproducible experiments possible on the existing
53–83-bit fields. Such experiments do not by themselves estimate a complete
ECDLP attack: relation collection, matrix rank, linear algebra modulo the
subgroup order, and individual-log recovery must also be implemented and timed
before making that claim.

## Exact systems and Semaev equations

For a prime field of characteristic greater than three and
`E: y^2=x^3+A*x+B`, use

```
f3(x1,x2,x3) = (x1-x2)^2*x3^2
              -2*((x1+x2)*(x1*x2+A)+2*B)*x3
              +(x1*x2-A)^2-4*B*(x1+x2).
```

For the ordinary binary model `y^2+x*y=x^3+a*x^2+b`, use

```
f3(x1,x2,x3) = (x1*x2+x1*x3+x2*x3)^2+x1*x2*x3+b.
```

The binary formula is independent of a. This does not mean that its two
quadratic twists have the same rational point group or relation yield: their
point lifts and signed sums still differ. For a two-summand relation to a
specified finite target Q, substitute `x3=x(Q)` **before** recording polynomial
degrees or constructing the solver input.

The vanishing of f3 describes an x-coordinate condition over the algebraic
closure. It does not alone prove that the proposed x coordinates lift to points
over the required field, that the points belong to the chosen factor base or
subgroup, or that their signed sum is the specified Q rather than -Q. All those
conditions require a separate verification stage. For example, on
`y^2=x^3-x` over F7 with target `(0,0)`, `f3(2,3,0)=0`; neither x=2 nor x=3 has
a rational ordinate. Treating that polynomial root as a relation is incorrect.

Record the exact affine generators, variable order, coefficient field,
monomial order, normalization/preprocessing, factor-base equations, and
field equations. Hash their canonical sparse representation. A sparse term
has a canonical coefficient and an exponent vector of the declared variable
length. Coefficients and large moduli are decimal strings; exponents and small
structural metadata can be bounded JSON integers. A solver must never replace
an 83-bit prime modulus with a nearby machine-word modulus.

## Binary Weil descent and Boolean normalization

For a binary field, preserve the dataset's exact polynomial basis and
irreducible modulus. Write each unknown factor-base x coordinate as an F2
linear combination of explicitly listed basis elements. If the factor-base
subspace has dimension d, two unknown x coordinates introduce 2d Boolean
variables. Expanding f3 in the recorded field basis produces m coefficient
equations over F2. The equality is an equality of all m coordinates, not just
one integer bit or one trace equation.

Perform coefficient-field reduction modulo the irreducible binary polynomial.
Boolean normalization may replace each variable power v^e, e>0, by v in the
Semaev-derived equations, because the affine system also includes `v^2-v=0`.
Keep those Boolean field equations explicitly: Boolean-normalizing the field
equations themselves to zero and dropping them loses both the finite-domain
constraint and the stated degree-of-regularity definition.

With a fixed target and two unknown x coordinates, the binary f3 equations
become degree at most two after this declared Boolean normalization. Squaring
in the extension field is F2-linear; it is not generic integer squaring of a
coordinate bitmask. State whether zero equations, constant equations, linear
elimination, or variable substitution have been removed or performed. Those
choices affect the exact generator list and can affect the regularity metric.

The polynomial coefficient field for this descent is F2. The group-relation
coefficient modulus is the selected subgroup order r, usually an odd prime.
These are different moduli. In particular, `-P=(x,x+y)` is usually different
from P, and a group-relation coefficient -1 must not be changed to +1 merely
because the curve's field has characteristic two.

## Factor bases and target populations

A case must identify its factor-base **x domain** and its admissible **point
pool**, including whether both signs are retained or one point per x is chosen
with a recorded sign convention. A domain of B x values is not automatically a
factor base of 2B distinct points: some x values have no rational lift, and
order-two points have only one lift. For a binary subspace, report both its
cardinality 2^d and the number of admissible rational point lifts.

Prime-field membership can be encoded by `product(x-a for a in X_base)=0`.
Its roots are exactly the listed, distinct base-field x values. If this
polynomial splits over Fp, the separate field equation `x^p-x=0` is redundant
for that variable; record whether it is included. For a binary subspace,
parameterizing x by independent Boolean coordinates restricts membership
without introducing a large univariate membership polynomial.

Planted cases and random cases answer different questions and must have
separate result tables and aggregates:

- A planted target is constructed as Q=P1+P2 from admissible factor-base
  points. Retain the planted point pair and verify that it satisfies every
  recorded affine equation and the declared duplicate/sign policy. Such a
  case has a known relation; it tests solver recovery, not random-target yield.
  Sampling planted sums biases the target distribution toward targets with
  more decompositions.
- A random target must state its population and sampling procedure. Uniform
  nonzero scalars times the chosen subgroup generator give uniform nonzero
  subgroup targets. A point drawn from x-coordinate scans is not a uniform
  random curve point. Multiples of the maximal-exponent point R sample only
  `<R>`; on noncyclic records they do not sample the whole rational group.

The factor-base points must also lie in the selected cyclic subgroup if the
results are described as ECDLP subgroup relations. Otherwise label them
full-curve relation experiments. Projecting a point by a cofactor/exponent
multiplier changes its x coordinate and can remove it from the declared
factor base; membership has to be checked again. The dataset's rank-two
selected-r cases require exact membership in the chosen `<G>`, not only an
r-annihilation test.

State whether infinity, the coordinate x=0, order-two points, repeated x
coordinates, repeated signed points, inverse pairs, and duplicate targets are
permitted. These are separate restrictions. Infinity has no finite x
coordinate and therefore cannot simply be inserted into the displayed f3.
An empty admissible target population must be reported as such, not handled
by an infinite rejection-sampling loop. Record input draws and policy
rejections separately from actual target-solving attempts.

## Relation yield and multiplicities

Record at least these counts for each attempted target:

| Count | Exact interpretation |
|---|---|
| Polynomial roots | Ordered finite-domain assignments satisfying the polynomial system. Boolean coordinate representations must be unique. |
| Verified x tuples | Ordered x tuples for which at least one admissible rational point lift yields the requested signed sum. |
| Ordered signed decompositions | Distinct ordered point pairs (P1,P2) satisfying P1+P2=Q and the factor-base policy. |
| Unordered signed decompositions | The preceding point pairs modulo swapping the factors. |
| Unique relation rows | Distinct normalized factor-base coefficient rows relative to the specified target; repeated factors have combined coefficients. |

Do not divide a root count by two or four without checking orbit sizes.
A repeated point has a swap orbit of size one; an order-two point has one sign
rather than two. Distinct signed pairs have swap orbits of size two. Repeated
x coordinates do not always imply repeated signed points. If duplicate factors
are forbidden, enforce that policy before counting relations and retain the
number excluded. Coefficient rows over a subgroup are interpreted modulo r;
a sign-quotiented x factor base needs an explicit representative for each x.

`any_relation_success_rate` is the number of attempted targets with at least
one verified relation divided by **all attempted targets**. A separate
`complete_solve_success_rate` requires both a completed solve and a verified
relation. `completion_rate` counts completed solves, including completed cases
with no relation. `relations_per_attempt` is the total number of verified,
uniquely counted relations divided by all attempted targets; it is a yield and
can exceed one. Label its counting convention.

Timeouts, solver errors, and attempted cases yielding no relation remain in
the denominator. Do not report a conditional rate over completed successful
runs as the primary success rate. A timed-out case with no verified relation
has an unknown mathematical yield, not a proof of zero relations. Partial
verified relations can be reported as a lower bound with a partial status.
Use null, not NaN or infinity, when a denominator is zero. If relation storage
is capped, record the cap and that the count may be a lower bound.

## Degree of regularity: the measured definition

Let F be the exact, declared affine generator list after the recorded
preprocessing, including factor-base and field equations. For each nonzero
f in F, take its highest-total-degree homogeneous component `top(f)`. Define

```
S = coefficient_field[v1,...,vn]
J = <top(f) : f in F, f != 0>.
```

This benchmark's `degree_of_regularity` is the first degree d at which the
Hilbert function `dim(S/J)_d` is zero, **provided J is Artinian**. For a proper
Artinian J, compute a Gröbner basis of J and its standard monomials. Then

```
degree_of_regularity = 1 + max(total_degree(M) for standard monomials M).
```

If J is the unit ideal, the Hilbert function is already zero in degree zero,
so the value is 0. That condition concerns J: an inconsistent affine ideal
I=<F> with Gröbner basis [1] can still have a proper, nonunit highest-parts
ideal J. Do not assign regularity zero merely because the affine solve finds
no roots. If J is not Artinian, report the finite metric as undefined with a
reason, rather than guessing a bound.

The listed Boolean field equations contribute `v_i^2` to J, so every standard
monomial is squarefree and the quotient is finite. Prime factor-base
membership equations of degree B contribute `x_i^B`. These pure powers give
explicit exponent bounds for exact standard-monomial enumeration at small n.
Compute the Hilbert distribution by rejecting monomials divisible by a leading
monomial of the homogeneous Gröbner basis. Preserve the distribution, quotient
dimension, bound used, and a completion status. If the calculation times out,
its exact degree is unavailable; a bound is a separately labeled quantity.

For the two-variable prime f3 system with target x fixed and equal factor-base
cardinality B>=2, the highest part of f3 is `x1^2*x2^2`, while the membership
polynomials contribute `x1^B,x2^B`. In this basic formulation the largest
standard-monomial degree is B, and the defined regularity is B+1. It therefore
has little curve/target variation. This is a useful correctness control and
also a reason not to treat the number as measured attack complexity.

This definition is tied to the exact generator/preprocessing list. Different
affine generating sets for the same ideal can have different highest-parts
ideals. Record both the original and preprocessed generator hashes when using
elimination, substitution, or other transformations.

Keep three distinct quantities:

| Metric | Meaning |
|---|---|
| `degree_of_regularity` | The exact Hilbert-function quantity for J defined above. |
| `max_final_groebner_degree` | Largest total degree in the final affine Gröbner basis; [1] has maximum degree zero by the stated convention. |
| `algorithm_peak_degree` / a documented solving-degree trace | Largest degree processed by the actual algorithm, with an explicit instrumentation definition covering intermediate polynomials or Macaulay matrices. |

Neither regularity nor the maximum degree of the final basis is automatically
the actual solving degree. SymPy's ordinary public Gröbner API does not expose
its complete intermediate degree trace. If it is not instrumented, report that
metric as unavailable; do not relabel the final basis maximum as solving degree.

## Runtime, completion, and backend identity

Measure monotonic wall time around the complete target pipeline, including
case generation, polynomial construction/descent, normalization, affine
Gröbner computation, finite-field root extraction, point lifting, and relation
verification. Retain stage times and the enclosing total. Account for
unattributed orchestration/serialization time; stage sums need not equal the
outer measurement exactly. If regularity analysis is an additional stage,
report its time separately and state whether the displayed total includes it.
A core relation-pipeline time can exclude that diagnostic work only if an
additional whole-case total includes it and the distinction is explicit.

Record solver/backend name, exact algorithm, version, coefficient domain,
variable and monomial order, worker/thread count, seed, hardware/runtime
metadata, repetitions, and timeout budget. SymPy `method='buchberger'` and
`method='f5b'` are actual different algorithms. The latter is F5 in Buchberger's
style; it should not be described as an F4 matrix implementation or assumed to
match a separately implemented F5 solver.

A meaningful timeout must interrupt the actual worker, not merely test elapsed
time after a blocking solve returns. The deadline must cover whichever total
pipeline the configuration claims. Retain the timeout status, configured
budget, elapsed parent wall time, stages completed, and partial evidence if
available. A terminated computation is not a completed basis/regularity solve.
Use distinct statuses for completed/no-relation, completed/relations-found,
timeout, backend error, unsupported capability, and invalid case.

Completed-run medians and distributions are conditional statistics. Publish
the completion fraction and full attempted-case table alongside them. A timeout
is censored data, so do not treat its budget as an exact solve duration when
computing a median or as an omitted observation when claiming speedups.
Repeated targets or seeds should be identifiable. A planted-recovery rate and
a random-target relation rate must not share one unqualified aggregate.

Computing a Gröbner basis and then exhaustively enumerating the small finite
domain is a valid measured hybrid pipeline if named accurately. Record whether
roots are obtained from the basis, from original equations, or by an
independent exhaustive check. Characteristic-zero SymPy root routines must not
be used silently to solve the F2 or Fp system.

## Optional Sage/Singular adapter contract

An adapter receives a canonical system descriptor: coefficient field, ordered
variables, monomial order, sparse affine generators, domain equations,
preprocessing metadata, seed/case identifier, and timeout request. It returns
its actual backend/algorithm/version, completion status, canonical basis,
finite-domain roots or its root-extraction method, stage timing, and any
available degree trace. A different system/order requires a different case
identity, even if two runs share the same target curve.

Map order names explicitly (`grevlex` versus Singular `dp`/Sage `degrevlex`).
Reject unsupported large-prime coefficient domains instead of truncating the
modulus or moving to characteristic zero. The binary descended systems have
coefficient field F2, while the original extension-field basis remains part of
the system-construction provenance. An adapter is optional until actually
installed and executed; do not include an unrun backend in measured results.

Independently checking that input generators reduce to zero modulo a returned
basis G proves `I subset <G>`. Buchberger S-pair checks prove G is a Gröbner
basis of its own ideal. Together these still do not prove `<G>=I`: G={1}
passes both checks for any I. To claim an independently certified basis, also
prove every returned basis polynomial belongs to I through combination
witnesses or a trusted independent computation. Otherwise label it
backend-computed and independently verify candidate roots and actual curve
relations against the original case.

## Endomorphism families and a genuine GLV/GLS four-way construction

A two-way decomposition uses a verified endomorphism eigenvalue lambda on the
specified prime-order cyclic subgroup, not a guessed root of an endomorphism
polynomial. Test the point identity `phi(G)=[lambda]G`. The relation lattice is

```
L2 = {(a,b) in Z^2 : a+b*lambda = 0 modulo r}.
```

For four-way GLV/GLS, take an ordinary j=0 curve
`E0/Fp: y^2=x^3+B`, with p=1 modulo 3. Use an explicitly represented quadratic
extension K=Fp^2, an independently checked nonsquare eta in K, and the twist

```
E_eta/K: y^2=x^3+B*eta^3.
```

Let zeta be a nontrivial cube root of unity in Fp. Put
`alpha=eta^((1-p)/2)` in K and define

```
phi(x,y) = (zeta*x,y)
psi(x,y) = (alpha^2*x^p, alpha^3*y^p).
```

Both are endomorphisms of this twist, and they commute because zeta is fixed
by p-power Frobenius. Phi has degree one and order three. Psi has inseparable
degree p. Globally, as a geometric endomorphism,

```
psi^2 = -pi_(p^2),
```

where pi_(p^2) is the degree-p^2 Frobenius. On K-rational points that Frobenius
acts as the identity, so `psi^2(P)=-P`. It is incorrect to assert
`psi^2=[-1]` as a global geometric-map identity: the degrees disagree.

If t is E0's trace over Fp, the twist has exact K-rational order
`(p-1)^2+t^2`. Psi also satisfies `psi^2-t*psi+[p]=0`. On a stable prime-order
subgroup where t is invertible modulo r, its eigenvalue is therefore
`mu=(p-1)/t modulo r`. Verify `psi(G)=[mu]G` and `mu^2=-1 modulo r`. Choose the
correct root lambda of `lambda^2+lambda+1=0` by checking phi(G), not merely by
selecting the smaller root. If a subgroup is not unique, establish its
stability under both maps and verify the eigenvalue identities explicitly.

The four scalar images are `G,phi(G),psi(G),phi(psi(G))`, with eigenvalues
`1,lambda,mu,lambda*mu`. The actual four-dimensional relation lattice is

```
L4 = {a in Z^4 : a0+a1*lambda+a2*mu+a3*lambda*mu = 0 modulo r}.
```

A full-rank basis must have determinant of absolute value r, and each row
must satisfy the congruence. Reduced-basis/Babai decomposition must reconstruct
the requested scalar modulo r, and point reconstruction must match `[k]G`.
Record the basis, reduction method, coefficient sizes, maps, parameters,
eigenvalue identities, and independent test vectors. The four images are not
four linearly independent points over Z/rZ: the subgroup is cyclic. Nor does
this construction assert a quartic geometric endomorphism algebra; ordinary
elliptic-curve endomorphism algebras remain quadratic.

A generic LLL calculation alone does not prove a uniform r^(1/4) coefficient
bound. Record the actual coefficient lengths and any construction-specific
bound used. Timing a four-term reconstruction with four independent ordinary
scalar multiplications need not accelerate anything; benchmark the actual
joint multiplication algorithm and include decomposition and endomorphism
application costs.

For binary Koblitz models, the F2 Frobenius tau satisfies
`tau^2-mu0*tau+2=0`, with mu0=2a-1. Its action and two-way decomposition are
useful separate cases. Simply using `1,tau,tau^2,tau^3` does not create a
balanced four-way GLV/GLS construction: the displayed quadratic relation gives
constant-size integer syzygies, leaving only two effective scalar directions.
For general binary coefficients, F2 Frobenius maps to a conjugate curve and is
not automatically an endomorphism of the same curve. Label those relative
isogeny cases separately.

## Sources and scope

- Igor Semaev, *Summation polynomials and the discrete logarithm problem on
  elliptic curves*, IACR ePrint 2004/031: https://eprint.iacr.org/2004/031 .
- SymPy polynomial reference and internal Gröbner algorithms:
  https://docs.sympy.org/latest/modules/polys/reference.html and
  https://docs.sympy.org/latest/modules/polys/internals.html .
- Yao Sun and Dingkang Wang, *The F5 Algorithm in Buchberger's Style*:
  https://arxiv.org/abs/1006.5299 .
- Patrick Longa and Francesco Sica, *Four-Dimensional Gallant-Lambert-Vanstone
  Scalar Multiplication*, full version: https://eprint.iacr.org/2011/608 .

These sources motivate the constructions. Exact case equations, proof
boundaries, measured backend executions, and independently verified relation
and endomorphism vectors remain the evidence for this benchmark release.
