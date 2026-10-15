# Polynomial relation-search measurement

`ecbench.algebra.run_suite` runs each case in a fresh Python subprocess. The
optional SymPy backend provides a portable correctness baseline, with actual
`buchberger` or `f5b` computation over the declared prime coefficient field and
`lex`, `grlex`, or `grevlex` monomial order. Binary extension-field cases use
their supplied Weil-restriction equations over GF(2). No timing is predicted,
extrapolated, or substituted for an unfinished computation.

Each case has a hard wall-time budget, enforced by killing its subprocess.
The budget includes interpreter startup, SymPy import, the solve, and the
optional degree-of-regularity diagnostic. Operating-system scheduling and
process teardown can make observed wall time slightly exceed the deadline.
A solve that finishes before an optional diagnostic times out remains a
completed solve, with `diagnostic_status: timeout` and no regularity value.

These measurements solve **prebuilt fixtures**. Factor-base selection, target
construction, Semaev-polynomial construction, Weil descent, and golden-truth
enumeration took place before each backend benchmark and are outside its
recorded time. A separate fixture-generation benchmark may report preparation
cost; it does not retroactively become part of a backend solve measurement.
`fixture_generation_seconds` and
`generation_inclusive_total_pipeline_seconds` are null, explicitly marked
unmeasured for these solver runs. No contiguous generation-inclusive total
runtime can be recovered from these runs.

The main solve interval starts at prebuilt-polynomial normalization and parsing,
then includes polynomial setup, Groebner basis computation, extraction of
all roots in the declared finite coordinate domains, and independent group
verification of the resulting point decompositions. It excludes Python and
SymPy startup, parent-side input serialization, and the homogeneous-ideal
diagnostic, all distinguished in output. The suite also records observed
subprocess wall times, including this overhead. The enclosing suite wall time
includes parent serialization and orchestration over the prebuilt cases;
it excludes fixture generation and initial suite JSON loading. Ground-truth roots are read
only after the computed basis has supplied solver roots, for comparison.

Root extraction evaluates the **computed Groebner basis** over the finite
coordinate grid. This baseline does not claim a fast triangular or FGLM root
solver. Its stage timing makes the domain-enumeration cost visible. Field
and factor-base equations must be included by the case generator to match
the intended ideal; the finite domain additionally defines which roots are
eligible for relation verification. A unit ideal is a completed, empty solve.

For each computed root, the runner decodes the two factor-base x coordinates,
tries the supplied actual curve points with those coordinates, and verifies
`P + Q == target` using the independent EC reference implementation. It
records algebraic roots, roots yielding relations, ordered point
decompositions, and decompositions modulo swapping the two points. Repeated
points are allowed and have an orbit of size one. Factor-base identity points
are excluded; the current fixtures also use nonidentity targets. The reference
group operation can handle an identity target, but a fixture for that target
must define the applicable polynomial construction and domain. These policies
are recorded in each result. These exact point checks prevent Semaev
x-coordinate solutions or the opposite sign of the target from being
mistaken for verified relations. Supplied factor-base points and targets are
checked on the curve and independently checked for membership in the selected
cyclic subgroup, including the rank-two torsion distinction. Membership checks
are cached within a case, including the mathematically equivalent membership
of a point and its negative, and their cost is included in relation
verification time.

`groebner_basis.max_total_degree` means the greatest total degree in the final
basis. It is **not** a solving degree. SymPy does not expose the greatest
degree reached during its computation, so `backend.solving_degree` is null,
with an explicit reason.

The separate `degree_of_regularity` diagnostic computes the ideal generated
by the highest-total-degree homogeneous component of every input polynomial,
including finite-domain or Boolean field equations. It computes a homogeneous
Groebner basis and enumerates standard monomials in its leading ideal. The
value is the first degree at which the quotient's Hilbert function is zero.
This exact definition can differ from regularity conventions for the
homogenized ideal and from an instrumented solver's maximum working degree.
Pure-power leading monomials bound enumeration. A nonzero-dimensional
homogeneous ideal or more than one million candidate monomials returns a
documented null instead of an estimate. The Hilbert function and its
computation time accompany successful diagnostics.

Planted and independently sampled subgroup targets are summarized in separate
cohorts. Planted targets test recovery and do not estimate random-target
success probability. Unfinished solves are censored: a timeout is never
reported as an empty relation set. The primary `all_attempt_success_fraction`
is verified successful targets divided by **all attempted targets**, and
`relations_per_attempt` is verified distinct unordered signed decompositions
divided by all attempted targets. Both are lower bounds if any solve is
unfinished; the unknown yields of those cases remain censored. The companion
upper bound on success allows every unfinished target to succeed. Relations
per attempt can exceed one.

Completed-only success fractions and relations per solve second remain
available with names ending in `completed`. Complete-case rates can be biased
when hard cases time out, so they are supplementary metrics. The timeout and
error counts remain visible. Stratified summaries retain target sampling,
basis construction, and polynomial model variant.

Results include backend version, method, monomial order, Python/platform,
case hash, variable count, polynomial count, total input term count, maximum
input degree, finite-domain grid size, basis size/terms/degree, stage wall
times, exact relation counts, and ground-truth comparisons. The case hash
is a SHA-256 of compact sorted-key JSON, independent of solver settings.
Mathematical assignment values are canonical decimal strings in result JSON;
this preserves 83-bit coordinates in clients that represent JSON numbers as
IEEE 754 doubles. Internal solving and truth comparison use exact integers.
