# Measured polynomial benchmark baseline

The 69 prebuilt two-factor relation fixtures cover all seven fields. The main
55-case battery retains its original measurements, and a separate 14-case
coverage batch adds the degree-53/79 binary fields and the other three prime
fields. Each batch was run once with SymPy 1.14.0 `buchberger`, followed by
`f5b`, using `grevlex`. Cases and methods ran sequentially, in isolated
Python 3.12.14 processes on x86_64 Linux. Each case
had a 2.5-second subprocess wall budget including startup and the optional
regularity diagnostic. The recorded run completed on 2026-10-12 UTC.

The machine-readable sources are the main
[Buchberger](../build/algebra_buchberger.json) and
[F5b](../build/algebra_f5b.json) results, the coverage
[Buchberger](../build/algebra_coverage_buchberger.json) and
[F5b](../build/algebra_coverage_f5b.json) results, and the compact
[summary CSV](../build/algebra_summary.csv). Every completed solve matches its
golden finite-domain roots and exact signed point decompositions. Independent
curve arithmetic also verifies selected cyclic-subgroup membership.

| Batch | Method | Attempts | Completed solves | Censored solve timeouts | Backend errors | Additional diagnostic timeouts | Observed process wall, all attempts |
|---|---|---:|---:|---:|---:|---:|---:|
| Main | Buchberger | 55 | 53 | 2 | 0 | 1 | 45.74 s |
| Main | F5b | 55 | 51 | 4 | 0 | 10 | 60.58 s |
| Coverage | Buchberger | 14 | 14 | 0 | 0 | 0 | 4.49 s |
| Coverage | F5b | 14 | 14 | 0 | 0 | 0 | 4.60 s |
| Combined | Buchberger | 69 | 67 | 2 | 0 | 1 | 50.23 s |
| Combined | F5b | 69 | 65 | 4 | 0 | 10 | 65.18 s |

The CSV separates `main_55`, `coverage_14`, and `combined_69` rows and includes
field counts and cohort-specific denominators. The combined rows sum separate
measured batches; they do not imply a new contiguous generation/solve run.

A diagnostic timeout occurs **after a completed relation solve**. It leaves
that solve and its verified relations intact, with the diagnostic value null.
A solve timeout has an unknown relation yield and is censored; it is never
turned into a completed empty solution set.

Primary yield metrics use all attempted targets. The table counts distinct
unordered signed point decompositions; repeated signed factors are allowed
and have a swap orbit of size one.

| Scope | Method | Target cohort | Attempts | Verified successful targets | Verified success / all attempts | Verified relations / all attempts | Censored solve timeouts |
|---|---|---|---:|---:|---:|---:|---:|
| Main | Buchberger | Planted | 25 | 24 | 96% | 1.80 | 1 |
| Main | F5b | Planted | 25 | 22 | 88% | 1.64 | 3 |
| Main | Both methods | Independently seeded subgroup | 30 | 0 | 0% | 0 | 1 |
| Coverage | Both methods | Planted | 7 | 7 | 100% | 1.43 | 0 |
| Coverage | Both methods | Independently seeded subgroup | 7 | 0 | 0% | 0 | 0 |
| Combined | Buchberger | Planted | 32 | 31 | 96.875% | 1.71875 | 1 |
| Combined | F5b | Planted | 32 | 29 | 90.625% | 1.59375 | 3 |
| Combined | Both methods | Independently seeded subgroup | 37 | 0 | 0% | 0 | 1 |

Rows marked "Both methods" give identical counts separately for each method.

With unfinished solves, these verified counts divided by all attempts are
lower bounds on mathematical yield. The success upper bound allows every
unfinished target to succeed: planted recovery ranges from 96–100% for
Buchberger and 88–100% for F5b; independently seeded-target success ranges
from 0–3.33% for both methods in the main batch. The combined bounds are
96.875–100% and 90.625–100% for planted recovery, and 0–2.70% for independently
seeded-target success. All completed planted targets were recovered.
All 36 completed independently seeded targets across the two batches had verified empty relation
sets. These target scalars are deterministically sampled from a hash seed,
with reduction modulo the subgroup order; they are not exact uniform samples.

Planted recovery measures whether a solver finds existing relations. It does
not estimate random-target relation density. The mixed-cohort aggregate is
a benchmark recovery statistic and should not be interpreted as random-target
success probability. Five independently seeded cases deliberately use an
empty eligible factor base and exercise the unit-ideal/zero-yield path.

On the 51 main cases completed by **both** methods, summed core solve time was
19.14 seconds for Buchberger and 22.12 seconds for F5b. On all 65 cases
completed by both methods across the two batches, the corresponding totals
are 19.95 and 22.95 seconds. This comparison uses
the same fixtures, while the separate all-attempt wall totals also account
for startup, diagnostics, and time spent on terminated cases. There is one
timing observation per case and no repeated-run confidence interval; this
run supplies a baseline for future implementations, not a universal backend
ranking. `f5b` is SymPy's F5 in Buchberger's style, not an F4 implementation.

Among finished diagnostics, the largest reported homogeneous-top-part
regularity is 17. The largest degree in a completed final affine basis is 9.
These are different quantities. The maximum working degree reached by the
solver is unavailable and stays null. See
[the metric definitions](ALGEBRA_METRICS.md) and
[the mathematical model](ALGEBRA.md) for the exact ideal and proof boundaries.

The core solve timer includes normalization of prebuilt equations, polynomial
setup, affine Groebner computation, extraction of roots from the computed
basis over the declared finite coordinate grid, selected-subgroup checks,
and signed relation verification. The process wall timer additionally
includes interpreter/SymPy startup and regularity diagnostics. The enclosing
main-suite timers were 45.84 and 60.69 seconds, including parent serialization
and orchestration. Coverage-suite timers were 4.51 and 4.61 seconds.

Fixture generation, factor-base and target construction, Semaev construction,
Weil descent, and golden-truth enumeration were performed before these timed
backend runs. An independent
[55-fixture generation benchmark](../build/algebra_generation_benchmark.json)
measured a byte-identical rebuild at 5.04 seconds, including preparation,
golden-truth enumeration, and output serialization. That separate preparation
measurement is not added to the backend timings as a claimed contiguous
pipeline run. Generation-inclusive pipeline runtime remains null for these
solver results. The Sage/Singular adapter
is an executable export/replay path, but no Sage or Singular performance
measurement is claimed in these artifacts.
