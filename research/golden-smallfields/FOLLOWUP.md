# Open-lead follow-up — 15 October 2026

The [updated 28-page PDF](docs/ecbench_updated_report.pdf) appends six pages of
new research to the original report. The [six-page follow-up](docs/ecbench_open_leads_followup.pdf)
is also available separately. The original release and its historical measurements
are unchanged; all 108 files in its checksum inventory still match.

## Findings

1. **Exact preprocessing improves the difficult binary systems.** GF(2) row
   reduction exposes hidden linear constraints. All 57 binary fixtures have
   two-way row-span certificates proving equality of the original and derived
   polynomial ideals, plus exhaustive checks of their original finite domains.
   Variables, factor bases, targets, and signed relation policies are unchanged.
2. **The improvement is measured.** Four difficult cases × two methods × two
   variants × three repetitions produced 48 fresh runs. All 24 transformed runs
   completed; 21 of 24 original runs completed. For the seven case/method
   combinations with matched completed pairs, median paired wall speedups range
   from 1.38× to 9.86×. Wall time includes preprocessing, witness verification,
   startup, serialization, solving, and relation checks. The eighth combination
   has three original timeouts and three completed transformed runs; no finite
   speedup is inferred for it. These are local measurements on selected cases.
3. **Regularity is generator dependent.** In one exact row-equivalent pair,
   the highest-input-part regularity changes from 7 to 0; another changes from
   7 to 6. The affine polynomial ideal is identical. This diagnostic cannot
   classify curves independently of the recorded generators and preprocessing.
4. **The empty sampled-target results have very little statistical power.**
   An exact census of all signed pair sums in 27 distinct factor bases agrees
   with every golden relation count. Under a hypothetical uniform nonzero
   target model, the 37 seeded-case contexts yield only 1.80005e-6 expected
   successful targets in total. The existing targets are deterministic hash
   draws, so this number is a model calculation rather than a p-value.
5. **A represented isogeny path is now explicit.** A degree-two map and its
   dual connect `p61-E4-0` and `p61-E4-1`, with the exact coordinate scale for the
   second stored model. Their groups are `Z/2 × Z/2^60` and `Z/2^61` despite
   equal order. Symbolic curve and dual-composition identities, 26 point vectors,
   and 340 homomorphism checks pass.
6. **All 69 original systems now have complete final-basis metrics and exact
   highest-input-part Hilbert diagnostics.** The consolidated catalog links each
   value to an actual completed computation. No method-specific timing is
   borrowed to fill another method's timeout. Algorithm maximum working degree
   remains unavailable.

## Reproduce

Use the source checkout, with the optional algebra/report dependencies installed:

```sh
python -m pip install '.[algebra,report]'
python tools/follow_up_leads.py --phase all --timeout 20 --repetitions 3
python tools/verify_followup.py
python -m unittest discover -s tests -v
python tools/build_followup_report.py
```

Phases can run separately: `derive`, `recover`, `compare`, `regularity`,
`density`, `isogeny`, and `consolidate`. Run the timed comparison without
concurrent CPU-intensive work. Timings vary by machine; three repetitions do
not provide a confidence interval. PDF generation requires the system DejaVu
Sans fonts described in the original README. The bundled wheel is the frozen
v1.0.0 baseline; the new preprocessing module is in this source snapshot.

Results and certificates are in [followups/2026-10-15](followups/2026-10-15/).
`derived_cases.json` contains the new, separately hashed fixtures;
`row_equivalence_proofs.json` contains exact source-combination masks;
`complete_mathematical_metrics.json` covers the original 69 systems.
The `backend.solving_degree` field remains null where instrumentation is absent.

## Remaining leads

- Add uniform-target controls in smaller groups with enough pair-sum density
  to observe both successful and unsuccessful decompositions.
- Instrument actual intermediate degree, or add a measured Macaulay/F4 backend.
- Extend relation arity, then measure relation-matrix rank and linear algebra.
- Execute the existing Sage/Singular adapter and compare exact results.
- Compare isogeny- and endomorphism-aware factor bases with matched target and
  subgroup policies.

These experiments stop at relation search; they do not implement a complete
index-calculus discrete-logarithm solver.
