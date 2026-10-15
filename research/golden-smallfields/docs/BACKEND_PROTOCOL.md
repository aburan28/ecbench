# External solver export and replay

The primary runnable baseline is SymPy. `ecbench.backends` additionally offers
an actual Sage/Singular adapter, with an explicit sparse-basis replay contract.
Sage is optional and is not installed or timed by creating an export.

```python
from ecbench.backends import export_sage_case, replay_external_basis
import json

info = export_sage_case(case, "build/external/example", order="grevlex")
print(info["command"])
# Execute that argument vector: sage -python build/external/example/solver.py
# Capture its one JSON object, for example in build/external/example/basis.json.
result = json.load(open("build/external/example/basis.json"))
verified = replay_external_basis(case, result)
```

The adapter runs `ideal.groebner_basis(algorithm="singular:std")`, times
polynomial setup and basis computation with `time.perf_counter`, and emits
`ecbench-external-basis-v1` JSON. Sage/Singular's small prime-characteristic
path supports these GF(2) cases; this adapter explicitly rejects characteristics
larger than `2^31-1`. The 61–83-bit prime-field cases need an external backend
with arbitrary prime-characteristic support. They are handled by the primary
SymPy baseline. An unsupported export never substitutes GF(2), QQ, or a smaller
prime for the requested coefficient field.

The result contract is:

```json
{
  "schema_version": "ecbench-external-basis-v1",
  "case_sha256": "SHA-256 of compact sorted-key original case JSON",
  "variables": ["x1_0", "x2_0"],
  "characteristic": "2",
  "backend": {
    "name": "sage-singular",
    "sage_version": "the installed Sage version",
    "method": "singular:std",
    "monomial_order": "grevlex",
    "solving_degree": null
  },
  "stage_times": {
    "polynomial_setup_seconds": 0.001,
    "groebner_seconds": 0.001
  },
  "polynomials": [
    {"terms": [{"coefficient": "1", "exponents": [1, 0]}]}
  ]
}
```

The numeric times above illustrate the JSON format only. They are not benchmark
results. Coefficients are decimal integers; exponent vectors follow the
declared variable order. Other backends can emit this same contract with their
actual name, version, method, order, computed basis, and measured stages.

Replay checks the case hash and coefficient field, extracts roots from the
**supplied computed basis** over the declared finite domains, and independently
verifies signed point decompositions and selected-subgroup membership. It
compares the exact finite-domain roots and point relations to golden truth.
It does not certify equality of ideals over an algebraic closure; a backend
requiring that stronger guarantee should supply transformation or ideal
membership certificates. Imported timings remain labeled as external
adapter measurements, while replay stages have independent local timings.
The caller must impose a subprocess timeout and record observed process wall
time for an external performance comparison; exporting or replaying alone
does not produce that end-to-end measurement.
