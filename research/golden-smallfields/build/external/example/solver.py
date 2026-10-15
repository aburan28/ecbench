# Generated ecbench adapter; execute: sage -python solver.py
import json
import time
from sage.all import GF, PolynomialRing
from sage.version import version
spec = json.loads('{"case_sha256": "20ea7a98672ea91bb26a51b8ffe5160431def3a28b1d2cc85f1a383ef93d3125", "characteristic": "2", "order": "grevlex", "polynomials": [[{"coefficient": "1", "exponents": [0, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}, {"coefficient": "1", "exponents": [0, 1, 0, 1]}], [{"coefficient": "1", "exponents": [0, 1, 1, 0]}, {"coefficient": "1", "exponents": [1, 0, 0, 1]}], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}, {"coefficient": "1", "exponents": [0, 1, 1, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}, {"coefficient": "1", "exponents": [1, 0, 0, 1]}], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}, {"coefficient": "1", "exponents": [1, 0, 1, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}], [{"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 1, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}, {"coefficient": "1", "exponents": [1, 0, 0, 1]}], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}], [], [{"coefficient": "1", "exponents": [1, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}, {"coefficient": "1", "exponents": [1, 0, 1, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}], [{"coefficient": "1", "exponents": [0, 1, 1, 0]}, {"coefficient": "1", "exponents": [1, 0, 0, 1]}], [{"coefficient": "1", "exponents": [2, 0, 0, 0]}, {"coefficient": "1", "exponents": [1, 0, 0, 0]}], [{"coefficient": "1", "exponents": [0, 2, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}], [{"coefficient": "1", "exponents": [0, 0, 2, 0]}, {"coefficient": "1", "exponents": [0, 0, 1, 0]}], [{"coefficient": "1", "exponents": [0, 0, 0, 2]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}], [{"coefficient": "1", "exponents": [0, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 1, 0, 0]}, {"coefficient": "1", "exponents": [1, 1, 0, 0]}], [{"coefficient": "1", "exponents": [0, 0, 0, 0]}, {"coefficient": "1", "exponents": [0, 0, 0, 1]}, {"coefficient": "1", "exponents": [0, 0, 1, 1]}]], "sage_order": "degrevlex", "schema_version": "ecbench-external-basis-v1", "variables": ["x1_0", "x1_1", "x2_0", "x2_1"]}')
started = time.perf_counter()
field = GF(int(spec["characteristic"]))
ring = PolynomialRing(field, names=tuple(spec["variables"]), order=spec["sage_order"])
polynomials = [ring({tuple(term["exponents"]): field(int(term["coefficient"])) for term in terms}) for terms in spec["polynomials"]]
setup_seconds = time.perf_counter() - started
started = time.perf_counter()
basis = ring.ideal(polynomials).groebner_basis(algorithm="singular:std")
groebner_seconds = time.perf_counter() - started
result = {
    "schema_version": spec["schema_version"], "case_sha256": spec["case_sha256"],
    "variables": spec["variables"], "characteristic": spec["characteristic"],
    "backend": {"name": "sage-singular", "sage_version": version, "method": "singular:std", "monomial_order": spec["order"], "solving_degree": None},
    "stage_times": {"polynomial_setup_seconds": setup_seconds, "groebner_seconds": groebner_seconds},
    "polynomials": [{"terms": [{"coefficient": str(int(coefficient)), "exponents": [int(e) for e in exponents]} for exponents, coefficient in polynomial.dict().items()]} for polynomial in basis],
}
print(json.dumps(result, sort_keys=True))
