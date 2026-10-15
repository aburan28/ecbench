"""Concrete external-solver export and finite-domain replay protocol.

These helpers do not claim a Sage installation or invent external timings.
The generated adapter actually calls Sage's Singular std implementation.
"""
import json
import math
import time
from pathlib import Path
from types import SimpleNamespace

from .algebra import _basis_roots, _integer, _normalize_case, _record_input, _verify_relations
import hashlib


PROTOCOL = "ecbench-external-basis-v1"


def case_digest(case):
    return hashlib.sha256(json.dumps(case, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def export_sage_case(case, destination, order="grevlex"):
    """Write an executable Sage/Python adapter and its exact case JSON.

    ``sage -python solver.py`` emits one sparse-basis protocol result.  Sage's
    Singular prime-field path has a small-characteristic limit; unsupported
    fields fail explicitly rather than silently changing coefficient fields.
    """
    orders = {"lex": "lex", "grlex": "deglex", "grevlex": "degrevlex"}
    if order not in orders:
        raise ValueError("Unsupported monomial order")
    normalized = _normalize_case(case)
    names, modulus, polynomials, _, _ = _record_input(normalized)
    if modulus > 2_147_483_647:
        raise ValueError("Sage/Singular std adapter supports coefficient characteristics <= 2^31-1; this case needs an arbitrary-characteristic backend")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    spec = {"schema_version": PROTOCOL, "case_sha256": case_digest(case),
            "variables": names, "characteristic": str(modulus), "order": order,
            "sage_order": orders[order],
            "polynomials": [[{"coefficient": str(coefficient), "exponents": list(exponents)}
                             for exponents, coefficient in terms.items()] for terms in polynomials]}
    script = '''# Generated ecbench adapter; execute: sage -python solver.py
import json
import time
from sage.all import GF, PolynomialRing
from sage.version import version
spec = json.loads(SPEC_JSON)
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
'''.replace("SPEC_JSON", repr(json.dumps(spec, sort_keys=True)))
    script_path = destination / "solver.py"
    case_path = destination / "case.json"
    script_path.write_text(script)
    case_path.write_text(json.dumps(case, sort_keys=True, indent=2) + "\n")
    return {"script": str(script_path), "case": str(case_path),
            "command": ["sage", "-python", str(script_path)], "protocol": PROTOCOL,
            "expected_stdout": "one JSON object containing a sparse computed Groebner basis and real backend stage timings"}


def replay_external_basis(case, result):
    """Check external computed roots/relations against this golden case.

    This verifies the finite-domain solution set and subgroup relations, not
    equality of polynomial ideals over an algebraic closure. External basis
    coefficients are never trusted as point relations without group checks.
    Imported backend timings retain their declared source and scope.
    """
    import sympy
    if result.get("schema_version") != PROTOCOL or result.get("case_sha256") != case_digest(case):
        raise ValueError("External result protocol or case hash mismatch")
    normalized = _normalize_case(case)
    names, modulus, _, domains, _ = _record_input(normalized)
    if result["variables"] != names or _integer(result["characteristic"]) != modulus:
        raise ValueError("External polynomial ring differs from the case")
    order = result["backend"]["monomial_order"]
    if order not in {"lex", "grlex", "grevlex"}:
        raise ValueError("External monomial order is unsupported")
    if any(not isinstance(value, (int, float)) or value < 0 or not math.isfinite(value)
           for value in result.get("stage_times", {}).values()):
        raise ValueError("External stage times must be finite nonnegative seconds")
    started = time.perf_counter()
    replay_case = {**normalized, "polynomials": result["polynomials"]}
    _, _, terms, _, _ = _record_input(replay_case)
    symbols = sympy.symbols(names)
    basis = SimpleNamespace(polys=[sympy.Poly.from_dict(polynomial, symbols, modulus=modulus) for polynomial in terms])
    roots = _basis_roots(basis, domains, modulus)
    roots_seconds = time.perf_counter() - started
    started = time.perf_counter()
    relations = _verify_relations(normalized, roots)
    verification_seconds = time.perf_counter() - started
    return {"schema_version": "ecbench-external-replay-v1", "case_id": case["id"],
            "case_sha256": case_digest(case), "backend": result["backend"],
            "status": "ok" if relations["truth_matches"] else "error",
            "relations": relations, "backend_stage_times": result.get("stage_times", {}),
            "replay_stage_times": {"basis_parse_and_root_extraction_seconds": roots_seconds,
                                   "relation_verification_seconds": verification_seconds},
            "verification_scope": "exact roots in declared finite domain and exact selected-subgroup point relations; no polynomial-ideal equality certificate",
            "timing_source": "external adapter self-measured stages; replay stages measured independently here; process wall time must be recorded by the invoking benchmark"}
