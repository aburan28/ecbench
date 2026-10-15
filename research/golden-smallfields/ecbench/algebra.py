"""Isolated, measured polynomial relation-search benchmarks.

SymPy is an optional, deliberately simple baseline.  Ground-truth enumerations
in a case are used for validation only, never as solver output.
"""
from __future__ import annotations

import collections
import datetime
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time


RESULT_VERSION = "ecbench-algebra-results-v1"
RELATION_POLICY = {
    "decomposition_arity": 2,
    "repeated_factors": "allowed; a repeated signed point has swap-orbit size one",
    "factor_base_identity": "excluded",
    "target_identity": "reference arithmetic permits identity; current golden fixtures use nonidentity targets",
    "signed_sum": "P + Q must equal the supplied target exactly",
    "relation_counting": "distinct unordered signed point decompositions, with ordered decompositions retained separately",
}


def _timing_scope():
    return {
        "measurement_mode": "prebuilt-fixture solving; fixture generation is not measured",
        "total_solve_seconds": "worker normalization and parsing of prebuilt polynomial equations + polynomial setup + Groebner basis + evaluation of computed basis over declared finite domain + selected-subgroup checks and signed point relation verification; fixture generation, Semaev construction/Weil descent, factor-base and target construction, golden-truth enumeration, SymPy import, process startup, input JSON serialization, and degree-regularity diagnostic excluded",
        "observed_subprocess_wall_seconds": "process launch through completion or forced kill; includes import/startup and diagnostic, excludes parent input serialization and fixture generation",
        "suite_wall_seconds": "iteration over prebuilt cases including input serialization and subprocess orchestration; excludes fixture generation and loading the suite JSON",
        "fixture_generation_seconds": None,
        "fixture_generation_status": "unmeasured; no generation-inclusive total pipeline time is claimed",
        "generation_inclusive_total_pipeline_seconds": None,
        "root_extraction": "evaluate computed Groebner basis over finite allowed domain; independent truth enumeration is never consumed as solver roots",
        "degree_of_regularity": "first zero graded component of homogeneous highest-degree input-part ideal; final basis maximum degree is reported separately; solving degree is unavailable",
    }


def _integer(value):
    return int(value, 0) if isinstance(value, str) and value.startswith("0x") else int(value)


def _polynomial_terms(polynomial):
    """Canonical sparse terms; also accepts a direct term list."""
    if isinstance(polynomial, dict):
        polynomial = polynomial["terms"]
    for term in polynomial:
        if isinstance(term, dict):
            yield _integer(term.get("coefficient", term.get("coeff"))), tuple(term["exponents"])
        else:
            yield _integer(term[0]), tuple(term[1])


def _normalize_case(case):
    """Adapt the standalone raw case-export contract without changing data."""
    if "polynomial_ring" in case:
        return case
    normalized = dict(case)
    names = case["variables"]
    normalized["polynomial_ring"] = {"variables": names, "characteristic": case["coefficient_modulus"]}
    normalized["variable_domains"] = case["domains"]
    groups = []
    for original in case["coordinate_groups"]:
        group = dict(original)
        group["variables"] = [names.index(value) if isinstance(value, str) else value for value in group["variables"]]
        group.setdefault("encoding", "binary_linear" if case["field"]["kind"] == "binary" else "prime_x")
        groups.append(group)
    normalized["coordinate_groups"] = groups
    factor_base = dict(case["factor_base"])
    if isinstance(factor_base["points_by_x"], dict):
        factor_base["points_by_x"] = [{"x": x, "points": points} for x, points in factor_base["points_by_x"].items()]
    normalized["factor_base"] = factor_base
    target = dict(case["target"])
    target.setdefault("sampling", target.get("kind", "unspecified"))
    normalized["target"] = target
    truth = case.get("expected", {})
    normalized["ground_truth"] = {
        **({"root_assignments": [[_integer(value) for value in assignment] for assignment in truth["algebraic_assignments"]],
            "polynomial_roots": len(truth["algebraic_assignments"])} if "algebraic_assignments" in truth else {}),
        **({"ordered_relations": _integer(truth["ordered_count"])} if "ordered_count" in truth else {}),
        **({"unordered_relations": _integer(truth["swap_orbit_count"])} if "swap_orbit_count" in truth else {}),
        **({"ordered_point_decompositions": truth["ordered_point_decompositions"]} if "ordered_point_decompositions" in truth else {}),
    }
    return normalized


def _record_input(case):
    modulus = _integer(case["polynomial_ring"]["characteristic"])
    names = case["polynomial_ring"]["variables"]
    if modulus < 2 or not names or len(set(names)) != len(names):
        raise ValueError("Invalid polynomial ring")
    polynomials = []
    maximum_degree = 0
    role_metrics = {}
    for polynomial in case["polynomials"]:
        normalized = collections.defaultdict(int)
        for coefficient, exponents in _polynomial_terms(polynomial):
            if len(exponents) != len(names) or any(type(e) is not int or e < 0 for e in exponents):
                raise ValueError("Invalid monomial exponent vector")
            normalized[exponents] = (normalized[exponents] + coefficient) % modulus
        terms = {e: c for e, c in normalized.items() if c}
        if terms:
            maximum_degree = max(maximum_degree, max(map(sum, terms)))
        polynomials.append(terms)
        role = polynomial.get("role", "unspecified") if isinstance(polynomial, dict) else "unspecified"
        metrics = role_metrics.setdefault(role, {"polynomials": 0, "terms": 0, "max_total_degree": 0})
        metrics["polynomials"] += 1
        metrics["terms"] += len(terms)
        metrics["max_total_degree"] = max(metrics["max_total_degree"], max(map(sum, terms), default=0))
    domains = case["variable_domains"]
    if isinstance(domains, dict):
        domains = [domains[name] for name in names]
    if len(domains) != len(names) or any(not values for values in domains):
        raise ValueError("Each variable needs a nonempty finite domain")
    domains = [tuple(_integer(value) for value in values) for values in domains]
    if any(len(set(values)) != len(values) or any(value < 0 or value >= modulus for value in values) for values in domains):
        raise ValueError("Noncanonical or duplicated variable-domain elements")
    return names, modulus, polynomials, domains, {
        "variables": len(names), "polynomials": len(polynomials),
        "max_total_degree": maximum_degree,
        "terms": sum(map(len, polynomials)),
        "domain_grid_size": math.prod(map(len, domains)),
        "polynomial_roles": role_metrics,
    }


def _eval_sparse(terms, values, modulus):
    answer = 0
    for exponents, coefficient in terms:
        value = coefficient
        for x, exponent in zip(values, exponents):
            if exponent:
                value = value * pow(x, exponent, modulus) % modulus
        answer = (answer + value) % modulus
    return answer


def _basis_roots(basis, domains, modulus):
    terms = [[(exponents, int(coefficient) % modulus) for exponents, coefficient in polynomial.terms()] for polynomial in basis.polys]
    return [list(values) for values in itertools.product(*domains)
            if all(_eval_sparse(polynomial, values, modulus) == 0 for polynomial in terms)]


def _degree_regularity(sympy, polynomials, symbols, modulus, method, order):
    """Exact Hilbert function of the homogeneous top-part ideal.

    Enumeration uses pure-power bounds in its leading ideal, rather than
    treating an observed final basis degree as a solving degree.
    """
    top = []
    for terms in polynomials:
        if not terms:
            continue
        maximum = max(map(sum, terms))
        selected = {exponents: coefficient for exponents, coefficient in terms.items() if sum(exponents) == maximum}
        top.append(sympy.Poly.from_dict(selected, symbols, modulus=modulus).as_expr())
    basis = sympy.groebner(top, *symbols, modulus=modulus, method=method, order=order)
    leading = [tuple(p.LM(order=basis.order).exponents) for p in basis.polys]
    if any(not any(exponents) for exponents in leading):
        return {"value": 0, "hilbert_function": [0], "homogeneous_basis_size": len(basis.polys),
                "definition": "first zero graded component of quotient by homogeneous highest-degree input parts"}
    bounds = []
    for i in range(len(symbols)):
        pure = [e[i] for e in leading if e[i] and sum(e) == e[i]]
        if not pure:
            return {"value": None, "hilbert_function": None,
                    "reason": "homogeneous_top_part_ideal_not_zero_dimensional"}
        bounds.append(min(pure))
    if math.prod(bounds) > 1_000_000:
        return {"value": None, "hilbert_function": None, "reason": "standard_monomial_enumeration_limit"}
    histogram = collections.Counter()
    for exponents in itertools.product(*(range(bound) for bound in bounds)):
        if not any(all(a >= b for a, b in zip(exponents, monomial)) for monomial in leading):
            histogram[sum(exponents)] += 1
    hilbert = [histogram[d] for d in range(max(histogram, default=-1) + 2)]
    return {"value": next(d for d, dimension in enumerate(hilbert) if dimension == 0),
            "hilbert_function": hilbert, "quotient_dimension": sum(hilbert),
            "homogeneous_basis_size": len(basis.polys),
            "definition": "first zero graded component of quotient by homogeneous highest-degree input parts"}


def _decode_x(case, values):
    """Decode factor-base coordinates from the solver's variable assignments."""
    groups = case["coordinate_groups"]
    decoded = []
    for group in groups:
        if group["encoding"] == "prime_x":
            decoded.append(values[group["variables"][0]])
        elif group["encoding"] == "binary_linear":
            result = 0
            for index, basis_element in zip(group["variables"], group["basis"]):
                if values[index]:
                    result ^= _integer(basis_element)
            decoded.append(result)
        else:
            raise ValueError("Unsupported coordinate encoding")
    if len(decoded) != 2:
        raise ValueError("This baseline requires two factor-base points")
    return decoded


def _verify_relations(case, roots):
    """Check decompositions with independent group arithmetic after GB solving."""
    from .reference import Curve, Field
    field = Field(case["field"])
    curve = Curve(case["curve"], field)
    target = curve.parse_point(case["target"]["point"])
    points_by_x = {}
    membership = {}

    def subgroup_contains(point):
        if point not in membership:
            membership[point] = curve.subgroup_contains(point)
            # Membership in a subgroup is invariant under taking negatives.
            membership[curve.neg(point)] = membership[point]
        return membership[point]

    for row in case["factor_base"]["points_by_x"]:
        x = _integer(row["x"])
        points = [curve.parse_point(point) for point in row["points"]]
        for point in points:
            if point is None or point[0] != x or not curve.on_curve(point):
                raise ValueError("Invalid factor-base point")
            if not subgroup_contains(point):
                raise ValueError("Factor-base point is outside selected subgroup")
        points_by_x[x] = points
    if not curve.on_curve(target):
        raise ValueError("Target is not on the curve")
    if not subgroup_contains(target):
        raise ValueError("Target is outside selected subgroup")
    ordered = set()
    unordered = set()
    successful_assignments = []
    for values in roots:
        x0, x1 = _decode_x(case, values)
        found = False
        for p, q in itertools.product(points_by_x.get(x0, []), points_by_x.get(x1, [])):
            if curve.add(p, q) == target:
                ordered.add((p, q))
                unordered.add(tuple(sorted((p, q))))
                found = True
        if found:
            successful_assignments.append(values)
    expected = case.get("ground_truth", {})
    checks = {}
    comparisons = {
        "polynomial_roots": len(roots),
        "ordered_relations": len(ordered),
        "unordered_relations": len(unordered),
        "relation_root_assignments": len(successful_assignments),
    }
    for name, value in comparisons.items():
        if name in expected:
            checks[name] = value == expected[name]
    if "root_assignments" in expected:
        checks["root_assignments"] = sorted(roots) == sorted(expected["root_assignments"])
    if "ordered_point_decompositions" in expected:
        reference = {(curve.parse_point(row["P"]), curve.parse_point(row["Q"]))
                     for row in expected["ordered_point_decompositions"]}
        checks["ordered_point_decompositions"] = ordered == reference
    return {
        "polynomial_roots": len(roots), "root_assignments": [[str(value) for value in assignment] for assignment in roots],
        "relation_root_assignments": len(successful_assignments),
        "successful_root_assignments": [[str(value) for value in assignment] for assignment in successful_assignments],
        "ordered_relations": len(ordered), "unordered_relations": len(unordered),
        "target_success": bool(ordered),
        "truth_checks": checks,
        "truth_matches": all(checks.values()) if checks else None,
        "selected_subgroup_membership_verified": True,
        "subgroup_membership_cache_entries": len(membership),
        "ordered_relation_points": [[curve.point(p), curve.point(q)] for p, q in sorted(ordered)],
    }


def _worker(case, method, order, diagnostics):
    import sympy
    solve_start = time.perf_counter()
    stages = {}

    def emit(event, **payload):
        print(json.dumps({"event": event, **payload}, sort_keys=True), flush=True)

    started = time.perf_counter()
    case = _normalize_case(case)
    names, modulus, polynomials, domains, input_metrics = _record_input(case)
    symbols = sympy.symbols(names)
    expressions = [sympy.Poly.from_dict(terms, symbols, modulus=modulus).as_expr() for terms in polynomials]
    stages["parse_and_polynomial_setup_seconds"] = time.perf_counter() - started
    backend = {"name": "sympy", "version": sympy.__version__, "method": method, "monomial_order": order,
               "characteristic": str(modulus), "solving_degree": None,
               "solving_degree_reason": "SymPy does not expose maximal degree reached during computation"}
    emit("prepared", input=input_metrics, backend=backend, stage_times=stages)
    started = time.perf_counter()
    basis = sympy.groebner(expressions, *symbols, modulus=modulus, method=method, order=order)
    stages["groebner_seconds"] = time.perf_counter() - started
    basis_metrics = {"size": len(basis.polys),
                     "max_total_degree": max((p.total_degree() for p in basis.polys), default=0),
                     "terms": sum(len(p.terms()) for p in basis.polys),
                     "zero_dimensional": bool(basis.is_zero_dimensional),
                     "unit_ideal": any(p.total_degree() == 0 and not p.is_zero for p in basis.polys)}
    # A unit ideal has no roots; SymPy's is_zero_dimensional is false for it.
    emit("groebner_complete", groebner_basis=basis_metrics, stage_times=stages)
    started = time.perf_counter()
    roots = [] if basis_metrics["unit_ideal"] else _basis_roots(basis, domains, modulus)
    stages["root_extraction_seconds"] = time.perf_counter() - started
    emit("roots_complete", polynomial_roots=len(roots), stage_times=stages)
    started = time.perf_counter()
    relation_metrics = _verify_relations(case, roots)
    stages["relation_verification_seconds"] = time.perf_counter() - started
    stages["total_solve_seconds"] = time.perf_counter() - solve_start
    result = {"status": "ok", "input": input_metrics, "backend": backend, "groebner_basis": basis_metrics,
              "relations": relation_metrics, "stage_times": dict(stages),
              "degree_of_regularity": {"value": None, "reason": "diagnostics_disabled"}}
    emit("solve_complete", result=result)
    if diagnostics:
        started = time.perf_counter()
        result["degree_of_regularity"] = _degree_regularity(sympy, polynomials, symbols, modulus, method, order)
        result["stage_times"]["degree_regularity_diagnostic_seconds"] = time.perf_counter() - started
    emit("complete", result=result)


def run_case(case, timeout=5.0, method="buchberger", order="grevlex", diagnostics=True):
    """Run one case with a hard subprocess wall-time budget (seconds)."""
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be a positive finite number")
    if method not in {"buchberger", "f5b"} or order not in {"lex", "grlex", "grevlex"}:
        raise ValueError("Unsupported SymPy method or monomial order")
    preparation = time.perf_counter()
    payload = json.dumps({"case": case, "method": method, "order": order, "diagnostics": bool(diagnostics)}, separators=(",", ":"))
    environment = dict(os.environ)
    package_root = str(Path(__file__).resolve().parent.parent)
    environment["PYTHONPATH"] = package_root + os.pathsep + environment.get("PYTHONPATH", "")
    preparation_seconds = time.perf_counter() - preparation
    started = time.perf_counter()
    process = subprocess.Popen([sys.executable, "-m", "ecbench.algebra", "--worker"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, env=environment)
    timed_out = False
    try:
        stdout, stderr = process.communicate(payload, timeout=max(1e-9, timeout - (time.perf_counter() - started)))
    except subprocess.TimeoutExpired:
        timed_out = True
        process.kill()
        stdout, stderr = process.communicate()
    finally:
        for stream in (process.stdin, process.stdout, process.stderr):
            if stream is not None and not stream.closed:
                try:
                    stream.close()
                except OSError:
                    pass
    observed_wall = time.perf_counter() - started
    events = []
    for line in stdout.splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    finished = next((event["result"] for event in reversed(events) if event.get("event") == "complete"), None)
    solved = next((event["result"] for event in reversed(events) if event.get("event") == "solve_complete"), None)
    if finished:
        result = finished
    elif solved:
        # Preserve a completed solve if only the optional diagnostic timed out.
        result = solved
        result["degree_of_regularity"] = {"value": None, "reason": "diagnostic_timeout" if timed_out else "diagnostic_error"}
        result["diagnostic_status"] = "timeout" if timed_out else "error"
    else:
        result = {"status": "timeout" if timed_out else "error", "stage_times": {}, "relations": None,
                  "groebner_basis": None, "degree_of_regularity": {"value": None, "reason": "solve_incomplete"}}
        for event in events:
            for name in ("input", "backend", "groebner_basis", "stage_times"):
                if name in event:
                    result[name] = event[name]
        if not timed_out:
            result["error"] = stderr.strip()[-4000:] or f"Worker exited {process.returncode} before completion"
    if "backend" not in result:
        try:
            version = importlib.metadata.version("sympy")
        except importlib.metadata.PackageNotFoundError:
            version = None
        result["backend"] = {"name": "sympy", "version": version, "method": method, "monomial_order": order,
                             "solving_degree": None,
                             "solving_degree_reason": "SymPy does not expose maximal degree reached during computation"}
    if result.get("relations") and result["relations"]["truth_matches"] is False:
        result["status"] = "error"
        result["error"] = "Computed roots or relations disagree with case ground truth"
    result.update({"case_id": case["id"], "curve_id": case.get("curve_id"),
                   "field_id": case.get("field_id") or case.get("field", {}).get("id") or case.get("curve", {}).get("field_id"),
                   "cohort": case.get("target", {}).get("sampling", case.get("target", {}).get("kind", case.get("cohort", "unspecified"))),
                   "timeout_seconds": timeout, "observed_subprocess_wall_seconds": observed_wall,
                   "censored": result["status"] == "timeout", "worker_exit_code": process.returncode,
                   "input_serialization_seconds": preparation_seconds,
                   "case_sha256": hashlib.sha256(json.dumps(case, sort_keys=True, separators=(",", ":")).encode()).hexdigest()})
    for name in ("model_variant", "basis_kind", "factor_base_dimension", "subspace_dimension"):
        if name in case:
            result[name] = case[name]
    if "factor_base" in case:
        by_x = case["factor_base"]["points_by_x"]
        result["factor_base_counts"] = {
            "x_coordinates": len(by_x),
            "points": sum(len(points) for points in by_x.values()) if isinstance(by_x, dict)
            else sum(len(row["points"]) for row in by_x),
        }
    result["cohort_dimensions"] = {
        "target_sampling": result["cohort"],
        "basis_kind": case.get("basis_kind", case.get("factor_base", {}).get("basis_kind", "unspecified")),
        "model_variant": case.get("model_variant", "unspecified"),
    }
    result["relation_policy"] = dict(RELATION_POLICY)
    return result


def _aggregate(rows):
    completed = [row for row in rows if row["status"] == "ok"]
    relation_rows = [row["relations"] for row in completed]
    elapsed = [row["stage_times"]["total_solve_seconds"] for row in completed]
    solve_sum = sum(elapsed)
    successful = sum(row["target_success"] for row in relation_rows)
    ordered_relations = sum(row["ordered_relations"] for row in relation_rows)
    unordered_relations = sum(row["unordered_relations"] for row in relation_rows)
    unfinished = len(rows) - len(completed)
    return {
        "cases": len(rows), "completed": len(completed),
        "completion_fraction": len(completed) / len(rows) if rows else None,
        "timeouts_censored": sum(row["status"] == "timeout" for row in rows),
        "errors": sum(row["status"] == "error" for row in rows),
        "diagnostic_timeouts": sum(row.get("diagnostic_status") == "timeout" for row in rows),
        "successful_targets_completed": successful,
        "all_attempt_success_fraction": successful / len(rows) if rows else None,
        "all_attempt_success_fraction_upper_bound": (successful + unfinished) / len(rows) if rows else None,
        "relations_per_attempt": unordered_relations / len(rows) if rows else None,
        "ordered_relations_per_attempt": ordered_relations / len(rows) if rows else None,
        "primary_metrics_definition": "verified successful targets and unordered signed decompositions / all attempted targets; lower bounds when unfinished cases exist; unfinished cases remain censored",
        "target_success_fraction_completed": successful / len(completed) if completed else None,
        "lower_bound_yield_per_attempted_target": successful / len(rows) if rows else None,
        "lower_bound_yield_definition": "verified successful completed targets / all attempted targets; unfinished targets remain censored",
        "ordered_relations_completed": ordered_relations,
        "unordered_relations_completed": unordered_relations,
        "unordered_relations_per_solve_second_completed": unordered_relations / solve_sum if solve_sum else None,
        "total_solve_seconds_completed": solve_sum,
        "median_solve_seconds_completed": statistics.median(elapsed) if elapsed else None,
        "observed_subprocess_wall_seconds_all": sum(row["observed_subprocess_wall_seconds"] for row in rows),
        "yield_denominator": "primary success/yield metrics use all attempted targets; completed-only rates are explicitly suffixed completed; unfinished solves stay censored and do not prove zero mathematical yield",
    }


def run_suite(cases, timeout=5.0, method="buchberger", order="grevlex", max_cases=None, diagnostics=True,
              progress=None):
    """Run deterministic cases and report planted/random cohorts separately.

    ``progress`` receives each completed case result; it is not required.
    ``cases`` may be a list or a case-suite object containing ``cases``.
    """
    suite = cases if isinstance(cases, dict) else {}
    selected = suite.get("cases", cases)
    selected = list(selected)
    if max_cases is not None:
        if max_cases < 1:
            raise ValueError("max_cases must be positive")
        selected = selected[:max_cases]
    rows = []
    started = time.perf_counter()
    for case in selected:
        row = run_case(case, timeout, method, order, diagnostics)
        rows.append(row)
        if progress:
            progress(row)
    cohorts = collections.defaultdict(list)
    strata = collections.defaultdict(list)
    for row in rows:
        cohorts[row["cohort"]].append(row)
        strata[tuple(sorted(row["cohort_dimensions"].items()))].append(row)
    return {"schema_version": RESULT_VERSION, "suite_id": suite.get("suite_id", suite.get("id", suite.get("schema"))),
            "suite_seed": suite.get("seed"),
            "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "environment": {"python": platform.python_version(), "platform": platform.platform(), "machine": platform.machine()},
            "configuration": {"backend": "sympy", "method": method, "order": order,
                              "per_case_wall_timeout_seconds": timeout, "degree_regularity_diagnostics": bool(diagnostics)},
            "timing_scope": _timing_scope(),
            "primary_metrics": ["all_attempt_success_fraction", "relations_per_attempt"],
            "suite_wall_seconds": time.perf_counter() - started,
            "summary": _aggregate(rows), "cohorts": {key: _aggregate(values) for key, values in sorted(cohorts.items())},
            "stratified_cohorts": [{"cohort": dict(key), **_aggregate(values)} for key, values in sorted(strata.items())],
            "results": rows}


def _main():
    if len(sys.argv) != 2 or sys.argv[1] != "--worker":
        raise SystemExit("Use ecbench algebra-bench; this module's --worker entry is internal")
    payload = json.load(sys.stdin)
    try:
        _worker(payload["case"], payload["method"], payload["order"], payload["diagnostics"])
    except Exception as error:
        print(json.dumps({"event": "error", "error": f"{type(error).__name__}: {error}"}), flush=True)
        raise


if __name__ == "__main__":
    _main()
