"""Known-answer tests for measured solving, relation signs, and timeouts."""
import importlib.util
import unittest

from ecbench.algebra import run_case, run_suite


def toy_case():
    # On y^2=x^3+1/F17, P=(0,1) has order three, so P+P=-P.
    field = {"kind": "prime", "cardinality": "17", "element_bits": 5,
             "coordinate_bytes": 1, "representation": {"modulus": "0x11"}}
    positive, negative = {"x": "0x00", "y": "0x01"}, {"x": "0x00", "y": "0x10"}
    curve = {"coefficients": {"a4": "0x00", "a6": "0x01"},
             "subgroup": {"order": "3", "generator": positive, "r_torsion_rank": 1}}
    return {"id": "toy-sign-relation", "field_id": "F17", "curve_id": "toy",
            "field": field, "curve": curve, "coefficient_modulus": "17", "variables": ["x0", "x1"],
            "polynomials": [{"terms": [{"coefficient": "1", "exponents": [1, 0]}]},
                            {"terms": [{"coefficient": "1", "exponents": [0, 1]}]}],
            "domains": {"x0": ["0", "1"], "x1": ["0", "1"]},
            "coordinate_groups": [{"variables": ["x0"]}, {"variables": ["x1"]}],
            "factor_base": {"points_by_x": {"0x00": [positive, negative]}},
            "target": {"kind": "planted", "point": negative},
            "expected": {"algebraic_assignments": [["0", "0"]], "ordered_count": 1, "swap_orbit_count": 1}}


@unittest.skipUnless(importlib.util.find_spec("sympy"), "optional SymPy backend not installed")
class AlgebraTests(unittest.TestCase):
    def test_groebner_roots_then_group_verification_and_hilbert_degree(self):
        result = run_case(toy_case(), timeout=10)
        self.assertEqual(result["status"], "ok", result)
        self.assertEqual(result["relations"]["root_assignments"], [["0", "0"]])
        self.assertEqual(result["relations"]["ordered_relations"], 1)
        self.assertTrue(result["relations"]["truth_matches"])
        self.assertEqual(result["degree_of_regularity"]["hilbert_function"], [1, 0])
        self.assertEqual(result["degree_of_regularity"]["value"], 1)
        self.assertIsNone(result["backend"]["solving_degree"])
        self.assertGreater(result["stage_times"]["total_solve_seconds"], 0)

    def test_unit_ideal_is_completed_zero_yield_not_timeout(self):
        case = toy_case()
        case["id"] = "toy-unsatisfiable"
        case["polynomials"].append({"terms": [{"coefficient": "1", "exponents": [0, 0]}]})
        case["expected"] = {"algebraic_assignments": [], "ordered_count": 0, "swap_orbit_count": 0}
        result = run_case(case, timeout=10)
        self.assertEqual(result["status"], "ok", result)
        self.assertTrue(result["groebner_basis"]["unit_ideal"])
        self.assertEqual(result["relations"]["polynomial_roots"], 0)
        self.assertEqual(result["degree_of_regularity"]["value"], 0)

    def test_timeout_censored_and_excluded_from_success_denominator(self):
        result = run_suite([toy_case()], timeout=0.000001)
        self.assertEqual(result["results"][0]["status"], "timeout")
        self.assertTrue(result["results"][0]["censored"])
        self.assertIsNone(result["results"][0]["relations"])
        self.assertEqual(result["summary"]["timeouts_censored"], 1)
        self.assertIsNone(result["summary"]["target_success_fraction_completed"])
        self.assertEqual(result["summary"]["all_attempt_success_fraction"], 0)
        self.assertEqual(result["summary"]["all_attempt_success_fraction_upper_bound"], 1)
        self.assertEqual(result["summary"]["relations_per_attempt"], 0)
        self.assertIsNone(result["timing_scope"]["fixture_generation_seconds"])
        self.assertIsNone(result["timing_scope"]["generation_inclusive_total_pipeline_seconds"])

    def test_primary_yield_uses_all_attempts_including_completed_empty_case(self):
        empty = toy_case()
        empty["id"] = "toy-completed-empty"
        empty["polynomials"].append({"terms": [{"coefficient": "1", "exponents": [0, 0]}]})
        empty["expected"] = {"algebraic_assignments": [], "ordered_count": 0, "swap_orbit_count": 0}
        result = run_suite([toy_case(), empty], timeout=10, diagnostics=False)
        self.assertEqual(result["summary"]["completed"], 2)
        self.assertEqual(result["summary"]["all_attempt_success_fraction"], 0.5)
        self.assertEqual(result["summary"]["relations_per_attempt"], 0.5)
        self.assertEqual(result["summary"]["all_attempt_success_fraction_upper_bound"], 0.5)
        self.assertEqual(result["results"][0]["relation_policy"]["factor_base_identity"], "excluded")

    def test_false_truth_is_rejected_after_real_solve(self):
        case = toy_case()
        case["expected"]["ordered_count"] = 99
        result = run_case(case, timeout=10, diagnostics=False)
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["relations"]["ordered_relations"], 1)

    def test_on_curve_point_outside_selected_subgroup_is_rejected(self):
        case = toy_case()
        # (1,6) lies on the curve and has order nine, outside <(0,1)>.
        case["factor_base"]["points_by_x"]["0x01"] = [{"x": "0x01", "y": "0x06"}]
        result = run_case(case, timeout=10, diagnostics=False)
        self.assertEqual(result["status"], "error")
        self.assertIn("outside selected subgroup", result["error"])

    def test_external_replay_enforces_case_hash_and_golden_roots(self):
        from ecbench.backends import PROTOCOL, case_digest, replay_external_basis
        case = toy_case()
        external = {"schema_version": PROTOCOL, "case_sha256": case_digest(case),
                    "variables": case["variables"], "characteristic": "17",
                    "backend": {"name": "known-answer", "monomial_order": "grevlex"},
                    "stage_times": {}, "polynomials": case["polynomials"]}
        self.assertEqual(replay_external_basis(case, external)["status"], "ok")
        external["polynomials"] = [{"terms": [{"coefficient": "1", "exponents": [0, 0]}]}]
        self.assertEqual(replay_external_basis(case, external)["status"], "error")
        external["case_sha256"] = "wrong"
        with self.assertRaisesRegex(ValueError, "hash mismatch"):
            replay_external_basis(case, external)


if __name__ == "__main__":
    unittest.main()
