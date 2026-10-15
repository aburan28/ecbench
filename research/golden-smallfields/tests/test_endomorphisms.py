"""Conformance, cyclic-subgroup selection and corruption tests for GLV2."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

from ecbench.endomorphisms import (apply_endomorphism, benchmark_cases, joint_mul,
                                  scalar_mul, verify_cases)
from ecbench.reference import Curve, Field
from ecbench.validate import load_dataset

ROOT = Path(__file__).resolve().parents[1]


class EndomorphismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_dataset()
        cls.catalog = json.loads((ROOT / 'build' / 'endomorphism_cases.json').read_text())
        cls.cases = {row['curve_id']: row for row in cls.catalog['cases']}
        fields = {record['id']: Field(record) for record in cls.data['fields']}
        cls.curves = {record['id']: Curve(record, fields[record['field_id']]) for record in cls.data['curves']}

    def test_entire_catalog_verifies_480_signed_vectors(self):
        report = verify_cases(self.data, self.catalog)
        self.assertEqual(report['status'], 'pass')
        self.assertEqual(report['curves_classified'], 56)
        self.assertEqual(report['supported_cases'], 30)
        self.assertEqual(report['unsupported_cases'], 26)
        self.assertEqual(report['binary_frobenius_cases'], 6)
        self.assertEqual(report['scalar_vectors_checked'], 480)

    def test_rational_torsion_does_not_establish_selected_cyclic_stability(self):
        for curve_id in ('p64-E3-2', 'p64-E3-3'):
            with self.subTest(curve=curve_id):
                curve, case = self.curves[curve_id], self.cases[curve_id]
                image = apply_endomorphism(curve, case['map'], curve.G)
                self.assertTrue(curve.r_torsion(image))
                self.assertFalse(curve.subgroup_contains(image))
                self.assertEqual(case['characteristic_roots_mod_r'], [])
                self.assertEqual(case['status'], 'unsupported')
                self.assertEqual(case['reason'], 'selected_cyclic_subgroup_not_stable')

    def test_real_signed_components_reconstruct_prime_and_binary_points(self):
        for curve_id in ('p83-E4-2', 'p83-E3-3', 'K83-a0'):
            with self.subTest(curve=curve_id):
                curve, case = self.curves[curve_id], self.cases[curve_id]
                vector = next(row for row in case['scalar_vectors'] if row['name'] == 'fixed_128_bit')
                scalar = int(vector['scalar'])
                components = list(map(int, vector['components']))
                self.assertTrue(any(value < 0 for value in components))
                self.assertTrue(all(value != 0 for value in components))
                image = apply_endomorphism(curve, case['map'], curve.G)
                expected = curve.mul(scalar, curve.G)
                self.assertEqual(joint_mul(curve, components, [curve.G, image]), expected)
                self.assertEqual(curve.add(curve.mul(components[0], curve.G), curve.mul(components[1], image)), expected)
                basis = [list(map(int, row)) for row in case['lattice_basis']]
                for signed_scalar in (-scalar, -1, 0, curve.r, scalar + curve.r):
                    self.assertEqual(scalar_mul(curve, case['map'], basis, signed_scalar), curve.mul(signed_scalar, curve.G))

    def test_false_eigenvalue_and_lattice_index_are_rejected(self):
        for property_name, error in [('eigenvalue', 'Incorrect selected eigenvalue'),
                                     ('lattice_determinant', 'Incorrect lattice determinant')]:
            with self.subTest(property=property_name):
                damaged = deepcopy(self.catalog)
                case = next(row for row in damaged['cases'] if row['status'] == 'supported')
                case[property_name] = str(int(case[property_name]) + 1)
                with self.assertRaisesRegex(ValueError, error):
                    verify_cases(self.data, damaged)

    def test_sign_corruption_is_rejected_before_point_reconstruction(self):
        damaged = deepcopy(self.catalog)
        case = next(row for row in damaged['cases'] if row['curve_id'] == 'p83-E4-2')
        vector = next(row for row in case['scalar_vectors'] if row['name'] == 'fixed_128_bit')
        self.assertLess(int(vector['components'][0]), 0)
        vector['components'][0] = str(-int(vector['components'][0]))
        with self.assertRaisesRegex(ValueError, 'Incorrect scalar decomposition'):
            verify_cases(self.data, damaged)

    def test_source_parameters_and_unstable_classification_cannot_be_silently_changed(self):
        damaged_data = deepcopy(self.data)
        damaged_data['curves'][0]['subgroup']['order'] = '3'
        with self.assertRaisesRegex(ValueError, 'source hash mismatch'):
            verify_cases(damaged_data, self.catalog)
        damaged = deepcopy(self.catalog)
        case = next(row for row in damaged['cases'] if row['curve_id'] == 'p64-E3-2')
        case['status'] = 'supported'
        with self.assertRaisesRegex(ValueError, 'Incorrect unstable-subgroup classification'):
            verify_cases(self.data, damaged)

    def test_timing_scopes_distinguish_precomputed_and_full_paths(self):
        result = benchmark_cases(self.data, self.catalog, ids={'p61-E3-2'}, samples=1, iterations=1)
        operations = {row['operation'] for row in result['results']}
        self.assertEqual(operations, {'endomorphism', 'scalar_decompose', 'scalar_mul_subgroup_baseline',
                                      'glv2_joint_mul_precomputed', 'glv2_scalar_mul_full'})
        self.assertEqual(result['preflight']['status'], 'pass')
        self.assertIn('excludes decomposition and endomorphism', result['precomputed_scope'])
        self.assertTrue(all(row['median'] > 0 for row in result['results']))


if __name__ == '__main__':
    unittest.main()
