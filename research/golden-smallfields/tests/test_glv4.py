"""Proof and conformance tests for the genuine GLV+GLS supplement."""
import copy
import unittest

from ecbench.extension import Fp2
from ecbench.glv4 import _load,verify_cases,benchmark_cases


class GLV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.case=_load()[0]

    def test_quadratic_field_reference(self):
        field=Fp2(7,3)
        for value in range(1,field.q):
            self.assertEqual(field.mul(value,field.inverse(value)),1)
            self.assertEqual(field.power(value,7),field.frobenius(value))
            root=field.sqrt(field.square(value))
            self.assertEqual(field.square(root),field.square(value))
        with self.assertRaises(ValueError):field.inverse(0)
        with self.assertRaises(ValueError):field.unpack(field.q)

    def test_full_case_and_uniform_component_bounds(self):
        summary=verify_cases(self.case)['cases'][0]
        self.assertEqual(summary['status'],'pass')
        self.assertEqual(summary['subgroup_bits'],80)
        self.assertEqual(summary['guaranteed_glv4_component_bits'],20)
        self.assertEqual(summary['guaranteed_glv2_component_bits'],40)
        self.assertEqual(summary['point_count_proof']['incorrect_traces_eliminated'],5)

    def test_wrong_gls_eigenvalue_rejected(self):
        case=copy.deepcopy(self.case)
        row=case['endomorphisms']['psi']
        row['subgroup_eigenvalue']=str(int(row['subgroup_eigenvalue'])+1)
        with self.assertRaises(ValueError):verify_cases(case)

    def test_lattice_and_decomposition_tamper_rejected(self):
        case=copy.deepcopy(self.case)
        row=case['lattices']['glv4']['basis'][0]
        row[0]=str(int(row[0])+1)
        with self.assertRaises(ValueError):verify_cases(case)
        case=copy.deepcopy(self.case)
        row=case['decomposition_vectors'][-1]['glv4_components']
        row[0]=str(int(row[0])+1)
        with self.assertRaises(ValueError):verify_cases(case)

    def test_primality_witness_tamper_rejected(self):
        case=copy.deepcopy(self.case)
        case['certificates'][case['curve']['subgroup_order']]['witness']='1'
        with self.assertRaises(ValueError):verify_cases(case)

    def test_benchmark_measures_complete_equal_paths(self):
        result=benchmark_cases(self.case,samples=1,iterations=1)
        sample=result['cases'][0]['samples'][0]
        self.assertTrue(sample['outputs_equal'])
        for name in ('ordinary_total_ns','glv2_decomposition_ns','glv2_maps_ns',
                     'glv2_joint_multiplication_ns','glv4_decomposition_ns',
                     'glv4_maps_ns','glv4_joint_multiplication_ns'):
            self.assertGreater(sample['timings_ns_per_iteration'][name],0)
        with self.assertRaises(ValueError):benchmark_cases(self.case,samples=0)


if __name__=='__main__':unittest.main()
