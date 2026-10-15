import copy
import unittest

from ecbench.algebra import _normalize_case, _record_input
from ecbench.cli import catalog
from ecbench.preprocess import derive_case, digest, verify_transformation


class PreprocessingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = {case['id']: case for case in catalog('algebra_cases')['cases']}
        cls.case = cls.cases['K83-a0-f3-low_monomials-d6-planted']

    def test_exposes_hidden_linear_equations_with_two_way_proof(self):
        changed, proof = derive_case(self.case)
        self.assertEqual(verify_transformation(self.case, changed, proof)['status'], 'pass')
        originals = _record_input(_normalize_case(self.case))[2]
        original_linear = sum(max(map(sum, poly), default=0) == 1 for poly in originals)
        self.assertGreater(proof['new_linear_generators'], original_linear)
        self.assertEqual(proof['new_linear_generators'], 4)

    def test_forged_unit_rejected_even_with_updated_case_hash(self):
        changed, proof = derive_case(self.case)
        changed['polynomials'][0]['terms'] = [{'coefficient': '1', 'exponents': [0] * len(changed['variables'])}]
        proof['derived_case_sha256'] = digest(changed)
        proof['unit_generator_found'] = True
        proof['new_linear_generators'] -= 1 if max(sum(term['exponents']) for term in derive_case(self.case)[0]['polynomials'][0]['terms']) == 1 else 0
        with self.assertRaises(ValueError):
            verify_transformation(self.case, changed, proof)

    def test_source_combination_witness_must_be_correct(self):
        changed, proof = derive_case(self.case)
        proof['output_xor_input_masks'][0] = '0x0'
        with self.assertRaisesRegex(ValueError, 'source combination'):
            verify_transformation(self.case, changed, proof)

    def test_missing_generator_rejected_by_reverse_inclusion(self):
        changed, proof = derive_case(self.case)
        # Preserve every listed witness, but throw away a necessary row.
        changed['polynomials'] = changed['polynomials'][:-1]
        proof['output_xor_input_masks'] = proof['output_xor_input_masks'][:-1]
        proof['output_generator_count'] -= 1
        proof['derived_case_sha256'] = digest(changed)
        terms = _record_input(_normalize_case(changed))[2]
        proof['monomial_columns'] = len(set().union(*(set(poly) for poly in _record_input(_normalize_case(self.case))[2])))
        proof['new_linear_generators'] = sum(max(map(sum, poly), default=0) == 1 for poly in terms)
        with self.assertRaisesRegex(ValueError, 'outside derived row span'):
            verify_transformation(self.case, changed, proof)

    def test_prime_characteristic_cannot_be_silently_changed(self):
        with self.assertRaisesRegex(ValueError, 'characteristic two'):
            derive_case(self.cases['p83-S-0-f3-fb4-planted'])

    def test_changed_nonpolynomial_semantics_rejected(self):
        changed, proof = derive_case(self.case)
        changed['target'] = copy.deepcopy(changed['target'])
        changed['target']['point'] = None
        proof['derived_case_sha256'] = digest(changed)
        with self.assertRaisesRegex(ValueError, 'non-polynomial'):
            verify_transformation(self.case, changed, proof)


if __name__ == '__main__': unittest.main()
