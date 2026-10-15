import copy
import unittest

from ecbench.cli import catalog, select_algebra_cases
from ecbench.fixtures import verify_cases


class FixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.suite=catalog('algebra_cases')

    def test_original_prime_and_binary_equations_match_signed_truth(self):
        selection=select_algebra_cases(self.suite,case_ids=[
            'p83-S-0-f3-fb4-planted','K83-a0-f3-low_monomials-d2-planted'])
        result=verify_cases(selection)
        self.assertEqual(result['cases_checked'],2)
        self.assertEqual([row['ordered_relations'] for row in result['results']],[1,4])

    def test_changed_polynomial_cannot_keep_old_golden_truth(self):
        selection=copy.deepcopy(select_algebra_cases(self.suite,['p83-S-0-f3-fb4-planted']))
        case=selection['cases'][0]
        case['polynomials'].append({'id':'tampered_unit','role':'invalid_fixture',
            'terms':[{'coefficient':'1','exponents':[0]*len(case['variables'])}]})
        with self.assertRaisesRegex(ValueError,'Fixture truth mismatch'):
            verify_cases(selection)

    def test_selection_cannot_silently_ignore_unknown_case(self):
        with self.assertRaisesRegex(ValueError,'Unknown polynomial case'):
            select_algebra_cases(self.suite,['nonexistent'])


if __name__=='__main__':unittest.main()
