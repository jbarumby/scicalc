import json
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / 'P231_PREFIX4_DATA.json').read_text())


def matrix(key):
    return [[F(x) for x in row] for row in DATA[key]]


class P231PrefixTests(unittest.TestCase):
    def test_provenance_and_p230_positive_margins(self):
        self.assertEqual(DATA['sector'], 'even')
        self.assertEqual(DATA['matrix_dimension'], 4)
        self.assertEqual(DATA['P193_center_gzip_sha256'],
                         'c105996b75a5782de338dc5b9711b97839f7093f190367bfb67e43cc418b0c65')
        self.assertGreater(F(DATA['P230_center_spd_lower_bound_exact']), 0)
        self.assertGreater(F(DATA['P230_200_mode_target_lambda_lower_bound_exact']), 0)
        self.assertEqual(DATA['checks']['all_16_reflected_boxes_within_epsilon'], True)

    def test_all_sixteen_p227_boxes_fit_coarse_epsilon(self):
        eps = F(DATA['entrywise_epsilon_exact'])
        rows = DATA['P227_target_boxes_exact_upper_triangle_and_reflection']
        self.assertEqual(len(rows), 16)
        self.assertEqual({(int(r[0]), int(r[1])) for r in rows},
                         {(i, j) for i in range(4) for j in range(4)})
        for i, j, lo, hi, c, dev in rows:
            lo, hi, c, dev = map(F, (lo, hi, c, dev))
            self.assertLessEqual(lo, hi)
            self.assertEqual(dev, max(abs(lo-c), abs(hi-c)))
            self.assertLessEqual(dev, eps)

    def test_center_is_symmetric_and_shifted_ldl_is_exact(self):
        C = matrix('rounded_center_matrix_exact')
        L = matrix('shifted_LDL_unit_lower_exact')
        d = [F(x) for x in DATA['shifted_LDL_diagonal_exact']]
        shift = F(DATA['center_shift_exact'])
        self.assertEqual(len(C), 4)
        self.assertEqual(len(L), 4)
        self.assertEqual(len(d), 4)
        for i in range(4):
            self.assertGreater(d[i], 0)
            for j in range(4):
                self.assertEqual(C[i][j], C[j][i])
                target = C[i][j] - (shift if i == j else 0)
                reconstructed = sum(L[i][k] * d[k] * L[j][k] for k in range(4))
                self.assertEqual(target, reconstructed)
        self.assertEqual(F(DATA['quadratic_error_rho_exact']),
                         4 * F(DATA['entrywise_epsilon_exact']))
        self.assertLess(F(DATA['quadratic_error_rho_exact']), shift)

    def test_generated_lean_matches_data_and_has_no_placeholders(self):
        source = (ROOT.parent / 'NUMERIA_E8K_FORMAL' / 'P231_P230Prefix4.lean').read_text()
        import hashlib
        digest = hashlib.sha256((ROOT / 'P231_PREFIX4_DATA.json').read_bytes()).hexdigest()
        self.assertIn(digest, source)
        for marker in ('sorry', 'admit', 'axiom'):
            self.assertNotIn(marker, source)
        for theorem in ('factor_entries4', 'pivots_pos4', 'center_lower4', 'matrix_bridge', 'P231_box_positive'):
            self.assertIn(theorem, source)


if __name__ == '__main__':
    unittest.main()
