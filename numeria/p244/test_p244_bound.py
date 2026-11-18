import shutil
import tempfile
import unittest
from pathlib import Path

from p244_replay_bound import PACKAGE, run_bound, verify_static_pins
from p244_replay_fail_closed import validate_replays

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / 'p237' / 'data'
CERT = PACKAGE / 'P242_AUDIT_CERT.json'


class P244BoundTests(unittest.TestCase):
    def test_static_pins(self):
        pins = verify_static_pins(CERT)
        self.assertEqual(pins['verified_exporter_sha256'], '58be39fcfd6b424cdec9af9b20b5a3b9d78e6fa3f0a40d4127916f8c97b92fac')
        self.assertEqual(pins['verified_p242_certificate_sha256'], '9f89c17339c9daca7d1e2e020dbfdbecc83fc47a79301b6ebde544dad0099447')

    def test_tampered_inverse_bound_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            work = Path(temp)
            for sector in ('even', 'odd'):
                name = f'P242_PIVOTS_{sector}_1024.txt'
                shutil.copy2(PACKAGE / name, work / name)
            path = work / 'P242_PIVOTS_even_1024.txt'
            lines = path.read_text().splitlines()
            lines = [('# inverse_frobenius_squared_upper 1 0' if line.startswith('# inverse_frobenius_squared_upper ') else line) for line in lines]
            path.write_text('\n'.join(lines) + '\n')
            with self.assertRaisesRegex(ArithmeticError, 'regenerated trace SHA mismatch'):
                validate_replays(work, PACKAGE, CERT)

    def test_check_existing_reports_scope_and_full_matrix_margin(self):
        with tempfile.TemporaryDirectory() as temp:
            work = Path(temp)
            for sector in ('even', 'odd'):
                name = f'P242_PIVOTS_{sector}_1024.txt'
                shutil.copy2(PACKAGE / name, work / name)
            result = run_bound(DATA, PACKAGE, CERT, work, recompute=False)
            self.assertFalse(result['recomputed_using_compiled_MPFR_C_in_this_run'])
            self.assertEqual(result['trace_check_mode'], 'preexisting_trace_check_only')
            self.assertEqual(result['positive_pivots'], 754)
            self.assertEqual(result['margin_lower_exact'], '1/10^18')
            for sector in ('even', 'odd'):
                self.assertTrue(result['MPFR_traces'][sector]['trace_digest_and_bytes_match'])
                self.assertNotIn('identical_reexecution', result['MPFR_traces'][sector])


if __name__ == '__main__':
    unittest.main(verbosity=2)
