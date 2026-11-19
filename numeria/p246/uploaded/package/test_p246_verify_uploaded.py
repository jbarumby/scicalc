import json
import shutil
import tempfile
import unittest
from pathlib import Path
from p246_verify_uploaded import PINS, verify

SOURCE_DIR=Path(__file__).resolve().parent/'sources'

class TestP246(unittest.TestCase):
    def test_original_uploads_verify(self):
        d=verify(SOURCE_DIR)
        self.assertEqual(d['archived_p245d_inclusions_verified_rationally'],78)
        self.assertEqual(d['archived_high_degree_coarse_cases'],18)
        self.assertEqual(d['p237_analytic_amplitudes_independently_reconstructed'],754)
        self.assertFalse(d['p245d_gl60_contributions_reexecuted_from_integrand'])

    def mutate(self, filename, mutator, bypass_hash=False):
        with tempfile.TemporaryDirectory() as td:
            dest=Path(td)
            for name in PINS:
                shutil.copy2(SOURCE_DIR/name,dest/name)
            p=dest/filename
            obj=json.loads(p.read_text())
            mutator(obj)
            p.write_text(json.dumps(obj))
            with self.assertRaises((ArithmeticError,ValueError,KeyError)):
                verify(dest,check_hashes=not bypass_hash)

    def test_reject_modified_hash(self):
        self.mutate('P245D_SELECTED_25_50.json', lambda d:d.update({'RH_proved':True}))

    def test_reject_fake_stored_enclosure_even_without_hash(self):
        def mutate(d):
            r=d['comparisons'][0]['selected_entry_comparisons'][0]['replay_interval']
            r['lower']=str(int(r['upper'])+1)
        self.mutate('P245D_SELECTED_500_525.json',mutate,bypass_hash=True)

    def test_reject_fake_ball_fraction_even_without_hash(self):
        def mutate(d):
            r=d['comparisons'][0]['selected_entry_comparisons'][0]
            r['source_interval_lower_fraction']='0'
        self.mutate('P245D_SELECTED_975_1000.json',mutate,bypass_hash=True)

    def test_reject_wrong_analytic_amplitude_even_without_hash(self):
        def mutate(d):
            d['sectors']['even']['max_basis_bound_exact']='83'
        self.mutate('P237_GL60_ANALYTIC_CERT.json',mutate,bypass_hash=True)

    def test_reject_wrong_analytic_geometry_even_without_hash(self):
        def mutate(d):
            d['t_imag_abs_upper']='1/2'
        self.mutate('P237_ANALYTIC_GATES.json',mutate,bypass_hash=True)

if __name__=='__main__':unittest.main()
