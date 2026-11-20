"""Fail-closed P255 checks. No Lean kernel/compiler is invoked."""
import hashlib
import tempfile
import unittest
from pathlib import Path
from .audit_statements import audit, check_source, check_headlines, check_mathlib, PINNED, EXPECTED, require

ROOT=Path(__file__).with_name('fixtures')

class P255Tests(unittest.TestCase):
    def test_authenticated_three_source_hashes(self):
        for name,digest in PINNED.items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),digest)

    def test_exactly_five_named_types(self):
        r=audit(ROOT)
        self.assertEqual({x['name']:(x['type'],x['depends_immediately_on']) for x in r['five_headline_types_verified_from_exact_pinned_sources']},EXPECTED)

    def test_explicit_limitation_flags(self):
        r=audit(ROOT)
        self.assertFalse(r['Lean_kernel_invoked_in_this_audit'])
        self.assertFalse(r['solution_full_export_replayed_in_this_audit'])
        self.assertFalse(r['RH_proved_independently_by_NUMERIA'])
        self.assertTrue(r['theorem_statement_verification_is_not_a_proof'])

    def test_title_was_not_redefined(self):
        main=(ROOT/'Main.lean').read_text()
        self.assertNotIn('def RiemannHypothesis',main)
        self.assertIn('theorem riemannHypothesis : RiemannHypothesis := rhOfAll',main)

    def test_pinned_modification_rejected(self):
        raw=(ROOT/'Main.lean').read_bytes()
        with self.assertRaisesRegex(ValueError,'SHA-256'):
            check_source('Main.lean',raw+b'\n')

    def test_statement_weakening_rejected(self):
        raw=(ROOT/'Statements.lean').read_bytes()
        with self.assertRaisesRegex(ValueError,'SHA-256'):
            check_source('Statements.lean',raw.replace(b'0 < s.re',b'1/4 < s.re'))

    def test_main_wrong_type_rejected(self):
        x=(ROOT/'Main.lean').read_text().replace('theorem riemannHypothesis : RiemannHypothesis',
                                                   'theorem riemannHypothesis : True')
        with self.assertRaisesRegex(ValueError,'types'):
            check_headlines(x,(ROOT/'Statements.lean').read_text(),(ROOT/'MainTheorems.lean').read_text())

    def test_missing_allchi_rejected(self):
        x=(ROOT/'MainTheorems.lean').read_text().replace('everyCharacter','everyCharacterFake')
        with self.assertRaisesRegex(ValueError,'all-character'):
            check_headlines((ROOT/'Main.lean').read_text(),(ROOT/'Statements.lean').read_text(),x)

    def test_modified_mathlib_byte_rejected(self):
        with self.assertRaisesRegex(ValueError,'blob mismatch'):
            check_mathlib(b'def RiemannHypothesis : Prop := True\n')

    def test_empty_input_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                audit(Path(directory))

    def test_require_not_python_assert(self):
        with self.assertRaisesRegex(ValueError,'adversarial'):
            require(False,'adversarial')

if __name__=='__main__':unittest.main()
