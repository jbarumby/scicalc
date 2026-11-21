"""Small-geometry and adversarial tests; heavy matrix is in separate workflow step."""
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from flint import arb_mat
from .gate_b import (SOURCE_SHA, CUTOFF2, FLOOR2, validate_parameters,
                     check_p258_archive, prefix, rational_negative_witness)

class GateBTests(unittest.TestCase):
    def test_parameters(self):
        self.assertEqual(validate_parameters()['log_850_504_lower'],'346/677')
        self.assertEqual((CUTOFF2,FLOOR2),(850,F(1,2)))

    def test_p258_sha_is_fixed(self):
        self.assertEqual(len(SOURCE_SHA),64)

    def test_reject_modified_source(self):
        with TemporaryDirectory() as folder:
            path=Path(folder)/'bad.json'
            path.write_text('{"schema":"fake"}\n')
            with self.assertRaisesRegex(ArithmeticError,'SHA'):
                check_p258_archive(path)

    def test_prefix_stable(self):
        a=arb_mat([[1,2,3],[2,4,5],[3,5,8]])
        q=prefix(a,2)
        self.assertEqual((q.nrows(),q.ncols()),(2,2))
        self.assertEqual((q[0,1],q[1,1]),(2,4))
        with self.assertRaises(ValueError):
            prefix(a,4)

    def test_negative_detection_does_not_claim_unproved_extension(self):
        matrix=arb_mat([[1,0],[0,-1]])
        witness=rational_negative_witness(matrix,1,{'analytical':F(0)})
        self.assertIsNotNone(witness)
        self.assertEqual(witness['scope'],'negative at exact P258 rational halfwidth only')
        self.assertEqual(witness['exact_rational_replay']['witness_strict_sign'],'negative')

    def test_no_false_negative_on_positive(self):
        m=arb_mat([[2,0],[0,1]])
        self.assertIsNone(rational_negative_witness(m,None,{}))

if __name__ == '__main__': unittest.main()
