import unittest,tempfile,json,ast
from pathlib import Path
from fractions import Fraction
from fixed_order_replay import interval,compare,sample_pairs,DEGREES,HERE
from verify_fixed import run
class TestFixedOrder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.b=json.loads((HERE/'sample_reference.json').read_text())
    def test_modes_preserved(self):
        self.assertEqual(tuple(sum(v) for v in DEGREES.values()),(24,26,44))
    def test_sample_counts(self):
        self.assertEqual([len(self.b['entries'][w]) for w in DEGREES],[8,8,14])
    def test_sample_partition(self):
        for w in DEGREES:
            a=sample_pairs(self.b,w,0);b=sample_pairs(self.b,w,1)
            self.assertFalse(set(a)&set(b))
            self.assertEqual(len(a)+len(b),len(self.b['entries'][w]))
    def test_fixed_not_adaptive(self):
        s=(Path(__file__).parent/'fixed_order_replay.py').read_text()
        ast.parse(s)
        self.assertIn('n=110',s)
        self.assertIn('M=arb(600)',s)
        self.assertNotIn('while not bool(err.upper()<target)',s)
        self.assertIn('arb.legendre_p_root',s)
        self.assertNotIn('acb.integral(',s)
    def test_exact_binary_endpoints(self):
        self.assertEqual(interval(['3','-1','1','-2']),(Fraction(5,4),Fraction(7,4)))
    def test_disjoint_detected(self):
        self.assertEqual(compare(['0','0','1','-3'],['1','0','1','-3']),'DISJOINT')
    def test_fail_closed_missing(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(run(d)['status'],'FAIL_CLOSED')
if __name__=='__main__':unittest.main()
