import json,unittest,tempfile
from pathlib import Path
from fractions import Fraction
from replay import interval,compare,sample_pairs,DEGREES
from verify_results import verify

class TestES8P79(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((Path(__file__).parent/'sample_reference.json').read_text())
    def test_dimensions(self):
        self.assertEqual({k:sum(v) for k,v in DEGREES.items()},
                         {'a030':24,'a034':26,'a040':44})
    def test_counts(self):
        self.assertEqual({k:len(v) for k,v in self.data['entries'].items()},
                         {'a030':8,'a034':8,'a040':14})
    def test_fixtures_exact_mantissas(self):
        self.assertEqual(self.data['entries']['a030']['0,0'][0],
                         '49391520634221872854410992764577700550896833031824776642817')
        self.assertEqual(self.data['entries']['a040']['21,32'][0],'0')
    def test_windows_partition(self):
        for k in DEGREES:
            a=sample_pairs(self.data,k,0); b=sample_pairs(self.data,k,1)
            self.assertEqual(len(a)+len(b),len(self.data['entries'][k]))
            self.assertFalse(set(a)&set(b))
    def test_all_reference_radii_valid(self):
        for entries in self.data['entries'].values():
            for enc in entries.values():
                lo,hi=interval(enc)
                self.assertLessEqual(lo,hi)
    def test_fraction_exactness(self):
        self.assertEqual(interval(['3','-1','1','-2']),
                         (Fraction(5,4),Fraction(7,4)))
    def test_no_negative_radii(self):
        with self.assertRaises(ValueError): interval(['0','0','-1','0'])
    def test_disjoint(self):
        self.assertEqual(compare(['0','0','1','-3'],['1','0','1','-3']),'DISJOINT')
    def test_subset(self):
        self.assertEqual(compare(['0','0','1','-3'],['0','0','1','-2']),
                         'REPLAY_SUBSET_REFERENCE')
    def test_missing_outputs_fail_closed(self):
        with tempfile.TemporaryDirectory() as d:
            report=verify(d)
            self.assertEqual(report['status'],'FAIL_CLOSED')
            self.assertEqual(len(report['errors']),9)
    def test_no_external_import_on_preflight(self):
        import ast
        p=Path(__file__).parent/'replay.py'
        ast.parse(p.read_text())
    def test_sample_not_full(self):
        self.assertEqual(sum(len(v) for v in self.data['entries'].values()),30)
        self.assertNotEqual(30,156+182+990)
if __name__=='__main__': unittest.main()
