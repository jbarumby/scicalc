"""ES8-P73 preflight tests. Does not require python-flint."""
import json
import unittest
from pathlib import Path
from fractions import Fraction
from replay_arb import endpoints, relation, selected_pairs

HERE=Path(__file__).parent

class TestReplayConfiguration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=json.loads((HERE/"historical_sample.json").read_text())
    def test_source_provenance(self):
        self.assertEqual(self.data["source_file_sha256"],"46457041797e5489745b8e34288fcb13c03ab194e566c50ec2b91d05e8b3a765")
    def test_dimensions(self):
        self.assertEqual(self.data["degrees"],[14,22,14])
        self.assertEqual(sum(self.data["degrees"]),50)
    def test_19_pairs_complete(self):
        self.assertEqual(len(self.data["pairs"]),19)
        self.assertEqual(len(self.data["entries"]),19)
    def test_shards_partition(self):
        samples=[selected_pairs(self.data,k) for k in range(4)]
        self.assertEqual(sum(map(len,samples)),19)
        self.assertEqual(len(set(p for group in samples for p in group)),19)
    def test_all_historical_balls_finite(self):
        for key, enc in self.data["entries"].items():
            lo,hi=endpoints(enc)
            self.assertLessEqual(lo,hi,key)
    def test_point_and_radius(self):
        self.assertEqual(endpoints([3,-1,1,-2]),(Fraction(5,4),Fraction(7,4)))
    def test_overlap(self):
        self.assertEqual(relation([0,0,1,0],[1,0,1,0]),"OVERLAP")
    def test_disjoint(self):
        self.assertEqual(relation([0,0,1,-3],[1,0,1,-3]),"DISJOINT")
    def test_subset(self):
        self.assertEqual(relation([0,0,1,-3],[0,0,1,-2]),"REPLAY_SUBSET_HISTORICAL")
    def test_no_negative_radii(self):
        with self.assertRaises(ValueError): endpoints([0,0,-1,0])

if __name__=="__main__": unittest.main()
