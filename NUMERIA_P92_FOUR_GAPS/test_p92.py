"""Exact P92 theorem checks; no floating point, including with python -O."""
import unittest
from fractions import Fraction as Q
from itertools import product
from .certificate import (
    SCHEMA, ENVELOPES, ENDPOINT_A, ENDPOINT_D,
    certify, d4, interior, dot, region_planes, region_vertices,
    tetrahedra, numerator, peval, check_barycentric_identity,
)

class P92FourGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=certify(full=False)

    def test_full_census_and_scope(self):
        c=self.cert
        self.assertEqual(c["schema"],SCHEMA)
        self.assertEqual((c["region_count"],c["tetra_count"],
                          c["positive_exact_coefficients"],
                          c["negative_coefficients"],c["half_simplex_volume"]),
                         (35,246,133363,0,"1/12"))
        self.assertTrue(c["all_four_gap_budget_proved"])
        for key in ("all_m_ge_5_budget_proved","full_Preissmann_Lemma5_proved",
                    "full_Preissmann_Lemma6_proved","Riemann_Hypothesis_proved"):
            self.assertFalse(c[key],key)

    def test_uniform_exact_sharp_equality(self):
        g=(Q(1,4),)*4
        self.assertEqual(interior(g),Q(164,9))
        self.assertEqual(d4(Q(1,4)),Q(82,9))
        self.assertEqual(2*d4(Q(1,4)),interior(g))
        self.assertEqual(d4(Q(1,3)),Q(27,4))
        self.assertEqual(d4(Q(1,2)),Q(4))
        self.assertEqual(d4(Q(1)),Q(0))

    def test_chord_ranges_join_exactly(self):
        for j in range(4):
            alpha,beta=ENVELOPES[j]
            lo,hi=ENDPOINT_D[j:j+2]
            for s in (lo,hi):
                if s>0:
                    self.assertEqual(d4(s),alpha/s+beta-3*s)

    def test_coverage_independent_sum(self):
        cells=0
        vol=Q(0)
        num_tetra=0
        for bits in product((0,1),repeat=3):
            for i in range(3):
                for j in range(4):
                    planes=region_planes(bits,i,j)
                    verts=region_vertices(planes)
                    pieces=list(tetrahedra(planes,verts))
                    if pieces:
                        cells+=1
                    for tet,size in pieces:
                        num_tetra+=1
                        self.assertGreater(size,0)
                        for v in tet:
                            self.assertTrue(all(dot(p,v)>=0 for p in planes))
                        vol+=size
        self.assertEqual((cells,num_tetra,vol),(35,246,Q(1,12)))

    def test_every_exponent_and_identity_in_summary(self):
        for cell in self.cert["regions"]:
            bits=tuple(map(int,cell["min_branches"]))
            i,j=cell["a_interval"],cell["d_interval"]
            for row in cell["tetrahedra"]:
                tet=tuple(tuple(map(Q,v)) for v in row["vertices_a_b_d"])
                p=numerator(bits,i,j,tet)
                self.assertEqual(row["positive_monomials"],len(p))
                self.assertTrue(p and all(c>0 for c in p.values()))
                self.assertTrue(all(sum(e)==13 for e in p))
                check_barycentric_identity(p,tet)
                self.assertEqual(str(Q(row["volume"])),row["volume"])

    def test_rational_finite_regressions_do_not_replace_proof(self):
        for n in range(4,13):
            for a in range(1,n-2):
                for b in range(1,n-a-1):
                    for c in range(1,n-a-b):
                        d=n-a-b-c
                        g=tuple(Q(k,n) for k in (a,b,c,d))
                        self.assertGreaterEqual(d4(g[0])+d4(g[-1])-interior(g),0)

    def test_invalid_gaps(self):
        with self.assertRaises(ValueError):
            d4(0)
        with self.assertRaises(ValueError):
            d4(Q(3,2))
        for g in ((Q(1,4),Q(1,4),Q(1,4),Q(0)),
                  (Q(1),Q(1),Q(1),Q(1))):
            with self.assertRaises(ValueError):interior(g)

if __name__=="__main__":
    unittest.main()
