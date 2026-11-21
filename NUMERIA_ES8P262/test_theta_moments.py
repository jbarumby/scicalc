"""P262: positive theta moments, certified tails, and fail-closed KPS sign logic."""
from fractions import Fraction as Q
import unittest
from flint import arb, ctx

from .theta_moments import (
    U_END, M_MAX, SOURCE, need, tail_rational, endpoints,
    moment, paper_G, paper_phi,
)

class ThetaMomentsTests(unittest.TestCase):
    def test_k_zero_positive_tail(self):
        self.assertGreater(tail_rational(0),0)
        self.assertLess(tail_rational(0),Q(1,10**14))
    def test_k_2_tail_small(self):
        self.assertTrue(0<tail_rational(2)<Q(1,10**15))
    def test_k_8_tail_small(self):
        self.assertTrue(0<tail_rational(8)<Q(1,10**14))
    def test_invalid_k_fails(self):
        for k in (-1,31,1.5):
            with self.assertRaises(ValueError):tail_rational(k)
    def test_invalid_precision_fails(self):
        with self.assertRaises(ValueError):moment(0,precision_bits=64)
        with self.assertRaises(ValueError):moment(1,precision_bits=192)
    def test_source_is_exact_pdf(self):
        self.assertEqual(U_END,5)
        self.assertEqual(M_MAX,3)
        self.assertIn("aif.3600.pdf",SOURCE)
    def test_G_algebra_exact(self):
        with ctx.workprec(192):
            x={0:arb(1),2:arb(1),4:arb(1)}
            self.assertEqual(paper_G(1,x),arb(63))
            self.assertEqual(paper_G(2,x),arb(383))
    def test_negative_printed_phi_for_positive_G(self):
        with ctx.workprec(192):
            g={1:arb(63),2:arb(383)}
            self.assertLess(paper_phi(1,g),0)
            self.assertTrue(paper_phi(1,g).overlaps(-arb(63)/(16*383)))
    def test_kps_printed_gate_not_equivalent_to_RH(self):
        g={1:arb(63),2:arb(383)}
        self.assertTrue(paper_phi(1,g)<0)
        with self.assertRaises(ArithmeticError):
            need(paper_phi(1,g)>0,'Bernstein positivity would be violated')
    def test_dyadic_endpoint_schema(self):
        e=endpoints(arb(1)/arb(3))
        self.assertEqual(len(e),2)
        self.assertEqual(len(e[0]),2)
        self.assertTrue(all(type(x[1]) is int for x in e))
    def test_moment_F0_positive(self):
        a, err=moment(0,precision_bits=160)
        self.assertTrue(a>0)
        self.assertLess(err,Q(1,10**14))
        self.assertTrue(a<1)
    def test_moment_F2_positive(self):
        a,err=moment(2,precision_bits=160)
        self.assertTrue(0<a<1)
        self.assertLess(err,Q(1,10**14))

if __name__=="__main__":unittest.main()
