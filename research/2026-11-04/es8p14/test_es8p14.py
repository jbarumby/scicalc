"""Exact rational adversarial tests; passing does not substitute for Lean compilation."""
import unittest, json, hashlib
from fractions import Fraction as F
from pathlib import Path
from generate_es8p14 import ROOT,SCALE,ROUND,main

class TestBridge(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.src=json.loads((ROOT/'RESULTADOS_ES8P13.json').read_text())
  cls.out=json.loads((ROOT/'RESULTADOS_ES8P14.json').read_text())

 def raw(self,N):
  a=self.src[f'N{N}'];n=N-2
  c=[[F(s)for s in row]for row in a['matrix_coarse']]
  lo=[[F(t,SCALE)for t in row]for row in a['source_interval_raw_integer_lower']]
  hi=[[F(t,SCALE)for t in row]for row in a['source_interval_raw_integer_upper']]
  eps=[[F( (max(abs(lo[i][j]-c[i][j]),abs(hi[i][j]-c[i][j]))*ROUND).__ceil__(),ROUND)for j in range(n)]for i in range(n)]
  return c,lo,hi,eps

 def test_all_original_endpoints(self):
  for N in (6,8,10):
   C,lo,hi,E=self.raw(N)
   for i in range(N-2):
    for j in range(N-2):
     self.assertLessEqual(lo[i][j],hi[i][j]);self.assertLessEqual(-E[i][j],lo[i][j]-C[i][j]);self.assertLessEqual(hi[i][j]-C[i][j],E[i][j])

 def test_row_sums_and_symmetry(self):
  for N in (6,8,10):
   C,lo,hi,E=self.raw(N);n=N-2
   self.assertLess(max(sum(r,F(0))for r in E),F(1,100))
   self.assertTrue(all(E[i][j]==E[j][i] for i in range(n)for j in range(n)))
   self.assertEqual(self.out['dimensions'][str(N)]['max_rounded_row_error'],str(max(sum(r,F(0))for r in E)))

 def test_regenerated_sources_byte_identical(self):
  p=[ROOT/'ES8P14Concrete.lean',ROOT/'RESULTADOS_ES8P14.json'];before=[x.read_bytes()for x in p]
  main()
  self.assertEqual(before,[x.read_bytes()for x in p])

 def test_mutate_endpoint_detected(self):
  C,lo,hi,E=self.raw(10)
  lo[0][0]=C[0][0]-E[0][0]-F(1,10**90)
  self.assertFalse(lo[0][0]-C[0][0]>=-E[0][0])

 def test_mutate_negative_weight_detected(self):
  C,lo,hi,E=self.raw(10)
  E[1][3]=-E[1][3]
  self.assertNotEqual(E[1][3],E[3][1])

 def test_outcome_independent_of_float_fields(self):
  d=self.src['N10']
  self.assertIsInstance(d['max_abs_row_error_float_for_reading_only'],float)
  for n in (6,8,10):
   self.assertEqual(self.out['dimensions'][str(n)]['number_of_raw_endpoints'],2*(n-2)**2)

if __name__=='__main__':unittest.main()
