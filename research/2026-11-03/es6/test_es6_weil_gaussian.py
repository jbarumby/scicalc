"""ES6 tests separate exact-rational certificates and exploratory numeric crosschecks."""
import math
from fractions import Fraction
import mpmath as mp
import pytest
from es6_weil_gaussian import (
 SCALE,ceildiv,pi_bounds,log_int_bounds,unit_exp_bounds,exp_neg_bounds,
 positive_product,sqrt_scaled,inverse_positive,spf_sieve,power_base,
 prime_gaussian_certificate,prime_powers_mpmath,analytic_parts,run_exploration
)

def test_01_ceiling():
 assert ceildiv(10,3)==4 and ceildiv(-10,3)==-3 and ceildiv(0,3)==0

def test_02_pi_enclosure():
 l,u=pi_bounds()
 with mp.workdps(80): assert mp.mpf(l)/SCALE <mp.pi<mp.mpf(u)/SCALE
 assert l>3*SCALE and u<4*SCALE

def test_03_log_enclosure():
 with mp.workdps(80):
  for n in [1,2,3,5,16,47,64,128]:
   l,u=log_int_bounds(n)
   assert mp.mpf(l)/SCALE<=mp.log(n)<=mp.mpf(u)/SCALE

def test_04_exponent_enclosure_extreme():
 with mp.workdps(85):
  for x in [0,1,3,11,100,260]:
   l,u=exp_neg_bounds(x*SCALE,x*SCALE)
   target=mp.exp(-x)
   assert mp.mpf(l)/SCALE<=target<=mp.mpf(u)/SCALE

def test_05_exponent_interval_monotonicity():
 with mp.workdps(85):
  l,u=exp_neg_bounds(2*SCALE,3*SCALE)
  assert mp.mpf(l)/SCALE<=mp.exp(-3)<mp.exp(-2)<=mp.mpf(u)/SCALE

def test_06_positive_product_sqrt_inverse():
 with mp.workdps(85):
  l,u=sqrt_scaled((2*SCALE,2*SCALE));a=inverse_positive((l,u))
  assert mp.mpf(l)/SCALE<=mp.sqrt(2)<=mp.mpf(u)/SCALE
  assert mp.mpf(a[0])/SCALE<=1/mp.sqrt(2)<=mp.mpf(a[1])/SCALE
  c=positive_product(a,(3*SCALE,3*SCALE))
  assert mp.mpf(c[0])/SCALE<=3/mp.sqrt(2)<=mp.mpf(c[1])/SCALE

def test_07_power_sieve():
 spf=spf_sieve(100)
 assert [power_base(n,spf) for n in (2,4,8,9,12,16,17,25,27,36,49,64,81,100)] == [2,2,2,3,None,2,17,5,3,None,7,2,3,None]

def test_08_certified_prime_term_contains_reference():
 cert=prime_gaussian_certificate(64)
 lo,hi=map(int,cert['prime_term_all_n_exact_scaled'])
 assert lo<=hi and hi-lo<10**9
 with mp.workdps(85):
  target=prime_powers_mpmath(500,mp.mpf(1)/50)[0]
  assert mp.mpf(lo)/SCALE<=target<=mp.mpf(hi)/SCALE

def test_09_certification_refines_with_cutoff():
 a=prime_gaussian_certificate(64)
 b=prime_gaussian_certificate(128)
 alo,ahi=map(int,a['prime_term_all_n_exact_scaled'])
 blo,bhi=map(int,b['prime_term_all_n_exact_scaled'])
 assert max(alo,blo)<=min(ahi,bhi)
 assert bhi-blo<=ahi-alo+100

def test_10_reject_unsafe_cutoff():
 with pytest.raises(ValueError):prime_gaussian_certificate(63)
 with pytest.raises(ValueError):run_exploration(J=4)
 with pytest.raises(ValueError):run_exploration(J=7,N=10)

def test_11_guinand_weil_crosscheck_and_differential():
 # NUMERICAL ONLY: no zero-completeness certification, no interval quadrature.
 data=run_exploration(a_values=('0.02',),J=18,N=64,dps=43,T=96)
 c=data['exploration']['cases'][0]
 assert abs(mp.mpf(c['rhs_minus_known_zero_sum_numeric']))<mp.mpf('1e-28')
 assert abs(mp.mpf(c['derivative_discrepancy_numeric']))<mp.mpf('1e-27')
 assert abs(mp.mpf(c['ES4_minus_zeta_known_zero_sum_numeric']))>mp.mpf('1e-4')
 assert 'NON-CERTIFIED' in c['numeric_status']

def test_12_cannot_claim_global_certification():
 r=run_exploration(a_values=(),J=5,N=64,dps=35,T=64)
 assert 'NOT RH proof' in r['scientific_status']
 assert 'without RH' in r['convention']['formula'] if 'without RH' in r['convention']['formula'] else 'no RH assumed' in r['convention']['zero_variable']
 assert 'not interval-certified' in r['exploration']['limitations'][1]
