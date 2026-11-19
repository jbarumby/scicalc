#!/usr/bin/env python3
"""ES8-P237: rational verification of the E3 symbolic majorant.

The non-rational mathematical inputs are standard classical identities:
psi(z)=-gamma-1/z+sum(z/[n(n+z)]) on Re(z)>0;
conjugation psi(zbar)=conj(psi(z)); the harmonic-number integral bound;
gamma<1; modulus and triangle inequalities; and cosh-series expansion.
The script checks every *finite rational* inequality in the supplied proof.
No numerical sampling, Arb, float-based acceptance, or RH assertion.
"""
from fractions import Fraction as F
from math import factorial
import argparse,json
from p237_beta_tight import beta_tight
from functools import lru_cache

_cached_beta_tight = lru_cache(maxsize=1)(beta_tight)
from pathlib import Path


def exp_lower(x,n):
    x=F(x)
    return sum((x**k/F(factorial(k)) for k in range(n+1)),F(0))


def need(condition, description):
    if not condition: raise ArithmeticError(description)

@lru_cache(maxsize=1)
def check_certificate():
    # E3 panel geometry, extracted from exact affine mapping of [-1,1].
    half_axis_real = F(1,4)*(3+F(1,3))/2
    half_axis_imag = F(1,4)*(3-F(1,3))/2
    t_real_max = F(2*1007+1,4)+half_axis_real
    t_imag_max = half_axis_imag
    z_real_min=F(1,4)-t_imag_max/2
    z_real_max=F(1,4)+t_imag_max/2
    z_imag_max=t_real_max/2
    z_mod_upper=z_real_max+z_imag_max
    need(half_axis_real==F(5,12) and half_axis_imag==F(1,3), "rational gate failure")
    need(t_real_max==F(3025,6) and t_imag_max==F(1,3), "rational gate failure")
    need(z_real_min==F(1,12) and z_real_max==F(5,12), "rational gate failure")
    need(z_mod_upper==F(505,2) and z_mod_upper<253, "rational gate failure")
    # Exponential lower bounds yield logarithm upper bounds by monotonicity.
    need(exp_lower(1,2)>2, "exp(1)>2")  # log(2)<1
    need(exp_lower(2,3)>5, "exp(2)>5")  # log(5),log(4),log(3)<2
    need(exp_lower(6,7)>253, "exp(6)>253")  # log(253)<6
    # An additional *entirely rational* check of the H_253 bound,
    # independent of the integral comparison and log(253)<6 argument.
    H_253=sum((F(1,n) for n in range(1,254)),F())
    need(H_253 < 7, "harmonic-number bound failed")
    # pi<4 is a classical external fact; log pi <log4<2 follows above.
    # For positive Re z and |z|<253: |psi(z)|< gamma +1/|z|+H_253+253*sum(n^-2,n>253).
    # gamma<1, |z|>=1/12, H_253<=1+log(253)<7, tail<253/253=1.
    digamma_upper=F(1)+F(12)+F(7)+F(1)
    need(digamma_upper==21, "rational gate failure")
    # cosh(2/3) <= sum (2/3)^(2k)/2^k = 9/7 (using (2k)!>=2^k).
    cosh_majorant=F(1,1)/(F(1)-F(2,9))
    need(cosh_majorant==F(9,7), "rational gate failure")
    # A_L =2 log 2/sqrt2 + 2 log3/sqrt3 + log2 + 2log5/sqrt5 < 2+4+1+4 = 11.
    prime_coeff_majorant=F(11)
    prime_sum_majorant=cosh_majorant*prime_coeff_majorant
    # Certified independently using exactly rational Machin/atanh/sqrt code.
    beta_lo,beta_hi,_,_=_cached_beta_tight()
    need(0 < beta_lo < beta_hi < F(2), "beta rational certificate not positive/<2")
    final=digamma_upper+2+prime_sum_majorant+2
    need(final==F(274,7) and final<40 and final<50, "rational gate failure")
    return {'schema':'NUMERIA.ES8P237.SYMBOL_E3_RATIONAL_GATE.1',
      'status':'finite_rational_inequalities_verified_conditional_on_named_classical_identities',
      'L':'161/200','Tsharp':504,'panels':1008,'rho':3,
      'ellipse_t_real_bound':str(t_real_max),'ellipse_t_imag_bound':str(t_imag_max),
      'digamma_arguments_real_min':str(z_real_min),'digamma_arguments_real_max':str(z_real_max),
      'digamma_arguments_modulus_strict_less':253,
      'finite_exp_inequality_witnesses':{
        'exp1_sum_through_2':str(exp_lower(1,2)),
        'exp2_sum_through_3':str(exp_lower(2,3)),
        'exp6_sum_through_7':str(exp_lower(6,7))},
      'digamma_bound_strict_less':str(digamma_upper),
      'cosh_2over3_bound':str(cosh_majorant),
      'prime_coefficient_sum_bound_strict_less':str(prime_coeff_majorant),
      'beta_bound': '0 < beta_504 < 2; recalculated exactly via beta_tight()',
      'beta_rational_lo_gt_zero':beta_lo>0,'beta_rational_hi_lt_two':beta_hi<2,
      'symbol_minus_beta_bound_strict_less':str(final),
      'supports_GL60_factor_50_conditionally':True,
      'classical_inputs_not_formalized':['digamma partial fractions on Re(z)>0',
        'digamma conjugation', 'Euler gamma < 1', 'H_n <= 1+log n',
        'pi<4', 'cosh power series and (2k)! >= 2^k',
        'beta bound from P237 independent rational certificate'],
      'Weil_Zhu_identity_verified':False,'global_LDL_done':False,'RH_proved':False}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args();cert=check_certificate();a.out.write_text(json.dumps(cert,indent=2,ensure_ascii=False)+'\n')
    print('PASS E3 symbol bound |Psi-beta|<274/7<40<50 (standard digamma identity imported)')

if __name__=='__main__': main()
