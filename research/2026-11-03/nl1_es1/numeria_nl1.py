#!/usr/bin/env python3
"""NUMERIA-NL1/ES1, 2026-11-03.

This program has two deliberately distinct layers:
  (A) EXACT rational interval Krawczyk certificate for a 2x2 illustrative
      self-adjoint spectral pencil; this is a *true finite model* certificate.
  (B) 2D numerical Newton for actual zeta zeros with mpmath; the reported
      residuals are *NOT* interval-certified, nor a RH proof.

Python >=3.10; mpmath 1.3.0; Python stdlib.
"""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path
import platform
from typing import Iterable
import mpmath as mp


@dataclass(frozen=True)
class Interval:
    lo: Q
    hi: Q

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError('invalid interval')

    @staticmethod
    def point(x):
        x = Q(x)
        return Interval(x, x)

    def __add__(self, other):
        other = as_iv(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -as_iv(other)

    def __rsub__(self, other):
        return as_iv(other) - self

    def __mul__(self, other):
        other = as_iv(other)
        a = [x*y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Interval(min(a), max(a))

    __rmul__ = __mul__

    def strictly_inside(self, other):
        return other.lo < self.lo and self.hi < other.hi

    def maxabs(self):
        return max(abs(self.lo), abs(self.hi))

    def to_json(self):
        return {'lo': str(self.lo), 'hi': str(self.hi),
                'lo_approx': float(self.lo), 'hi_approx': float(self.hi)}


def as_iv(x):
    return x if isinstance(x, Interval) else Interval.point(x)


def spectral_F(a, b):
    # Eigenvalue zero and trace fixed for H(a,b)=[[a,1],[1,b]].
    return [a*b-Q(1), a+b-Q(5, 2)]


def spectral_J(a, b):
    return [[b, a], [as_iv(1), as_iv(1)]]


def invert_2x2(m):
    det = m[0][0]*m[1][1] - m[0][1]*m[1][0]
    if not det:
        raise ValueError('singular center Jacobian')
    return [[m[1][1]/det, -m[0][1]/det],
            [-m[1][0]/det, m[0][0]/det]]


def krawczyk_2d(center, width):
    """Sufficient strict Krawczyk inclusion and contraction via EXACT Q.

    K(X)=c-C F(c)+(I-C J(X))(X-c), C=J(c)^{-1}.
    A strict inclusion + infinity contraction <1 implies a unique zero in X.
    """
    c = [Q(z) for z in center]
    w = [Q(z) for z in width]
    if any(z <= 0 for z in w):
        raise ValueError('widths must be positive')
    X = [Interval(z-h, z+h) for z,h in zip(c,w)]
    C = invert_2x2([[Q(c[1]),Q(c[0])],[Q(1),Q(1)]])
    fc = spectral_F(*c)
    jx = spectral_J(*X)
    a = [[Interval.point(int(i==j)) - sum((C[i][k] * jx[k][j] for k in range(2)), Interval.point(0))
         for j in range(2)] for i in range(2)]
    k = [Interval.point(c[i]-sum(C[i][z]*fc[z] for z in range(2)))+
         sum((a[i][j]*(X[j]-c[j]) for j in range(2)), Interval.point(0))
         for i in range(2)]
    contraction = max(sum(a[i][j].maxabs() for j in range(2)) for i in range(2))
    inclusion = all(k[i].strictly_inside(X[i]) for i in range(2))
    certified = inclusion and contraction < 1
    return {
        'center': [str(t) for t in c], 'box': [x.to_json() for x in X],
        'krawczyk_image': [x.to_json() for x in k],
        'contraction_infinity_exact': str(contraction),
        'strict_inclusion': inclusion,
        'certified_unique_root_in_box': certified,
        'spectral_interpretation': 'H(a,b)=[[a,1],[1,b]], det(H)=0, trace(H)=5/2',
        'model_only_not_weil': True,
    }


def numerical_spectral_newton(a0, b0, maxiter=25):
    a,b = float(a0),float(b0)
    for i in range(maxiter):
        f1,f2 = a*b-1, a+b-2.5
        det = b-a
        if abs(det)<1e-13:
            raise ValueError('singular Jacobian along Newton path')
        da = (f1-a*f2)/det
        db = (-f1+b*f2)/det
        a,b = a-da,b-db
        if max(abs(da),abs(db))<1e-14:
            return a,b,i+1,max(abs(a*b-1),abs(a+b-2.5))
    raise ArithmeticError('spectral Newton did not converge')


def zeta_newton_2d(sigma, t, digits=85, maxiter=30):
    """Genuine 2-variable Newton. Not a validated interval zeta evaluator."""
    with mp.workdps(digits):
        z = mp.mpc(mp.mpf(str(sigma)), mp.mpf(str(t)))
        iters = 0
        for iters in range(1, maxiter+1):
            w = mp.zeta(z)
            d = mp.zeta(z, derivative=1)
            u,v = mp.re(w),mp.im(w)
            aa,bb = mp.re(d),mp.im(d)
            det = aa*aa+bb*bb
            if det == 0:
                raise ArithmeticError('zeta derivative vanishes')
            delta_sigma = (aa*u+bb*v)/det
            delta_t = (aa*v-bb*u)/det
            z -= mp.mpc(delta_sigma,delta_t)
            if abs(delta_sigma)+abs(delta_t) < mp.mpf(10)**(-(digits-12)):
                break
        residual = abs(mp.zeta(z))
        return {'start': {'sigma':str(sigma), 't':str(t)},
                'sigma':mp.nstr(mp.re(z), min(digits-8,70)),
                't':mp.nstr(mp.im(z), min(digits-8,70)),
                'residual_abs':mp.nstr(residual,15),
                'iterations':iters,
                'converged_numerically': bool(residual < mp.mpf(10)**(-(digits-14))),
                'certified_by_interval': False}


def v24_checks(v24_json):
    getcontext().prec = 75
    data = json.loads(Path(v24_json).read_text(encoding='utf-8'))
    pairs=[]
    for x in data['test_deltas']:
        delta = Decimal(Q(x['delta']).numerator) / Decimal(Q(x['delta']).denominator)
        bound = Decimal(x['max_row_sum_bound'])
        pairs.append({'delta':str(delta), 'matrix_model_row_sum_bound':str(bound),
                      'quadratic_coefficient':str(bound/(delta*delta)),
                      'under_illustrative_margin_0p001':bool(bound < Decimal('0.001'))})
    coefs=[Decimal(z['quadratic_coefficient']) for z in pairs]
    relative=(max(coefs)-min(coefs))/max(coefs)
    return {'source':Path(v24_json).name,'v24_schema':data['schema'],
            'precise_scope':'Unprojected 10-mode V24 model error; NOT a projected Weil positivity certificate',
            'observations':pairs, 'quadratic_coefficient_relative_spread':str(relative),
            'quadratic_scaling_consistent_to_1e-25':bool(relative < Decimal('1e-25')),
            'illustrative_margin_is_not_from_real_weil_matrix':True}


def perturbation_gate_exact(lower_spectral_margin, upper_error_norm, upper_projection_norm):
    """Purely conditional Weyl/norm gate; all arguments rational upper/lower bounds.

    If mu <= lambda_min(M), ||E||_2 <= eps, and ||B||_2 <= b,
    then lambda_min(B* M B + B* E B) requires a lower bound for
    lambda_min(B* M B), NOT merely mu if B has null directions.
    Here `lower_spectral_margin` MUST already be for the projected matrix.
    """
    mu=Q(lower_spectral_margin)
    eps=Q(upper_error_norm)
    b=Q(upper_projection_norm)
    if eps < 0 or b < 0:
        raise ValueError('norm bounds must be nonnegative')
    allowance=b*b*eps
    return {'projected_margin_lower_bound':str(mu),
            'projected_error_upper_bound':str(allowance),
            'positive_definiteness_preserved_conditional':bool(mu>allowance),
            'requires_verified_same_matrix_and_projection':True}


def write_result(output, v24_json=None):
    spectral=[]
    for start in [(1.8, .7), (.7,1.8)]:
        a,b,it,res=numerical_spectral_newton(*start)
        # Round to exact rational centers; these intervals give exact certificates.
        c=[Q(str(round(a,10))),Q(str(round(b,10)))]
        proof=krawczyk_2d(c,[Q(1,1000),Q(1,1000)])
        spectral.append({'initial_guess':list(start),'newton_root':[a,b],
                         'iterations':it,'approximate_residual':res,'exact_krawczyk':proof})
    zeta=[zeta_newton_2d(.62,14.10),zeta_newton_2d(.56,21.06),zeta_newton_2d(.44,25.05)]
    output_dict={'schema':'NUMERIA.NL1.ES1.2026-11-03.1',
                 'status':'EXACT_FINITE_MODEL_CERTIFICATE_AND_UNCERTIFIED_ZETA_EXPLORATION',
                 'system':'F(a,b)=(ab-1,a+b-5/2)',
                 'spectral_certifications':spectral,
                 'zeta_zeros_exploratory':zeta,
                 'finite_model_implies_RH':False,
                 'zeta_zeros_interval_certified':False,
                 'environment':{'python':platform.python_version(),'mpmath':mp.__version__,
                                'exact_interval_backend':'fractions.Fraction'},
                 'limitations':['The Krawczyk certificates concern an illustrative 2x2 Hermitian pencil, not the Weil form.',
                                'Zeta 2D Newton uses ordinary arbitrary-precision mpmath point evaluations; residual size is not a proof.',
                                'No nontrivial zeta zero off the critical line was found, but the search cannot exclude such zeros.',
                                'V24 projection B*(G+P)B and analytic assumptions have not been checked by this program.']}
    if v24_json:
        output_dict['v24_model_envelope_check']=v24_checks(v24_json)
        # Fictional margin is used ONLY for testing arithmetic of the implication.
        output_dict['hypothetical_projection_gate_demo']=perturbation_gate_exact(
            Q(1,1000),Q('0.000095495337151444193475216932075531708268069'),Q(1))
        output_dict['hypothetical_projection_gate_demo']['actual_weil_margin_verified']=False
    output = Path(output)
    output.write_text(json.dumps(output_dict,ensure_ascii=False,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    return output_dict


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default='RESULTADO_NL1_ES1.json')
    parser.add_argument('--v24',default=None,help='local copy of the official V24 JSON')
    args=parser.parse_args()
    data=write_result(args.out,args.v24)
    print('saved',args.out,'exact_spectral_certificates',
          sum(x['exact_krawczyk']['certified_unique_root_in_box'] for x in data['spectral_certifications']),
          'zeta_point_roots',sum(x['converged_numerically'] for x in data['zeta_zeros_exploratory']))

if __name__=='__main__':
    main()
