"""NUMERIA P92-D: analytic rank-two Delta bound for anchor + uniform cluster.

Exact rational regressions are not the proof for ALL real L/h and M:
the proof is elementary, complete, and written in README.md.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json

SCHEMA = 'NUMERIA.P92D.DELTA.STAR.ANALYTIC.1'


def star_norm_square(s, m):
    """Exact squared norm for scale ratio lambda=s^2, rational s>=1, m>=2."""
    s = Q(s)
    if s < 1 or not isinstance(m, int) or m < 2:
        raise ValueError('s>=1 and m>=2 required')
    la = s*s
    return (la*(s-1)**4/4)*sum((1/(la+k)**4 for k in range(m)), Q(0))


def square_sum_upper(s):
    """Decreasing-sum integral bound: sum (lambda+k)^-4 <=lambda^-4+1/(3lambda^3)."""
    s = Q(s)
    if s < 1:
        raise ValueError('s>=1')
    la = s*s
    return 1/la**4 + 1/(3*la**3)


def bound_via_integral(s):
    """Exact rational envelope (1-u)^4(1+3u^2)/12 <=1/12, u=1/s."""
    s=Q(s)
    if s<1:
        raise ValueError('s>=1')
    u=1/s
    return (1-u)**4*(1+3*u*u)/12


def star_weights(s,m):
    """Nonzero Delta[0,j] entries, all rational if s is rational."""
    s=Q(s)
    if s<1 or not isinstance(m,int) or m<2:
        raise ValueError('s>=1 and m>=2')
    return [s*(s-1)**2/(2*(s*s+k)**2) for k in range(m)]


def exact_matrix_sample(s,m):
    """Check the rank-two star eigenvalue squared by exact Fraction."""
    weights=star_weights(s,m)
    squared=sum((v*v for v in weights),Q(0))
    if squared!=star_norm_square(s,m):
        raise ArithmeticError('rank-two identity mismatch')
    if any(v<0 for v in weights):
        raise ArithmeticError('negative Delta entry')
    # Matrix [[0,w^T],[w,0]] acts on e0 by w and on w by ||w||^2 e0.
    # Hence its two nonzero eigenvalues are +/-||w||.
    return squared


def certify():
    rows=[]
    for s in (Q(1),Q(6,5),Q(3,2),Q(2),Q(3),Q(5),Q(10),Q(20)):
        for m in (2,3,5,10,30,60,120):
            exact=exact_matrix_sample(s,m)
            series=sum((1/(s*s+k)**4 for k in range(m)),Q(0))
            envelope=bound_via_integral(s)
            if not (series<=square_sum_upper(s) and exact<=envelope<=Q(1,12)):
                raise ArithmeticError(f'bad exact regression s={s},m={m}')
            rows.append({'sqrt_scale_ratio':str(s),'cluster_count':m,
                         'squared_norm_exact':str(exact),
                         'integral_ceiling_exact':str(envelope),
                         'margin_to_1_over_12_exact':str(Q(1,12)-exact)})
    return {'schema':SCHEMA,
            'theorem':'For ALL real L>=h>0 and all integers M>=2, ||Delta||_2 <= 1/sqrt(12)',
            'proof':'rank-two exact eigenvalues; decreasing-series integral; (1-u)^4(1+3u^2)<=1',
            'bound_sharp_asymptotically':True,
            'star_family_proved_analytically':True,
            'global_Delta_bound_proved':False,
            'full_Preissmann_Lemma6_proved':False,
            'Riemann_Hypothesis_proved':False,
            'rational_regression_count':len(rows),
            'rational_regressions':rows}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    d=certify()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('P92D PASS: 56 Fraction exact regressions; universal STAR-family bound proven analytically')
    print('General Delta bound, full Preissmann and RH remain OPEN')


if __name__=='__main__':
    main()
