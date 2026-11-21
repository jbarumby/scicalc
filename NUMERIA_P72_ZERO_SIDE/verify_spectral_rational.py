#!/usr/bin/env python3
"""NUMERIA P72: fixed-point rational checks for Suzuki Psi(t) on a finite interval.

This program does NOT validate the supplied ordinates as zeta zeros. For a fully
rigorous certificate, feed it the exact dyadic intervals exported by
export_zeros_flint.py (python-flint/Arb) and verify the imported theorems.
With preview_zeros.json the computation is a rigorous implication CONDITIONAL
on the listed preview intervals actually containing zeta zero ordinates.

No primes enumerated. No floating-point computations used in inequalities.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from math import factorial
import json, hashlib
from pathlib import Path

SCALING = 10**75
H = 3*10**12
POLY_ORDER = 40
REM_DEN = 10**50
if not 4**(2*POLY_ORDER+2)*REM_DEN < factorial(2*POLY_ORDER+2):
    raise ArithmeticError("Taylor remainder proof failed")


def floor_frac(v): return v.numerator//v.denominator

def ceil_frac(v): return -((-v.numerator)//v.denominator)

def scale_range(lo, hi): return (floor_frac(lo*SCALING),ceil_frac(hi*SCALING))

def plus(a,b): return (a[0]+b[0], a[1]+b[1])

def times(a,b):
    values=[x*y for x in a for y in b]
    return (min(values)//SCALING, -((-max(values))//SCALING))


def atan_range_inv(q,n):
    s=sum((Q((-1)**k,(2*k+1)*q**(2*k+1)) for k in range(n)),Q(0))
    r=Q((-1)**n,(2*n+1)*q**(2*n+1))
    return s+min(Q(0),r),s+max(Q(0),r)

@lru_cache(None)
def pi_iv():
    a=atan_range_inv(5,100); b=atan_range_inv(239,100)
    out=scale_range(16*a[0]-4*b[1],16*a[1]-4*b[0])
    if not 3*SCALING < out[0] < out[1] < 4*SCALING:
        raise ArithmeticError('invalid pi interval')
    return out

@lru_cache(None)
def coeffs():
    return [scale_range(Q((-1)**k,factorial(2*k)),Q((-1)**k,factorial(2*k))) for k in range(POLY_ORDER+1)]


def cos_iv(phase):
    if phase[0] < -4*SCALING or phase[1] >4*SCALING:
        raise ArithmeticError('bad angle reduction')
    z=times(phase,phase)
    acc=coeffs()[-1]
    for c in coeffs()[-2::-1]: acc=plus(c,times(z,acc))
    remainder=(SCALING+REM_DEN-1)//REM_DEN
    return max(-SCALING,acc[0]-remainder),min(SCALING,acc[1]+remainder)


def dyad(v):
    m,e=v; m=int(m); e=int(e)
    return Q(m * (1 << e), 1) if e>=0 else Q(m,1<<(-e))


def read_zero_intervals(data,n):
    rows=data['zeros'][:n]
    if len(rows)!=n: raise ValueError('insufficient zeros')
    result=[]
    previous=Q(0)
    for j,row in enumerate(rows,1):
        if int(row['index'])!=j: raise ValueError('nonsequential index')
        if 'lo_dyad' in row:
            lo=dyad(row['lo_dyad']); hi=dyad(row['hi_dyad'])
        else:
            lo=Q(row['lo_decimal']); hi=Q(row['hi_decimal'])
        if not previous < lo <= hi: raise ArithmeticError('unordered or overlapping zero intervals')
        result.append(scale_range(lo,hi))
        previous=hi
    return result


def exp_up_rational(x):
    # exact series plus a rigorous geometric majorant for positive rational x
    # for tail first omitted term after k=M; ratio thereafter <= x/(M+2).
    M=200
    partial=Q(0)
    term=Q(1)
    for k in range(M+1):
        if k: term=term*x/k
        partial+=term
    first=term*x/(M+1)
    if not x < M+2: raise ArithmeticError('invalid exp tail')
    bound=partial+first/(1-x/Q(M+2))
    # Directed round upwards to 1e-9 for compact external certificate.
    den=10**9
    upper=Q(ceil_frac(bound*den),den)
    if upper<=bound: raise ArithmeticError('exp majorant rounding failed')
    return upper


def positive_height_tail_bound():
    """Bellotti-Wong zero counting, with N(H) subtraction retained.

    For N(y)=M(y)+R(y), |R(y)|<= a log y + b loglog y + c.
    Using loglog y<=log y (y>=e), log H<30, 3<pi<22/7,
    log(2pi)>1 and a+b=.34536, c=8.08344:

    sum_{Im rho > H}1/(Im rho)^2 = 2 int_H^inf N(u)/u^3 du - N(H)/H^2
      <= (log(H/(2pi))+1)/(2pi H)
       + ((a+b)*(2log H+1/2)+2c)/H^2
      < 5/H + 38/H^2.
    """
    exp30_lo=sum((Q(30**j,factorial(j)) for j in range(90)),Q(0))
    if exp30_lo <= H: raise ArithmeticError('log H<30 failed')
    a_b=Q(10076+24460,100000)
    c=Q(808344,100000)
    error_bound=a_b*(2*30+Q(1,2))+2*c
    if error_bound>=38: raise ArithmeticError('Bellotti-Wong remainder >=38')
    # log(2pi)>1 since pi>3 and e<3; 2pi>6>e.
    # 1/(2pi)<1/6 since pi>3; hence main <30/(6H)=5/H.
    # Strict inequality comes from the explicit results above.
    return Q(5,H)+Q(38,H*H)


def certify(data,n,tmax_milli,step_milli=5,begin_milli=18000):
    if not 0<step_milli<=10: raise ValueError('step must be between 1 and 10 milliseconds')
    if not begin_milli<tmax_milli or (tmax_milli-begin_milli)%step_milli:
        raise ValueError('interval endpoints must align with grid')
    G=read_zero_intervals(data,n)
    pi=pi_iv()
    minimum=None; min_at=None; grid_sha=hashlib.sha256()
    for t in range(begin_milli,tmax_milli+1,step_milli):
        subtotal=0
        for g in G:
            # nearest multiple of 2*pi; the resulting phase inclusion is checked
            turns=((g[0]+g[1])*t +1000*(pi[0]+pi[1]))//(2000*(pi[0]+pi[1]))
            if turns<0: raise ArithmeticError('negative turns')
            phase=(g[0]*t//1000 - 2*turns*pi[1],
                   -((-g[1]*t)//1000) - 2*turns*pi[0])
            c_lo,c_hi=cos_iv(phase)
            term_lo=2*(SCALING-c_hi)*SCALING*SCALING//(g[1]*g[1])
            subtotal+=term_lo
        if minimum is None or subtotal < minimum:
            minimum=subtotal; min_at=t
        grid_sha.update(f'{t},{subtotal}\n'.encode('ascii'))

    # For S_N''(t)=2*sum cos(gamma_j*t)<=2N, piecewise linear interpolation
    # from endpoints loses at most N*h^2/4 (strictly better than Lipschitz).
    loss=Q(n*step_milli*step_milli,4*1000*1000)
    exp_cap=exp_up_rational(Q(tmax_milli,2000))
    tail=(3+exp_cap)*positive_height_tail_bound()
    grid=Q(minimum,SCALING)
    uniform=grid-loss-tail
    validated=bool(data.get('zero_certification')=='arb_enclosures')
    record={
        'schema':'NUMERIA.P72.ZERO_SIDE.1',
        'start_t':str(Q(begin_milli,1000)),
        'end_t':str(Q(tmax_milli,1000)),
        'step_t':str(Q(step_milli,1000)),
        'zero_count':n,
        'zero_source':data.get('zero_certification','unknown'),
        'certified_zeros_from_flint':validated,
        'number_of_grid_points':1+(tmax_milli-begin_milli)//step_milli,
        'minimum_grid_at_t':str(Q(min_at,1000)),
        'grid_min_scaled':str(minimum),
        'grid_sha256':grid_sha.hexdigest(),
        'interpolation_loss_rational':[str(loss.numerator),str(loss.denominator)],
        'exp_end_upper_rational':[str(exp_cap.numerator),str(exp_cap.denominator)],
        'zero_tail_upper_rational':[str(tail.numerator),str(tail.denominator)],
        'uniform_lower_rational':[str(uniform.numerator),str(uniform.denominator)],
        'numerical_summary':{
            'grid_min_lower':f'{float(grid):.12g}',
            'curvature_loss':f'{float(loss):.12g}',
            'high_zero_tail_upper':f'{float(tail):.12g}',
            'uniform_lower_bound':f'{float(uniform):.12g}',
        },
        'positive_given_input_intervals_are_valid':uniform>0,
        'certified_positive':validated and uniform>0,
        'dependency':'Suzuki (1.3); Platt-Trudgian RH for |Im rho|<=3e12; Bellotti-Wong N(T) 2025',
    }
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--zeros',type=Path,default=Path(__file__).with_name('preview_zeros.json'))
    p.add_argument('--n',type=int,default=20)
    p.add_argument('--to',type=str,default='40')
    p.add_argument('--step-ms',type=int,default=5)
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    end=Q(args.to)*1000
    if end.denominator!=1: raise ValueError('end time must have integer milliseconds')
    raw=args.zeros.read_bytes()
    data=json.loads(raw.decode('utf-8'))
    rec=certify(data,args.n,int(end),args.step_ms)
    rec['zeros_file_sha256']=hashlib.sha256(raw).hexdigest()
    print(json.dumps(rec,indent=2))
    if args.out:
        args.out.write_text(json.dumps(rec,indent=2)+'\n',encoding='utf-8')
    # Do not fail when the proof is open: the report will mark status honestly.

if __name__=='__main__': main()
