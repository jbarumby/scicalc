#!/usr/bin/env python3
"""P73: deterministic cell-by-cell rational implication for Suzuki's zero-side sum.

This does NOT authenticate the input as the first zeros of zeta. An external
FLINT/Arb generation + independent analytic theorem audit remains mandatory.
The algorithm never uses pseudo-randomness, floats, or primality sampling.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'NUMERIA_P72_ZERO_SIDE'))
import verify_spectral_rational as P72

S = P72.SCALING


def spectral_lower(t, ordinates, pi_interval):
    """Outward enclosure of S_N(t), returning its lower integer scaled by 10**75."""
    lower_sum = 0
    pi_lo, pi_hi = pi_interval
    for lo, hi in ordinates:
        # Obtain nearest 2*pi turn count using rational midpoint only.
        midpoint_over_two_pi = F((lo+hi), 2)*t / (pi_lo+pi_hi)
        turns = P72.floor_frac(midpoint_over_two_pi + F(1, 2))
        phase = (P72.floor_frac(F(lo)*t)-2*turns*pi_hi,
                 P72.ceil_frac(F(hi)*t)-2*turns*pi_lo)
        cos_lo, cos_hi = P72.cos_iv(phase)
        lower_sum += 2*(S-cos_hi)*S*S//(hi*hi)
    return F(lower_sum, S)


def certify(zeros, n=20, start=F(18), end=F('46.2'),
            bucket=F('0.05'), max_depth=19, max_cells=100000):
    if not (F(0) < start < end and F(0) < bucket <= F(1)):
        raise ValueError('Invalid positive start, end, or bucket')
    if n<1 or max_depth<1 or max_cells<1:
        raise ValueError('Invalid count/depth/max_cells')
    ordinates = P72.read_zero_intervals(zeros, n)
    pi_interval = P72.pi_iv()
    tail_sum_upper = P72.positive_height_tail_bound()
    zero_eval_cache = {}
    tail_cache = {}

    def point(t):
        if t not in zero_eval_cache:
            zero_eval_cache[t] = spectral_lower(t, ordinates, pi_interval)
        return zero_eval_cache[t]

    def cell_tail(upper_endpoint):
        ceiling_bucket = P72.ceil_frac(upper_endpoint / bucket)
        if ceiling_bucket not in tail_cache:
            majorant_t = ceiling_bucket * bucket
            exp_upper = P72.exp_up_rational(majorant_t/2)
            tail_cache[ceiling_bucket] = (3+exp_upper)*tail_sum_upper
        return tail_cache[ceiling_bucket]

    stack = [(start, end, 0)]
    complete = 0
    deepest = 0
    worst = None
    min_cell = None
    transcript=[]
    while stack:
        a,b,depth = stack.pop()
        h=b-a
        # |S''_N| <= 2*n; lower than secant by at most n*h^2/4.
        curve_loss=F(n,4)*h*h
        lower=min(point(a),point(b))-curve_loss-cell_tail(b)
        if lower>0:
            complete += 1
            deepest=max(deepest, depth)
            if worst is None or lower<worst:
                worst, min_cell=lower, (a,b)
            transcript.append((a,b,lower))
        else:
            if depth>=max_depth:
                raise ArithmeticError(f'Cannot prove cell [{a},{b}] at depth {depth}; lower={lower}')
            if complete+len(stack)+2>max_cells:
                raise ArithmeticError(f'Exceeded max_cells={max_cells}')
            mid=(a+b)/2
            stack.append((mid,b,depth+1))
            stack.append((a,mid,depth+1))

    transcript.sort(key=lambda z:z[0])
    if transcript[0][0]!=start or transcript[-1][1]!=end:
        raise ArithmeticError('Cells do not cover entire interval')
    digest=sha256()
    last=start
    for a,b,lower in transcript:
        if last!=a:raise ArithmeticError('Gap in cell coverage')
        digest.update(f'{a}|{b}|{lower}\n'.encode('ascii'))
        last=b
    if last != end:raise ArithmeticError('Incomplete coverage')
    return {
        'schema':'NUMERIA.P73.DETERMINISTIC_CELL_PROOF.1',
        'start_t':str(start),'end_t':str(end),'number_of_zeros':n,
        'cell_count':complete,'zero_side_evaluations':len(zero_eval_cache),
        'max_depth':deepest,'tail_bucket_width':str(bucket),
        'tail_bucket_evaluations':len(tail_cache),
        'minimum_cell': [str(min_cell[0]),str(min_cell[1])],
        'min_lower_rational':[str(worst.numerator),str(worst.denominator)],
        'approx_min_lower':format(float(worst),'.12g'),
        'cell_coverage_sha256':digest.hexdigest(),
        'strictly_positive_conditional_on_zero_data_and_external_theorems':True,
        'zero_provenance': 'UNAUTHENTICATED INPUT; reproduce zeta zeros with FLINT/Arb before mathematical certification',
        'external_dependencies': ['Suzuki spectral identity (1.3)', 'Platt-Trudgian RH below H=3e12', 'Bellotti-Wong zero counting', 'FLINT/Arb zero isolation'],
        'rigor_note':'Cells use exact Fractions and outward fixed-point cosine; random search never enters the certificate.',
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--zeros',type=Path,required=True)
    p.add_argument('--n',type=int,default=20)
    p.add_argument('--from-t',dest='from_t',type=F,default=F(18))
    p.add_argument('--to',type=F,default=F('46.2'))
    p.add_argument('--tail-bucket',type=F,default=F('0.05'))
    p.add_argument('--max-depth',type=int,default=19)
    p.add_argument('--max-cells',type=int,default=100000)
    p.add_argument('--expected-zeros-sha256',required=True)
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    raw=args.zeros.read_bytes()
    hashed=sha256(raw).hexdigest()
    if hashed!=args.expected_zeros_sha256:
        raise SystemExit(f'Zeros SHA mismatch: got {hashed}, expected {args.expected_zeros_sha256}')
    result=certify(json.loads(raw),args.n,args.from_t,args.to,args.tail_bucket,args.max_depth,args.max_cells)
    result['zeros_file_sha256']=hashed
    out=json.dumps(result,indent=2,ensure_ascii=False)+'\n'
    if args.out:args.out.write_text(out,encoding='utf-8')
    print(out)

if __name__=='__main__':main()
