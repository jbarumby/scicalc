#!/usr/bin/env python3
"""Independent exact audit of existing P109 fixed 500 and 1000 matrices.

All signs use exact integer Bareiss + an outward dyadic Frobenius radius.
The original numeric estimate proposes shifts only; never certifies signs.
"""
import argparse, json, math, hashlib
from fractions import Fraction as Q
from pathlib import Path

SCALE=1<<88
def dyadic(m,e):return Q(int(m))*Q(2)**int(e)
def ball(v):
    if len(v)!=4:raise ValueError('malformed dyadic Arb ball')
    m,e,r,f=v
    mid=dyadic(m,e);rad=dyadic(r,f)
    if rad<0:raise ValueError('negative radius')
    mi=round(mid*SCALE)
    err=rad+abs(mid-Q(mi,SCALE))
    err_up=(err.numerator*SCALE+err.denominator-1)//err.denominator
    return mi,err_up

def bareiss(A,shift):
    n=len(A)
    x=[[A[i][j]-(shift if i==j else 0) for j in range(n)]
       for i in range(n)]
    prev=1
    for k in range(n):
        pivot=x[k][k]
        if pivot<=0:return False,k
        row=x[k]
        for i in range(k+1,n):
            other=x[i];b=other[k]
            for j in range(i,n):
                v,r=divmod(pivot*other[j]-b*row[j],prev)
                if r:raise ArithmeticError('nonexact Bareiss division')
                other[j]=x[j][i]=v
            other[k]=0
        prev=pivot
    return True,n

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--bound',type=int,required=True,choices=(500,1000))
    a=p.parse_args()
    base=Path('.')
    compact=base/f'P109_{a.bound}_COMPACT_MATRIX.json'
    auditpath=base/f'P109_{a.bound}_AUDIT.json'
    data=json.loads(compact.read_text())
    old=json.loads(auditpath.read_text())
    centers=data['centers'];n=len(centers)
    if n != {500:96,1000:169}[a.bound]:raise ValueError('incorrect centre count')
    if data['max_prime']!=a.bound or data['epsilon']!='1/100000000000' or data['precision_bits']!=256:
        raise ValueError('incompatible experimental inputs')
    if len(data['entries'])!=n*(n+1)//2:raise ValueError('missing entries')
    M=[[0]*n for _ in range(n)]
    E=[[0]*n for _ in range(n)]
    for i,x in enumerate(centers):
        for j in range(i,n):
            k=f'{x},{centers[j]}'
            m,e=ball(data['entries'][k])
            M[i][j]=M[j][i]=m
            E[i][j]=E[j][i]=e
    errorsq=sum(x*x for row in E for x in row)
    rad=math.isqrt(errorsq)
    rad+=rad*rad<errorsq
    rho=Q(rad,SCALE)
    if str(rho)!=str(Q(old['full_error_radius_upper_rational'])):
        raise ValueError('Frobenius error does not match original exact audit')
    h=old['full_approx_lambda_min_not_certified']
    records=[]
    success=False
    for factor in (0.99999,0.9999,0.999,0.995,0.99):
        sigma=math.floor(h*factor*SCALE)
        if sigma<=rad:
            records.append({'factor':factor,'candidate_above_radius':False})
            continue
        is_pos,pivots=bareiss(M,sigma)
        result={'factor':factor,'candidate_above_radius':True,
                'integer_pivots_positive':pivots,'shift_passed':is_pos,
                'shift_over_error_bound_rational':str(Q(sigma-rad,SCALE))}
        records.append(result)
        if is_pos:
            success=True
            break
    output={
       'status':'FINITE_MATRIX_POSITIVE_CONDITIONAL' if success else 'INCONCLUSIVE_MATRIX_INTERVAL',
       'max_prime':a.bound,'matrix_dimension':n,'number_of_entries':n*(n+1)//2,
       'source_artifact_sha256':hashlib.sha256(compact.read_bytes()).hexdigest(),
       'original_audit_sha256':hashlib.sha256(auditpath.read_bytes()).hexdigest(),
       'lambda_numerical_estimate_is_only_heuristic':h,
       'Frobenius_error_radius_upper_rational':str(rho),
       'integer_scale':str(SCALE),'tests':records,
       'lambda_min_lower_rational':records[-1]['shift_over_error_bound_rational'] if success else None,
       'lambda_min_lower_display':float(Q(records[-1]['shift_over_error_bound_rational'])) if success else None,
       'scope':'Finite 96x96 or 169x169 Weil subspace only. Conditional on the Arb/Cauchy GL inclusion proofs; no RH claim.'
    }
    Path(f'P109_{a.bound}_REVISED_EXACT_AUDIT.json').write_text(json.dumps(output,indent=2)+'\n')
    print('REVISED_AUDIT',a.bound,output['status'],
          'lambda_lower',output['lambda_min_lower_display'],
          'rho',float(rho),
          'verified_integer_pivots',records[-1].get('integer_pivots_positive',0),flush=True)
    # The 500 test is expected to be certified; never silently regress to an
    # inconclusive green job after an apparent coding success.
    if a.bound==500 and not success:raise RuntimeError('expected 96x96 rational positivity was not reproduced')

if __name__=='__main__':main()
