#!/usr/bin/env python3
"""P115: complete 169x169 refined Arb matrix -> exact fast FLINT certificate.

Consumes all 24 ORIGINAL P111 JSON row shards, checks complete coverage,
identical source hash, per-pair proof consistency, correct profile and (if
provided) 14,365 exact rational overlaps with the P109 original.
Outputs exact integer Frobenius radius and C-level fmpz characteristic
polynomial sign proof. No certificate can be issued from incomplete shards.
"""
import argparse,hashlib,json,math,sys,time
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from flint import fmpz_mat

SCALE=1<<80
N=169
SHARDS=24
PRIME_MAX=1000
def q(m,e): return Q(int(m))*Q(2)**int(e)
def ball(v):
    if not isinstance(v,list) or len(v)!=4:raise ValueError('invalid Arb binary ball')
    m,e,r,f=v;mid=q(m,e);rad=q(r,f)
    if rad<0:raise ValueError('negative interval radius')
    return mid-rad,mid+rad
def encoded(v):
    lo,hi=ball(v)
    mid=(lo+hi)/2;rad=(hi-lo)/2
    im=round(mid*SCALE)
    er=rad+abs(mid-Q(im,SCALE))
    e=(er.numerator*SCALE+er.denominator-1)//er.denominator
    return im,e

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--rows',type=Path,default=Path('.'))
    ap.add_argument('--old-matrix',type=Path,default=None)
    ap.add_argument('--outprefix',default='P115_169')
    a=ap.parse_args()
    paths=sorted(a.rows.glob('P111_1000_ROWS_*.json'))
    if len(paths)!=SHARDS:raise RuntimeError(f'INCOMPLETE: expected {SHARDS} shards; found {len(paths)}')
    centers=[1]+[p for p in range(2,1001) if all(p%d for d in range(2,math.isqrt(p)+1))]
    if len(centers)!=N or centers[-1]!=997:raise ValueError('wrong prime centres')
    old=None
    if a.old_matrix:
        old=json.loads(a.old_matrix.read_text())['entries']
        if len(old)!=14365:raise ValueError('incomplete original P109 comparison matrix')
    seen={};sha_set=set();summaries=[];diag=[]
    overlap_count=0
    # Stream shards to avoid retaining bulky per-panel proofs for all pairs.
    for file in paths:
        d=json.loads(file.read_text())
        if (d.get('status')!='P109_ROW_SHARD_INTERVALS_CONDITIONAL_ON_GL_CAUCHY_ARITHMETIC'
            or d.get('epsilon')!='1/10000000000000000'
            or d.get('bits')!=256 or d.get('max_prime')!=1000
            or d.get('matrix_dimension')!=N
            or d.get('primepower_upper_integer')!=2500
            or d.get('shards')!=SHARDS or d.get('pilot_pair') is not None
            or d.get('centers')!=centers):
            raise ValueError('mismatched P111 settings in '+file.name)
        profile=d.get('test_profile',{})
        if (profile.get('z0')!='127/128'
            or profile.get('inner_GL_error_per_panel')!='1e-22'
            or profile.get('outer_GL_error_per_panel')!='1e-20'):
            raise ValueError('wrong quadrature profile')
        index=d.get('shard_index')
        if index is None or index not in range(SHARDS):
            raise ValueError('missing/out-of-range shard index')
        i=(index*N)//SHARDS; j=((index+1)*N)//SHARDS
        if (d['row_start'],d['row_end'])!=(i,j):
            raise ValueError('wrong integer row partition')
        pairs={f'{centers[k]},{centers[m]}' for k in range(i,j) for m in range(k,N)}
        if set(d['entries'])!=pairs:raise ValueError('missing/extra pair in '+file.name)
        if d.get('diagonal_full_proof_record',{}).get('entry') != d.get('diagonal_reference'):
            raise ValueError('lost diagonal analytic proof')
        diag.append(d['diagonal_reference'])
        sha_set.add(d['source_sha256'])
        for key,value in d['entries'].items():
            if key in seen:raise ValueError('duplicated pair '+key)
            if d['proof_records'].get(key,{}).get('entry')!=value:
                raise ValueError('entry/proof record mismatch '+key)
            ball(value) # reject negative radii
            if old is not None:
                lo,hi=ball(value)
                oldlo,oldhi=ball(old[key])
                if max(lo,oldlo)>min(hi,oldhi):
                    raise ArithmeticError('NEW AND OLD ARB ENCLOSURES DISJOINT: '+key)
                overlap_count+=1
            seen[key]=value
        summaries.append({
            'index':index, 'first_row':i,'last_row':j,'entries':len(d['entries']),
            'source_file_sha256':hashlib.sha256(file.read_bytes()).hexdigest(),
            'outer_panels':d.get('outer_panels_verified'),
            'inner_ellipses':d.get('inner_ellipse_count'),
            'inner_proof_sha256':d.get('inner_ellipse_proof_sha256')
        })
        print('VERIFIED_SHARD',index,len(d['entries']),flush=True)
    if len(sha_set)!=1 or len(seen)!=N*(N+1)//2:
        raise ValueError('inconsistent source hashes or incomplete matrix')
    # The repeated diagonal enclosures must overlap exactly.
    if max(ball(d)[0] for d in diag)>min(ball(d)[1] for d in diag):
        raise ArithmeticError('24 diagonal balls are contradictory')
    M=[[0]*N for _ in range(N)];E=[[0]*N for _ in range(N)]
    for i,x in enumerate(centers):
        for j in range(i,N):
            y=centers[j]
            mid,er=encoded(seen[f'{x},{y}'])
            M[i][j]=M[j][i]=mid
            E[i][j]=E[j][i]=er
    ss=sum(er*er for row in E for er in row)
    rho=math.isqrt(ss)
    rho+=rho*rho<ss
    # The float eigensolver proposes a shift only; the exact sign comes below.
    Mf=np.asarray(M,dtype=np.float64)/SCALE
    hint=float(np.linalg.eigvalsh(Mf)[0])
    if not math.isfinite(hint) or hint<=0:raise RuntimeError('no positive heuristic candidate')
    candidate=Q(repr(hint))*Q(9,10)*SCALE
    shift=candidate.numerator//candidate.denominator
    if shift<=0:raise RuntimeError('nonpositive shift')
    # C-level exact fmpz characteristic-polynomial method, proven for
    # symmetric matrices by Descartes and reality of all eigenvalues.
    started=time.monotonic()
    A=fmpz_mat(N,N,[int(M[i][j]-(shift if i==j else 0))
                   for i in range(N) for j in range(N)])
    co=[int(c) for c in A.charpoly().coeffs()]
    if len(co)!=N+1 or co[-1]!=1:raise ArithmeticError('invalid charpoly')
    transformed=[(-1)**(N+k)*c for k,c in enumerate(co)]
    coeffs_positive=all(c>0 for c in transformed)
    certified=bool(coeffs_positive and shift>rho)
    result={
        'status':'INTERVAL_169_POSITIVE_CONDITIONAL' if certified else 'INCONCLUSIVE_169',
        'centers':centers,'dimension':N,'number_of_entries':len(seen),
        'precision_bits':256,'epsilon':'1/10000000000000000',
        'P111_source_sha256':sha_set.pop(),'shards':summaries,
        'rational_overlaps_with_P109':overlap_count if old is not None else None,
        'source_matrix':'14,365 P111 Arb intervals with outward rounding to 2^-80 integers',
        'rho_Frobenius_upper_rational':str(Q(rho,SCALE)),
        'rho_Frobenius_upper_float':float(Q(rho,SCALE)),
        'numeric_eigenvalue_hint_only':hint,
        'exact_integer_shift':shift,
        'shift_rational':str(Q(shift,SCALE)),
        'shift_float':float(Q(shift,SCALE)),
        'charpoly_all_positive':coeffs_positive,
        'positive_coefficients':sum(c>0 for c in transformed),
        'required_positive_coefficients':N+1,
        'verified_min_eigenvalue_lower_rational':str(Q(shift-rho,SCALE)) if certified else None,
        'verified_min_eigenvalue_lower_float':float(Q(shift-rho,SCALE)) if certified else None,
        'charpoly_seconds':time.monotonic()-started,
        'charpoly_sha256':hashlib.sha256(("\n".join(map(str,transformed))+"\n").encode()).hexdigest(),
        'scope':'finite Weil test subspace; depends on analytic GL/Cauchy and FLINT/Arb correctness; not RH'
    }
    prefix=Path(a.outprefix)
    prefix.with_name(prefix.name+'_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    prefix.with_name(prefix.name+'_CHARPOLY_COEFFS.txt').write_text(
       "\n".join(map(str,transformed))+"\n")
    prefix.with_name(prefix.name+'_COMPACT_MATRIX.json').write_text(
       json.dumps({'centers':centers,'entries':seen,
           'P111_source_sha256':result['P111_source_sha256']})+'\n')
    print('P115_FULL_169',json.dumps({
        'status':result['status'],'rho':result['rho_Frobenius_upper_float'],
        'lambda_hint':hint,'shift':result['shift_float'],
        'coeffs_positive':result['positive_coefficients'],
        'lower_bound':result['verified_min_eigenvalue_lower_float'],
        'elapsed_charpoly':result['charpoly_seconds']
    }),flush=True)
    if not certified:raise RuntimeError('INCONCLUSIVE: exact 169x169 interval matrix positivity not certified')
if __name__=='__main__':
    sys.exit(main())
