#!/usr/bin/env python3
"""ES8-P111: exact dyadic 63/96/169-centre positivity + error-radius audit.

np.linalg.eigvalsh is used only to propose a shift and for trends.
CERTIFICATION uses 100% integer Bareiss pivots and Frobenius bounds.
If eigenvalue is smaller than enclosures: INCONCLUSIVE, NOT negative.
"""
import argparse, csv, glob, hashlib, json, math
from fractions import Fraction as Q
from pathlib import Path

SCALE=1<<88
BOUND_COUNTS={300:63,500:96,1000:169}

def get_centers(bound):
    return tuple([1]+[p for p in range(2,bound+1)
                      if all(p%d for d in range(2,math.isqrt(p)+1))])

def dyadic(m,e):
    return Q(int(m))*Q(2)**int(e)

def ceilq(q):
    return -((-q.numerator)//q.denominator) if q.denominator>0 else None

def outward(v):
    if not isinstance(v,list) or len(v)!=4:
        raise ValueError('invalid 4-word Arb dyadic interval encoding')
    m,e,r,f=v
    mid=dyadic(m,e);rad=dyadic(r,f)
    if rad<0:raise ValueError('negative Arb radius')
    scaled=round(mid*SCALE)
    d=abs(mid-Q(scaled,SCALE))+rad
    error=(d.numerator*SCALE+d.denominator-1)//d.denominator
    return scaled,error

def bareiss_pd(data,shift=0):
    """Fraction-free Sylvester. Return (boolean, positive pivot count)."""
    n=len(data);A=[[int(data[i][j])-(shift if i==j else 0)
                    for j in range(n)] for i in range(n)]
    prev=1
    for k in range(n):
        p=A[k][k]
        if p<=0:return False,k
        rowk=A[k]
        for i in range(k+1,n):
            ri=A[i];aik=ri[k]
            for j in range(i,n):
                numerator=p*ri[j]-aik*rowk[j]
                q,rem=divmod(numerator,prev)
                if rem:raise ArithmeticError('Bareiss division nonexact')
                ri[j]=A[j][i]=q
            ri[k]=0
        prev=p
    return True,n

def frob_error(E,n):
    ss=sum(E[i][j]*E[i][j] for i in range(n) for j in range(n))
    r=math.isqrt(ss)
    return r+(r*r<ss)

def read_shards(paths,bound,expected_shards):
    centers=get_centers(bound);n=len(centers)
    if n!=BOUND_COUNTS[bound]:raise ValueError('unexpected prime cardinality')
    entries={};rows=[];source=set();countinner=countouter=0
    meta=[];diagonals=[]
    if len(paths)!=expected_shards:
        raise ValueError(f'expected {expected_shards} shard artifacts, got {len(paths)}')
    for path in paths:
        obj=json.loads(Path(path).read_text())
        if obj.get('status')!='P109_ROW_SHARD_INTERVALS_CONDITIONAL_ON_GL_CAUCHY_ARITHMETIC':
            raise ValueError('uncompleted shard')
        if tuple(obj['centers'])!=centers or obj['matrix_dimension']!=n or obj['max_prime']!=bound:
            raise ValueError('wrong prime family/matrix dimension')
        if obj['bits']!=256 or obj['epsilon']!='1/10000000000000000' or obj['primepower_upper_integer']!=2500:
            raise ValueError('bad precision/epsilon/prime power bound')
        if obj['pilot_pair'] is not None or obj['shards']!=expected_shards:
            raise ValueError('pilot or incompatible shard')
        i,j=obj['row_start'],obj['row_end']
        if not (0<=i<j<=n):raise ValueError('invalid row range')
        rows+=list(range(i,j))
        exp={f'{centers[k]},{centers[l]}' for k in range(i,j) for l in range(k,n)}
        if set(obj['entries'])!=exp:raise ValueError('missing/extra row entries')
        for key,val in obj['entries'].items():
            if key in entries:raise ValueError('duplicate entry '+key)
            outward(val)
            check=obj['proof_records'].get(key)
            if check is None or check.get('entry')!=val:
                raise ValueError('interval does not match recorded proof '+key)
            entries[key]=val
        source.add(obj['source_sha256'])
        diagonals.append(outward(obj['diagonal_reference']))
        countinner+=obj['inner_ellipse_count']
        countouter+=obj['outer_panels_verified']
        meta.append({'filename':Path(path).name,'row_range':[i,j],
            'sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'inner_ellipse_proof_sha256':obj['inner_ellipse_proof_sha256'],
            'inner_ellipse_count':obj['inner_ellipse_count'],
            'outer_panels_verified':obj['outer_panels_verified']})
    if len(source)!=1 or sorted(rows)!=list(range(n)) or len(entries)!=n*(n+1)//2:
        raise ValueError('missing/overlapping rows, different source or missing entries')
    # Repeated diagonal enclosures in every shard should share a common exact point.
    low=max(Q(m-e,SCALE) for m,e in diagonals)
    high=min(Q(m+e,SCALE) for m,e in diagonals)
    if low>high:raise ValueError('diagonal references contradict each other')
    M=[[0]*n for _ in range(n)]
    E=[[0]*n for _ in range(n)]
    for i,a in enumerate(centers):
        for j in range(i,n):
            b=centers[j];v=entries[f'{a},{b}'];mid,err=outward(v)
            M[i][j]=M[j][i]=mid;E[i][j]=E[j][i]=err
    return centers,M,E,entries,meta,source.pop(),countinner,countouter

def audit(bound,shards):
    import numpy as np
    paths=sorted(glob.glob(f'P111_{bound}_ROWS_*.json'))
    centers,M,E,entries,meta,sha,inner,outer=read_shards(paths,bound,shards)
    n=len(centers)
    probes=sorted(set(x for x in (2,8,16,26,47,63,80,96,110,130,150,169,n) if x<=n))
    rows=[]
    for k in probes:
        radius=frob_error(E,k)
        heuristic=float(np.linalg.eigvalsh(np.asarray([row[:k] for row in M[:k]],
            dtype='float64')/SCALE)[0])
        if not math.isfinite(heuristic):
            raise ValueError('non-finite eigenvalue heuristic')
        # Shift selection is heuristic. Every proof below is exact integer.
        proposed=math.floor(max(0,heuristic)*0.95*SCALE)
        trial=proposed;certified=False;attempts=0;piv=0
        while trial>radius and attempts<16:
            attempts+=1
            good,piv=bareiss_pd([row[:k] for row in M[:k]],trial)
            if good:
                certified=True;break
            trial//=2
        # For inconclusive cases, checking only the MIDPOINT is not
        # sufficient to certify the unknown interval matrix!
        lower=(trial-radius) if certified and trial>radius else None
        rows.append({
            'n':k,'largest_center':centers[k-1],
            'lambda_min_numeric_heuristic':heuristic,
            'error_radius_upper_rational':str(Q(radius,SCALE)),
            'error_radius_upper_for_plot':float(Q(radius,SCALE)),
            'lambda_lower_rational':str(Q(lower,SCALE)) if lower is not None else None,
            'lambda_lower_for_plot':float(Q(lower,SCALE)) if lower is not None else None,
            'verified_positive_conditional':lower is not None and lower>0,
            'exact_bareiss_pivots':piv,'attempts':attempts
        })
        print('TARGET',bound,'PREFIX',k,'TO',centers[k-1],
              'LAMBDA_LOWER',rows[-1]['lambda_lower_for_plot'],
              'ERROR_RADIUS',rows[-1]['error_radius_upper_for_plot'],
              'APPROX',heuristic,'CERTIFIED',rows[-1]['verified_positive_conditional'],flush=True)
    certs=[r for r in rows if r['verified_positive_conditional']]
    status=('POSITIVE_CONDITIONAL' if rows[-1]['verified_positive_conditional']
            else 'INCONCLUSIVE_MATRIX_INTERVAL')
    result={'status':status,'max_prime':bound,'centers':centers,
       'n':n,'epsilon':'1/10000000000000000','bits':256,
       'primepower_upper_integer':2500,'distinct_entries':len(entries),
       'source_sha256':sha,'shards':meta,
       'total_inner_ellipses_verified_during_arb':inner,
       'total_outer_panels_verified_during_arb':outer,
       'tested_prefixes':rows,
       'full_lambda_lower_rational':rows[-1]['lambda_lower_rational'],
       'full_error_radius_upper_rational':rows[-1]['error_radius_upper_rational'],
       'full_approx_lambda_min_not_certified':rows[-1]['lambda_min_numeric_heuristic'],
       'trust_model':'Cauchy-GL ellipse proof + FLINT/Arb; only exact integer Bareiss confirms any certified sign',
       'not_RH_proof':True}
    Path(f'P111_{bound}_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n')
    compact={'max_prime':bound,'centers':centers,'epsilon':'1/10000000000000000',
       'precision_bits':256,'entry_encoding':'[mantissa,exponent,rad_mantissa,rad_exponent]',
       'entries':entries,'audit_file':f'P111_{bound}_AUDIT.json'}
    Path(f'P111_{bound}_COMPACT_MATRIX.json').write_text(json.dumps(compact,indent=2)+'\n')
    with open(f'P111_{bound}_TREND.csv','w',newline='') as fp:
        writer=csv.DictWriter(fp,fieldnames=['max_prime','n','largest_center',
            'lambda_min_numeric_heuristic','lambda_lower_for_plot',
            'error_radius_upper_for_plot','verified_positive_conditional'])
        writer.writeheader()
        for row in rows:
            writer.writerow({'max_prime':bound,**{key:row[key] for key in writer.fieldnames if key!='max_prime'}})
    Path(f'P111_{bound}_REPORT.md').write_text(
        f'# ES8-P111 — Primes through {bound}\n\n'
        f'Centres: {n}; distinct Arb entries: {len(entries)}; epsilon=1e-16; 256 bits.\n\n'
        f'**Status:** {status}\n\n'
        f'Full lambda lower: {rows[-1]["lambda_lower_rational"]}; '
        f'error radius: {rows[-1]["error_radius_upper_rational"]}.\n\n'
        f'Numerical lambda estimate (heuristic, not proof): {rows[-1]["lambda_min_numeric_heuristic"]}.\n\n'
        'The certified lower bound requires the exact integer Bareiss positive '
        'pivots and the rigorous Frobenius error upper bound. '
        'If the bound is inconclusive, this does not imply an eigenvalue <=0. '
        'Finite-dimensional results do not prove RH.\n')
    print('AUDIT_DONE',bound,status,'CENTRES',n,'ENTRIES',len(entries),flush=True)
    if status != 'POSITIVE_CONDITIONAL':
        raise RuntimeError('P111 169x169 remains inconclusive: no false positive certificate')

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-prime',required=True,type=int,choices=[300,500,1000])
    ap.add_argument('--shards',type=int,default=24)
    a=ap.parse_args()
    audit(a.max_prime,a.shards)
