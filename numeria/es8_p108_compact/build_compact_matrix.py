#!/usr/bin/env python3
"""Build small, independently verifiable P108 47x47 interval matrix from CI shards."""
import hashlib, json, glob, csv
from fractions import Fraction as Q
from pathlib import Path
from math import isqrt

FILES=sorted(glob.glob("P108_ROWS_*.json"))
EXPECTED=[1]+[p for p in range(2,200) if all(p%q for q in range(2,isqrt(p)+1))]
SCALE=1<<80

def exact(m,e): return Q(int(m))*Q(2)**int(e)
def mat_interval(v):
    m,e,r,f=v
    x=exact(m,e); y=exact(r,f)
    if y<0:raise ValueError("negative error")
    return x-y,x+y

if len(FILES)!=10:raise RuntimeError(f"expected ten shards, got {len(FILES)}")
entries={};sources=set();rows=[];piece_sha={};refs=[]
sumouter=0;suminner=0
for fn in FILES:
    p=Path(fn);b=p.read_bytes()
    piece_sha[p.name]=hashlib.sha256(b).hexdigest()
    doc=json.loads(b)
    if doc.get('epsilon')!='1/100000000000' or doc.get('bits')!=256:
        raise ValueError('shard experimental settings mismatch')
    if doc.get('status')!='P108_ARBITRARY_MATRIX_ENTRY_INTERVALS_ANALYTICALLY_CONDITIONAL':
        raise ValueError('nonfinite/incomplete shard')
    if doc.get('primepower_upper_integer')!=500 or doc.get('matrix_dimension')!=47 or doc.get('centers')!=EXPECTED:
        raise ValueError('source family mismatch')
    sources.add(doc['source_sha256'])
    i,j=doc['row_start'],doc['row_end']
    if not 0<=i<j<=47:raise ValueError('bad shard indices')
    rows+=list(range(i,j))
    expected_keys={f"{EXPECTED[ii]},{EXPECTED[jj]}" for ii in range(i,j) for jj in range(ii,47)}
    if set(doc['entries'])!=expected_keys:raise ValueError('missing/extra entries in shard')
    for k,v in doc['entries'].items():
        if k in entries:raise ValueError('duplicated entry: '+k)
        mat_interval(v)
        entries[k]=v
    sumouter+=doc['stats'].get('outer_complex_verified',0)
    suminner+=doc['stats'].get('inner_complex_verified',0)
    refs.append({'rows':[i,j], 'filename':p.name,'sha256':piece_sha[p.name]})
if len(sources)!=1 or sorted(rows)!=list(range(47)) or len(entries)!=1128:
    raise ValueError('incompatible source hash, row coverage or pair count')

report=json.load(open('P108_47X47_EXACT_AUDIT.json'))
if report['status']!='ALL_47_PREFIXES_POSITIVE_CONDITIONAL' or len(report['rows'])!=47:
    raise ValueError('P108 matrix audit not complete/positive')
L=report['rows'][-1]['lambda_lower_dyadic_numerator']
Den=report['rows'][-1]['lambda_lower_dyadic_denominator']
if L<=0 or Den!=SCALE or Q(L,Den)<=0:
    raise ValueError('missing strictly positive exact integer lambda certificate')

# Independently rebuild scaled mids and radii, then re-check the final Bareiss
# pivot sequence without trusting prior floating-point eigenvalues.
M=[[0]*47 for _ in range(47)]
R=[[0]*47 for _ in range(47)]
for i,a in enumerate(EXPECTED):
    for j in range(i,47):
        b=EXPECTED[j]
        v=entries[f'{a},{b}']; m,e,r,f=v
        mid=exact(m,e);rad=exact(r,f)
        im=round(mid*SCALE)
        err=rad+abs(mid-Q(im,SCALE))
        er=-(-err.numerator*SCALE//err.denominator)
        M[i][j]=M[j][i]=im
        R[i][j]=R[j][i]=er
ss=sum(u*u for row in R for u in row)
rho=isqrt(ss)+(isqrt(ss)**2<ss)
if rho!=report['rows'][-1]['Frobenius_error_upper_int']:
    raise ValueError('independent Frobenius mismatch')
shift=L+rho
if shift!=report['rows'][-1]['candidate_shift_int']:
    raise ValueError('independent shift mismatch')
A=[[M[i][j]-(shift if i==j else 0) for j in range(47)] for i in range(47)]
prev=1;pivots=0
for k in range(47):
    piv=A[k][k]
    if piv<=0:raise ValueError(f'nonpositive integer Bareiss pivot {k}')
    for i in range(k+1,47):
        aik=A[i][k]
        for j in range(i,47):
            num=piv*A[i][j]-aik*A[k][j]
            z,r=divmod(num,prev)
            if r:raise ValueError('nonexact fraction-free division')
            A[i][j]=A[j][i]=z
    for i in range(k+1,47):A[i][k]=A[k][i]=0
    prev=piv;pivots+=1

compact={
  'status':'47X47_COMPLETE_INTERVAL_MATRIX_AND_EXACT_POSITIVITY_CONDITIONAL',
  'centers':EXPECTED,'epsilon':'1/100000000000','bits':256,
  'primepower_upper_integer':500,'matrix_order':47,'distinct_entries':len(entries),
  'entry_encoding':'[mantissa, exponent, radius_mantissa, radius_exponent] each x=mantissa*2^exponent +/- radius_mantissa*2^radius_exponent; Arb binary balls',
  'entries':entries,'source_sha256':sources.pop(),
  'row_shard_provenance':refs,'outer_panels_checked':sumouter,
  'inner_ellipses_checked':suminner,
  'lambda_lower_rational':str(Q(L,Den)),
  'lambda_lower_decimal_for_display':float(Q(L,Den)),
  'Frobenius_radius_upper_rational':str(Q(rho,Den)),
  'fraction_free_integer_Bareiss_pivots_positive':pivots,
  'certification_scope':'finite 47x47 subspace conditional on correct Arb/Cauchy analytic bounds, not RH'
}
Path('P108_COMPACT_47X47.json').write_text(json.dumps(compact,indent=2)+'\n')
Path('P108_COMPACT_PROVENANCE.json').write_text(json.dumps({'shards':refs,'source_sha256':compact['source_sha256'],
  'verified_pivots':pivots,'exact_lambda_lower':str(Q(L,Den)),
  'github_workflow_arb':'https://github.com/jbarumby/scicalc/actions/runs/42489201995'},indent=2)+'\n')
print('SUCCESS: compact 47x47:',len(entries),'entries;',pivots,'positive integer Bareiss pivots; lambda >=',float(Q(L,Den)),
      '; outer panels',sumouter,'; inner ellipses',suminner)
