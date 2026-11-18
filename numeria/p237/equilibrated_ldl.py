#!/usr/bin/env python3
"""Congruence-scaled Arb LDL replay from the P237 full assembled matrix artifact.

Scaling uses exact positive powers of two chosen from each diagonal midpoint.
This is a fixed positive diagonal congruence, so it preserves definiteness.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json
from pathlib import Path
from flint import arb,ctx

N=377;UPPER=N*(N+1)//2

def dyadic(w):
    if not isinstance(w,list) or len(w)!=2:raise ValueError('bad dyadic wire')
    m,e=map(int,w)
    return (arb(m)*(arb(2)**e))

def ball(x):
    m,e=x.mid().man_exp();r,k=x.rad().upper().man_exp()
    return [[str(m),str(e)],[str(r),str(k)]]

def load_matrix(rows):
    if not isinstance(rows,list) or len(rows)!=UPPER:raise ValueError('bad full upper triangle')
    R=[[arb(0) for _ in range(N)] for _ in range(N)]
    last=(-1,-1)
    for row in rows:
        if not isinstance(row,list) or len(row)!=3:raise ValueError('bad matrix row')
        i,j,w=row
        if type(i)!=int or type(j)!=int or not (0<=i<=j<N) or (i,j)<=last:raise ValueError('bad order/index')
        if not isinstance(w,list) or len(w)!=2:raise ValueError('bad ball')
        c=dyadic(w[0]);r=dyadic(w[1])
        if r.lower()<0:raise ValueError('negative radius')
        z=c+arb(0,r.upper())
        if not z.is_finite():raise ValueError('nonfinite matrix input')
        R[i][j]=R[j][i]=z;last=(i,j)
    return R

def scale_matrix(R):
    exponents=[];S=[]
    for i in range(N):
        d=R[i][i]
        if d.lower()<=0:raise ArithmeticError(f'diagonal not strictly positive at {i}: {ball(d)}')
        m,e=d.mid().man_exp()
        m=int(m);e=int(e)
        if m<=0:raise ArithmeticError('nonpositive midpoint diagonal')
        b=m.bit_length()-1+e
        q=-(b//2)
        exponents.append(q);S.append(arb(2)**q)
    B=[[R[i][j]*S[i]*S[j] for j in range(N)] for i in range(N)]
    return B,exponents

def interval_ldl(A):
    L=[[arb(0) for _ in range(N)] for _ in range(N)]
    D=[];pivots=[]
    for k in range(N):
        p=A[k][k]
        for a in range(k):p-=L[k][a]*L[k][a]*D[a]
        if not p.is_finite() or p.lower()<=0:
            return False,k,pivots,ball(p) if p.is_finite() else None
        D.append(p);pivots.append(ball(p))
        for i in range(k+1,N):
            t=A[i][k]
            for a in range(k):t-=L[i][a]*L[k][a]*D[a]
            L[i][k]=t/p
        if (k+1)%25==0:print('SCALED LDL positive pivots',k+1,'/',N,flush=True)
    return True,N,pivots,None

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--bits',type=int,default=2048)
    a=p.parse_args()
    if a.bits<1152:raise ValueError('precision below prior certified assembly precision')
    ctx.prec=a.bits
    raw=a.input.read_bytes();src=json.loads(gzip.decompress(raw))
    if src.get('schema')!='NUMERIA.ES8P237.FULL_FINITE_INTERVAL_LDL.1' or src.get('full_quadrature_coverage_verified') is not True or src.get('all_142506_interval_sums_shadow_audited') is not True or src.get('global_GL60_error_added_once_per_matrix_entry') is not True:
        raise ValueError('input artifact lacks full-C shadow-gated assembly')
    sectors={}
    for sec in ('even','odd'):
        old=src['sectors'][sec]
        R=load_matrix(old['assembled_full_R_upper_triangle'])
        B,exponents=scale_matrix(R)
        ok,count,pivots,failed=interval_ldl(B)
        sectors[sec]={'source_full_C_sha256':old['source_full_C_sha256'],
          'positive_diagonal_congruence_power_of_two_exponents':exponents,
          'interval_LDL_positive_pivot_count':count,'interval_LDL_complete':ok,
          'first_inconclusive_pivot_index':None if ok else count,
          'first_inconclusive_pivot_ball':failed,'pivots':pivots}
        print('SCALED RESULT',sec,'complete',ok,'positive pivots',count,flush=True)
    report={'schema':'NUMERIA.ES8P237.POWER2_CONGRUENCE_LDL.1','precision_bits':a.bits,
      'input_artifact_sha256':hashlib.sha256(raw).hexdigest(),
      'input_shadow_audit_sha256':src['independent_shadow_exact_audit_sha256'],
      'same_full_matrices_and_GL60_beta_pole_enclosures':True,
      'positive_diagonal_congruence_preserves_definiteness':True,
      'both_377_mode_finite_matrices_certified_positive':all(v['interval_LDL_complete'] for v in sectors.values()),
      'infinite_Weil_positivity_proved':False,'RH_proved':False,'sectors':sectors}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_bytes(gzip.compress(json.dumps(report,separators=(',',':')).encode(),compresslevel=9,mtime=0))
    print('CERTIFICATE',report['both_377_mode_finite_matrices_certified_positive'],
          'sha256',hashlib.sha256(a.out.read_bytes()).hexdigest(),flush=True)

if __name__=='__main__':main()
