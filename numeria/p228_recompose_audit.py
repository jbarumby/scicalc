#!/usr/bin/env python3
"""Recompose P227 Arb08 quadrature-sum chunks and compare with P193.

All chunk endpoints, pole bounds, P193 radii and overlap comparisons are exact
Fractions. The sole transcendental interval input generated here is beta_T,
computed with mpmath.iv. P227's quadrature error is independently bounded by
the analytic majorant in p229_quadrature_bound.py. This remains an audit of
finite exported artefacts, not of the global Weil reduction.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
from mpmath import iv
from p229_quadrature_bound import total_error_bound

L=F(4,5)
S_POLE=10**75

def dyad(x):
    m,e=map(int,x)
    return F(m*(1<<e),1) if e>=0 else F(m,1<<(-e))

def iv_endpoints(x):
    a,b=x._mpi_
    def one(v):
        sign,man,exp,bits=v
        z=F((-1 if sign else 1)*man)
        return z*(2**exp) if exp>=0 else z/F(2**(-exp))
    return one(a),one(b)

def beta_interval(T,dps=150):
    iv.dps=dps
    two=iv.mpf(2); ln2=iv.ln(two); ln3=iv.ln(iv.mpf(3)); pi=iv.pi
    A=2*ln2/iv.sqrt(two)+2*ln3/iv.sqrt(iv.mpf(3))+ln2
    b=iv.ln(iv.mpf(T)/(2*pi))-iv.mpf(1)/T-A
    return iv_endpoints(b)

def floor_sqrt_scaled(x: F,scale: int)->int:
    return math.isqrt((x.numerator*scale*scale)//x.denominator)

def ceil_fraction(x: F)->int:
    return -((-x.numerator)//x.denominator)

def pole_component(n:int,terms=32):
    """Exact rational enclosure of sqrt(2L(2n+1))*i_n(L/2), via P214 series."""
    z=L/2
    t=z**n/F(math.prod(range(1,2*n+2,2)))
    partial=F(0)
    for k in range(terms):
        partial+=t
        t=t*z*z/F(2*(k+1)*(2*n+2*k+3))
    ratio=z*z/F(2*(terms+1)*(2*n+2*terms+3))
    lo_i,hi_i=partial,partial+t/(1-ratio)
    rad=2*L*(2*n+1)
    lo_root=F(floor_sqrt_scaled(rad,S_POLE),S_POLE)
    hi_root=F(floor_sqrt_scaled(rad,S_POLE)+1,S_POLE)
    low=lo_root*lo_i; high=hi_root*hi_i
    return F((low*S_POLE).__floor__(),S_POLE),F(ceil_fraction(high*S_POLE),S_POLE)

def read_chunks(chunks_dir:Path,sector:str,T:int):
    paths=sorted(chunks_dir.glob('**/chunk.json.gz'))
    found={}
    for path in paths:
        with gzip.open(path,'rt') as f:d=json.load(f)
        if d.get('sector')!=sector:continue
        if d.get('schema')!='NUMERIA.ES8P227.ARB08_NATIVE_50PANELS.1':raise ValueError(f'bad schema {path}')
        if d.get('version')!='0.8.0' or d.get('backend')!='python-flint/Arb 0.8.0' or d.get('bits',0)<1024:raise ValueError(f'bad engine/precision {path}')
        if d.get('quadrature_error_added') is not False or d.get('beta_and_poles_added') is not False:raise ValueError(f'unexpected assembled chunk {path}')
        a,b=map(int,d['panel_range'])
        if b-a!=50 or a%50 or not (0<=a<b<=2*T) or d.get('Tsharp')!=T or d.get('modes')!=200:raise ValueError(f'bad chunk geometry {path}')
        if a in found:raise ValueError(f'duplicate chunk start {sector} {a}')
        if len(d.get('entries',[]))!=20100:raise ValueError(f'entry count mismatch {path}')
        found[a]=(b,path,d)
    expected=list(range(0,2*T,50))
    if sorted(found)!=expected:raise ValueError(f'{sector} chunks missing/extra: {sorted(found)}')
    # Each cell accumulates exact interval sums from all 50-panel chunks.
    sums={(i,j):(F(0),F(0)) for i in range(200) for j in range(i,200)}
    digest=[]
    for a in expected:
        b,path,d=found[a];digest.append({'file':path.name,'artifact_directory':path.parent.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'panels':[a,b]})
        seen=set()
        for i,j,ball in d['entries']:
            i=int(i);j=int(j)
            if not (0<=i<=j<200) or (i,j) in seen:raise ValueError(f'bad/duplicate matrix index in {path}')
            seen.add((i,j))
            (mm,me),(rm,re)=ball
            center=dyad([mm,me]);radius=dyad([rm,re])
            if radius<0:raise ValueError('negative Arb radius')
            old=sums[(i,j)];sums[(i,j)]=(old[0]+center-radius,old[1]+center+radius)
        if len(seen)!=20100:raise ValueError(f'incomplete triangle in {path}')
    return sums,digest

def read_p193(p193_dir:Path,sector:str):
    cert=json.loads((p193_dir/f'P193_CERT_{sector}.json').read_text())
    gz=p193_dir/f'P193_matrix_{sector}.json.gz'
    if hashlib.sha256(gz.read_bytes()).hexdigest()!=cert['matrix_gz_sha256']:raise ValueError('P193 matrix SHA mismatch')
    raw=gzip.open(gz,'rb').read()
    if hashlib.sha256(raw).hexdigest()!=cert['matrix_content_sha256']:raise ValueError('P193 decompressed-content SHA mismatch')
    m=json.loads(raw)
    T=200 if sector=='even' else 150
    if m.get('N')!=200 or m.get('T')!=T or m.get('sector')!=sector:raise ValueError('P193 metadata mismatch')
    E=F(cert['E_GL_total']['exact'])
    ref={};centers={}
    for i in range(200):
        for j in range(i,200):
            mid=dyad(m['mid'][i][j]);rad=dyad(m['rad_upper'][i][j])
            centers[(i,j)]=mid
            ref[(i,j)]=(mid-rad-E,mid+rad+E)
    return ref,E,cert,centers

def audit_sector(chunks_dir,p193_dir,sector):
    T=200 if sector=='even' else 150
    sums,manifest=read_chunks(chunks_dir,sector,T)
    ref,E193,cert,centers=read_p193(p193_dir,sector)
    E227=total_error_bound(T)
    beta=beta_interval(T)
    parity=0 if sector=='even' else 1
    degrees=[2*k+parity for k in range(200)]
    poles=[pole_component(n) for n in degrees]
    sign=1 if sector=='even' else -1
    overlap=0;disjoint=[];maxmid=F(0);maxwidth=F(0);maxkey=None
    perturbation_rows=[F(0) for _ in range(200)]
    for i in range(200):
        for j in range(i,200):
            lo,hi=sums[(i,j)]
            if i==j:lo+=beta[0];hi+=beta[1]
            pl,ph=poles[i],poles[j]
            if sign>0:p_lo,p_hi=2*pl[0]*ph[0],2*pl[1]*ph[1]
            else:p_lo,p_hi=-2*pl[1]*ph[1],-2*pl[0]*ph[0]
            lo+=p_lo-E227;hi+=p_hi+E227
            center=centers[(i,j)]
            deviation=max(abs(lo-center),abs(hi-center))
            perturbation_rows[i]+=deviation
            if i!=j:perturbation_rows[j]+=deviation
            rlo,rhi=ref[(i,j)]
            if max(lo,rlo)<=min(hi,rhi):overlap+=1
            else:
                disjoint.append([i,j,str(lo),str(hi),str(rlo),str(rhi)])
            delta=abs((lo+hi-rlo-rhi)/2)
            if delta>maxmid:maxmid=delta;maxkey=[i,j]
            width=hi-lo
            if width>maxwidth:maxwidth=width
    row_bound=max(perturbation_rows)
    return {'sector':sector,'Tsharp':T,'dimension':200,'upper_triangle_entries':20100,'entries_intersecting_P193':overlap,'entries_disjoint':len(disjoint),'max_midpoint_difference_exact':str(maxmid),'max_midpoint_difference_decimal':f'{float(maxmid):.8e}','max_difference_entry_ij':maxkey,'max_reconstructed_width_exact':str(maxwidth),'P193_quadrature_error_exact':str(E193),'P227_independent_quadrature_error_exact':str(E227),'P227_to_P193center_infinity_norm_upper_exact':str(row_bound),'P227_to_P193center_infinity_norm_upper_decimal':f'{float(row_bound):.8e}','beta_interval_exact':[str(x) for x in beta],'pole_series_terms':32,'chunks':manifest,'P193_source_sha256':cert['matrix_gz_sha256'],'condition':'P227 Arb balls plus the independent M=100000 Bernstein/Gauss bound in P229; this exact row-sum bound is combined with an independent center SPD certificate in P230.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--chunks',type=Path,required=True);p.add_argument('--p193',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();result={'stage':'NUMERIA ES8-P228 local rational audit of P227','sectors':[audit_sector(a.chunks,a.p193,s) for s in ('even','odd')],'global_claim':'No global tail/coupling verification and no RH proof.'}
    a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
